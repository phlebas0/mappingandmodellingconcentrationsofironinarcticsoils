# This code fits and validate the final Random Forest models

library(randomForest)
library(CAST)
library(sf)
library(terra)
library(spdep)

# Set paths
input_csv <- file.path("data", "full_domain_model_data.csv")
predictor_dir <- "predictors"
boundary_path <- file.path("support", "study_region.shp")
output_dir <- file.path("results", "full_domain_final")
model_dir <- file.path(output_dir, "models")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# Model settings
set.seed(42)

domain_label <- "Full domain"
ntree_cv <- 500
ntree_final <- 1001
nodesize <- 4
random_folds <- 10
spatial_folds <- 5
block_size <- 50000
knndm_samplesize <- 1000
knndm_maxp <- 0.8

# Read modelling data
dat <- read.csv(input_csv, check.names = FALSE)
dat$row_id <- seq_len(nrow(dat))

required_base_fields <- c("Fe_pct", "Fe_react_pct", "Easting", "Northing", "LITH_ORIG")
missing_base_fields <- setdiff(required_base_fields, names(dat))

# Remove the dense sampling block
dat$cluster_block_id <- paste(
  floor(dat$Easting / block_size),
  floor(dat$Northing / block_size),
  sep = "_"
)

cluster_counts <- sort(table(dat$cluster_block_id), decreasing = TRUE)
dense_cluster_block <- names(cluster_counts)[1]

n_before_cluster_exclusion <- nrow(dat)
dat <- dat[dat$cluster_block_id != dense_cluster_block, , drop = FALSE]
n_cluster_removed <- n_before_cluster_exclusion - nrow(dat)
dat$cluster_block_id <- NULL

write.csv(
  data.frame(
    domain = domain_label,
    dense_cluster_block = dense_cluster_block,
    n_before_exclusion = n_before_cluster_exclusion,
    n_removed = n_cluster_removed,
    n_after_exclusion = nrow(dat)
  ),
  file.path(output_dir, "dense_cluster_exclusion_summary.csv"),
  row.names = FALSE
)

cat("Removed", n_cluster_removed, "observations from dense block", dense_cluster_block, "\n")

# Predictor sets
dat$LITH_ORIG <- factor(as.character(dat$LITH_ORIG))
lith_levels <- levels(dat$LITH_ORIG)

non_lithology_predictors <- c(
  "MAG_RESID",
  "ELEV_M", "SLOPE_DEG", "RUGGEDNESS", "TWI",
  "ALT",
  "MAP", "SNOWFALL", "MAT", "TDD", "MAX_SWE",
  "AET",
  "GST", "T1M", "T2M", "T5M", "T10M",
  "SM_0_7", "SM_7_28", "SM_28_100", "SM_100_289"
)

lithology_distances <- sort(grep("^DIST_LITH_", names(dat), value = TRUE))

required_fields <- unique(c(
  required_base_fields,
  non_lithology_predictors,
  lithology_distances
))

missing_fields <- setdiff(required_fields, names(dat))

model_specs <- list(
  M1_bedrock_lithology = list(
    label = "Bedrock lithology",
    geology = "Categorical bedrock lithology",
    predictors = c("LITH_ORIG", non_lithology_predictors)
  ),
  M2_lithology_distances = list(
    label = "Lithology distances",
    geology = "47 lithology-distance surfaces",
    predictors = c(lithology_distances, non_lithology_predictors)
  )
)

model_specifications <- do.call(rbind, lapply(names(model_specs), function(model_name) {
  spec <- model_specs[[model_name]]

  data.frame(
    domain = domain_label,
    model = model_name,
    model_label = spec$label,
    geology = spec$geology,
    n_predictors = length(spec$predictors),
    predictors = paste(spec$predictors, collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  model_specifications,
  file.path(output_dir, "model_specifications.csv"),
  row.names = FALSE
)

# Raster domain used by kNNDM
elevation_path <- file.path(predictor_dir, "elevation_250m_3978.tif")

model_domain <- rast(elevation_path)
boundary_vector <- vect(boundary_path)

if (!same.crs(model_domain, boundary_vector)) {
  boundary_vector <- project(boundary_vector, crs(model_domain))
}

model_domain <- crop(model_domain, boundary_vector)
model_domain <- mask(model_domain, boundary_vector)

performance_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]

  sse <- sum((observed - predicted)^2)
  sst <- sum((observed - mean(observed))^2)
  r2 <- if (length(observed) > 1 && sst > 0) 1 - sse / sst else NA_real_

  data.frame(
    n = length(observed),
    r2 = r2,
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    bias = mean(predicted - observed)
  )
}

moran_residuals <- function(data, observed, predicted, k = 8) {
  keep <- is.finite(observed) & is.finite(predicted)
  residuals <- observed[keep] - predicted[keep]
  coords <- as.matrix(data[keep, c("Easting", "Northing")])
  k_use <- min(k, nrow(coords) - 1)

  if (k_use < 1) {
    return(data.frame(moran_i = NA_real_, moran_p = NA_real_))
  }

  tryCatch({
    neighbours <- knn2nb(knearneigh(coords, k = k_use))
    weights <- nb2listw(neighbours, style = "W", zero.policy = TRUE)

    test <- moran.test(
      residuals,
      weights,
      randomisation = TRUE,
      zero.policy = TRUE,
      alternative = "greater"
    )

    data.frame(
      moran_i = unname(test$estimate["Moran I statistic"]),
      moran_p = test$p.value
    )
  }, error = function(e) {
    data.frame(moran_i = NA_real_, moran_p = NA_real_)
  })
}

make_mtry <- function(number_predictors) {
  max(1, floor(number_predictors / 3))
}

make_random_indices <- function(n) {
  set.seed(42)
  fold_vector <- sample(rep(seq_len(random_folds), length.out = n))
  folds <- sort(unique(fold_vector))

  list(
    train = lapply(folds, function(fold) which(fold_vector != fold)),
    test = lapply(folds, function(fold) which(fold_vector == fold)),
    fold = fold_vector
  )
}

make_block_indices <- function(data) {
  x_block <- floor(data$Easting / block_size)
  y_block <- floor(data$Northing / block_size)
  block_id <- paste(x_block, y_block, sep = "_")
  block_counts <- sort(table(block_id), decreasing = TRUE)

  fold_sizes <- rep(0, spatial_folds)
  block_fold <- integer(length(block_counts))
  names(block_fold) <- names(block_counts)

  for (block in names(block_counts)) {
    fold <- which.min(fold_sizes)
    block_fold[block] <- fold
    fold_sizes[fold] <- fold_sizes[fold] + block_counts[block]
  }

  fold_vector <- unname(block_fold[block_id])
  folds <- sort(unique(fold_vector))

  list(
    train = lapply(folds, function(fold) which(fold_vector != fold)),
    test = lapply(folds, function(fold) which(fold_vector == fold)),
    fold = fold_vector,
    block_id = block_id
  )
}

make_knndm_indices <- function(data) {
  points <- st_as_sf(
    data,
    coords = c("Easting", "Northing"),
    crs = 3978,
    remove = FALSE
  )

  set.seed(42)

  folds <- knndm(
    tpoints = points,
    modeldomain = model_domain,
    dist_space = "geographical",
    k = spatial_folds,
    maxp = knndm_maxp,
    clustering = "hierarchical",
    linkf = "ward.D2",
    samplesize = knndm_samplesize,
    sampling = "regular"
  )

  fold_vector <- rep(NA_integer_, nrow(data))

  for (fold in seq_along(folds$indx_test)) {
    fold_vector[folds$indx_test[[fold]]] <- fold
  }

  list(
    train = folds$indx_train,
    test = folds$indx_test,
    fold = fold_vector,
    W = if (!is.null(folds$W)) folds$W else NA_real_
  )
}

fit_rf <- function(data, response, predictors, ntree) {
  x <- data[, predictors, drop = FALSE]
  y <- data[[response]]

  if ("LITH_ORIG" %in% predictors) {
    x$LITH_ORIG <- factor(as.character(x$LITH_ORIG), levels = lith_levels)
  }

  randomForest(
    x = x,
    y = y,
    ntree = ntree,
    mtry = make_mtry(length(predictors)),
    nodesize = nodesize,
    importance = TRUE,
    na.action = na.omit
  )
}

run_cv <- function(data, response, predictors, indices, validation_name) {
  predictions <- rep(NA_real_, nrow(data))
  fold_used <- rep(NA_integer_, nrow(data))

  for (fold in seq_along(indices$test)) {
    cat(validation_name, "- fold", fold, "of", length(indices$test), "\n")

    train_index <- indices$train[[fold]]
    test_index <- indices$test[[fold]]

    train <- data[train_index, , drop = FALSE]
    test <- data[test_index, , drop = FALSE]

    if ("LITH_ORIG" %in% predictors) {
      test$LITH_ORIG <- factor(as.character(test$LITH_ORIG), levels = lith_levels)
    }

    model <- fit_rf(train, response, predictors, ntree_cv)

    prediction_object <- predict(
      model,
      newdata = test[, predictors, drop = FALSE],
      predict.all = TRUE
    )

    predictions[test_index] <- as.numeric(prediction_object$aggregate)
    fold_used[test_index] <- fold
  }

  metrics <- performance_metrics(data[[response]], predictions)
  moran <- moran_residuals(data, data[[response]], predictions)

  prediction_table <- data.frame(
    domain = domain_label,
    row_id = data$row_id,
    response = response,
    observed = data[[response]],
    predicted = predictions,
    residual = data[[response]] - predictions,
    fold = fold_used,
    validation = validation_name,
    Easting = data$Easting,
    Northing = data$Northing,
    stringsAsFactors = FALSE
  )

  list(
    metrics = cbind(metrics, moran),
    predictions = prediction_table
  )
}

predictor_group <- function(predictor) {
  if (predictor == "LITH_ORIG") return("Bedrock lithology")
  if (grepl("^DIST_LITH_", predictor)) return("Lithology distances")
  if (predictor == "MAG_RESID") return("Magnetics")
  if (predictor %in% c("ELEV_M", "SLOPE_DEG", "RUGGEDNESS", "TWI")) return("Terrain")
  if (predictor %in% c("MAP", "SNOWFALL", "MAT", "TDD", "MAX_SWE")) return("Climate")

  if (predictor %in% c("ALT", "GST", "T1M", "T2M", "T5M", "T10M")) {
    return("Permafrost / ground temperature")
  }

  if (predictor %in% c("SM_0_7", "SM_7_28", "SM_28_100", "SM_100_289")) {
    return("Soil moisture")
  }

  if (predictor == "AET") return("Evapotranspiration")
  "Other"
}

all_results <- list()
all_predictions <- list()
all_importance <- list()

result_counter <- 1
prediction_counter <- 1
importance_counter <- 1

for (response in c("Fe_pct", "Fe_react_pct")) {
  response_label <- if (response == "Fe_pct") "Total Fe" else "Reactive Fe"
  cat("\nRunning", response_label, "models\n")

  # M1 and M2 use the same complete-case observations
  all_model_predictors <- unique(unlist(
    lapply(model_specs, function(spec) spec$predictors),
    use.names = FALSE
  ))

  required <- unique(c(response, "Easting", "Northing", all_model_predictors))
  model_data <- dat[complete.cases(dat[, required, drop = FALSE]), , drop = FALSE]

  model_data$LITH_ORIG <- factor(
    as.character(model_data$LITH_ORIG),
    levels = lith_levels
  )

  cat("Complete-case samples:", nrow(model_data), "\n")

  if (nrow(model_data) < 100) {
    stop("Too few complete cases for ", response_label, ".")
  }

  # The same fold assignments are used for M1 and M2
  random_indices <- make_random_indices(nrow(model_data))
  knndm_indices <- make_knndm_indices(model_data)
  block_indices <- make_block_indices(model_data)

  fold_table <- data.frame(
    domain = domain_label,
    row_id = model_data$row_id,
    Easting = model_data$Easting,
    Northing = model_data$Northing,
    random_fold = random_indices$fold,
    knndm_fold = knndm_indices$fold,
    block_50km_fold = block_indices$fold,
    block_50km_id = block_indices$block_id
  )

  write.csv(
    fold_table,
    file.path(output_dir, paste0(response, "_fold_assignments.csv")),
    row.names = FALSE
  )

  for (model_name in names(model_specs)) {
    spec <- model_specs[[model_name]]
    predictors <- spec$predictors
    p <- length(predictors)
    mtry_value <- make_mtry(p)

    cat("\n", response_label, "-", spec$label, "\n")

    random_result <- run_cv(
      model_data, response, predictors, random_indices, "Random 10-fold CV"
    )

    knndm_result <- run_cv(
      model_data, response, predictors, knndm_indices, "kNNDM CV"
    )

    block_result <- run_cv(
      model_data, response, predictors, block_indices, "50 km block CV"
    )

    # Fit the final model to all complete observations after validation
    final_model <- fit_rf(model_data, response, predictors, ntree_final)

    oob_metrics <- performance_metrics(final_model$y, final_model$predicted)
    oob_moran <- moran_residuals(
      model_data,
      final_model$y,
      final_model$predicted
    )

    saveRDS(
      final_model,
      file.path(model_dir, paste0(response, "_", model_name, "_rf.rds"))
    )

    permutation_importance <- importance(final_model, type = 1, scale = FALSE)

    importance_table <- data.frame(
      domain = domain_label,
      response = response_label,
      model = model_name,
      model_label = spec$label,
      predictor = rownames(permutation_importance),
      permutation_importance = as.numeric(permutation_importance[, 1]),
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    importance_table$predictor_group <- vapply(
      importance_table$predictor,
      predictor_group,
      FUN.VALUE = character(1)
    )

    all_importance[[importance_counter]] <- importance_table
    importance_counter <- importance_counter + 1

    add_result <- function(validation, metrics, knndm_w = NA_real_) {
      data.frame(
        domain = domain_label,
        response = response_label,
        response_field = response,
        model = model_name,
        model_label = spec$label,
        geology = spec$geology,
        n = metrics$n,
        n_predictors = p,
        mtry = mtry_value,
        nodesize = nodesize,
        ntree_cv = ntree_cv,
        ntree_final = ntree_final,
        validation = validation,
        r2 = metrics$r2,
        rmse = metrics$rmse,
        mae = metrics$mae,
        bias = metrics$bias,
        moran_i = metrics$moran_i,
        moran_p = metrics$moran_p,
        knndm_W = knndm_w,
        stringsAsFactors = FALSE
      )
    }

    all_results[[result_counter]] <- add_result(
      "Random 10-fold CV",
      random_result$metrics
    )
    result_counter <- result_counter + 1

    all_results[[result_counter]] <- add_result(
      "kNNDM CV",
      knndm_result$metrics,
      knndm_indices$W
    )
    result_counter <- result_counter + 1

    all_results[[result_counter]] <- add_result(
      "50 km block CV",
      block_result$metrics
    )
    result_counter <- result_counter + 1

    all_results[[result_counter]] <- add_result(
      "OOB",
      cbind(oob_metrics, oob_moran)
    )
    result_counter <- result_counter + 1

    for (prediction_result in list(
      random_result$predictions,
      knndm_result$predictions,
      block_result$predictions
    )) {
      prediction_result$response_label <- response_label
      prediction_result$model <- model_name
      prediction_result$model_label <- spec$label

      all_predictions[[prediction_counter]] <- prediction_result
      prediction_counter <- prediction_counter + 1
    }
  }
}

validation_results <- do.call(rbind, all_results)
oof_predictions <- do.call(rbind, all_predictions)
importance_results <- do.call(rbind, all_importance)

rownames(validation_results) <- NULL
rownames(oof_predictions) <- NULL
rownames(importance_results) <- NULL

# Summarise permutation importance within broad predictor groups
group_keys <- unique(importance_results[, c(
  "domain", "response", "model", "model_label", "predictor_group"
)])

grouped_importance_rows <- vector("list", nrow(group_keys))

for (i in seq_len(nrow(group_keys))) {
  key <- group_keys[i, , drop = FALSE]

  rows <- importance_results$domain == key$domain &
    importance_results$response == key$response &
    importance_results$model == key$model &
    importance_results$predictor_group == key$predictor_group

  values <- importance_results$permutation_importance[rows]

  grouped_importance_rows[[i]] <- data.frame(
    domain = key$domain,
    response = key$response,
    model = key$model,
    model_label = key$model_label,
    predictor_group = key$predictor_group,
    n_predictors = length(values),
    sum_permutation_importance = sum(values, na.rm = TRUE),
    mean_permutation_importance = mean(values, na.rm = TRUE),
    median_permutation_importance = median(values, na.rm = TRUE),
    max_permutation_importance = max(values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

grouped_importance <- do.call(rbind, grouped_importance_rows)
rownames(grouped_importance) <- NULL

write.csv(
  validation_results,
  file.path(output_dir, "validation_results.csv"),
  row.names = FALSE
)

write.csv(
  oof_predictions,
  file.path(output_dir, "oof_predictions.csv"),
  row.names = FALSE
)

write.csv(
  importance_results,
  file.path(output_dir, "variable_importance.csv"),
  row.names = FALSE
)

write.csv(
  grouped_importance,
  file.path(output_dir, "grouped_variable_importance.csv"),
  row.names = FALSE
)

cat("\Outputs written to", output_dir, "\n")
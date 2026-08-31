# This code generates prediction rasters and maps from the fitted Random Forest models
# Predictions are made at 250 m and aggregated to 1 km for display

library(randomForest)
library(terra)
library(sf)
library(ggplot2)

# Set paths
predictor_dir <- "predictors"
support_dir <- "support"
results_dir <- "results"
manifest_path <- "predictor_manifest.csv"

regions <- list(
  full_domain = list(
    label = "Full domain",
    boundary = file.path(support_dir, "study_region.shp"),
    model_dir = file.path(results_dir, "full_domain_final", "models"),
    output_dir = file.path(results_dir, "full_domain_final", "maps_1km")
  ),
  reduced_domain = list(
    label = "Reduced domain",
    boundary = file.path(support_dir, "study_area_2.shp"),
    model_dir = file.path(results_dir, "reduced_domain_final", "models"),
    output_dir = file.path(results_dir, "reduced_domain_final", "maps_1km")
  )
)

# Mapping settings
aggregation_factor <- 4
figure_width <- 8.5
figure_height <- 7
figure_dpi <- 400

# Predictor manifest
manifest <- read.csv(
  manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

manifest$raster_path <- file.path(
  predictor_dir,
  basename(manifest$filename)
)

get_predictor_path <- function(predictor) {
  manifest$raster_path[
    match(predictor, manifest$predictor)
  ]
}

get_rf_predictors <- function(model) {
  names(model$forest$ncat)
}

get_rf_factor_levels <- function(model, variable) {
  predictor_index <- match(
    variable,
    names(model$forest$ncat)
  )

  if (is.na(predictor_index)) return(NULL)

  if (!is.null(model$forest$xlevels)) {
    xlevels <- model$forest$xlevels

    if (!is.null(names(xlevels)) &&
        variable %in% names(xlevels)) {
      return(as.character(xlevels[[variable]]))
    }

    if (length(xlevels) >= predictor_index) {
      candidate <- xlevels[[predictor_index]]

      if (length(candidate) > 0) {
        return(as.character(candidate))
      }
    }
  }

  if (!is.null(model$xlevels) &&
      !is.null(names(model$xlevels)) &&
      variable %in% names(model$xlevels)) {
    return(as.character(model$xlevels[[variable]]))
  }

  NULL
}

# Rebuild factor predictors before raster prediction
predict_rf_safe <- function(model, data) {
  output <- rep(NA_real_, nrow(data))
  factor_predictors <- names(model$forest$ncat)[
    model$forest$ncat > 1
  ]

  for (variable in factor_predictors) {
    if (!variable %in% names(data)) next

    factor_levels <- get_rf_factor_levels(
      model,
      variable
    )

    data[[variable]] <- factor(
      as.character(data[[variable]]),
      levels = factor_levels
    )
  }

  valid_rows <- complete.cases(data)

  if (any(valid_rows)) {
    output[valid_rows] <- predict(
      model,
      newdata = data[
        valid_rows,
        ,
        drop = FALSE
      ]
    )
  }

  output
}

get_response_name <- function(model_file) {
  file_name <- basename(model_file)

  if (grepl("^Fe_react_pct_", file_name)) {
    return("Fe_react_pct")
  }

  if (grepl("^Fe_pct_", file_name)) {
    return("Fe_pct")
  }

  NA_character_
}

get_model_name <- function(model_file, response) {
  file_name <- sub(
    "_rf\\.rds$",
    "",
    basename(model_file)
  )

  sub(
    paste0("^", response, "_"),
    "",
    file_name
  )
}

get_response_label <- function(response) {
  if (response == "Fe_pct") return("Total Fe")
  if (response == "Fe_react_pct") return("Reactive Fe")

  response
}

prediction_records <- list()
record_counter <- 1

for (region_name in names(regions)) {
  region <- regions[[region_name]]

  if (!file.exists(region$boundary) ||
      !dir.exists(region$model_dir)) {
    next
  }

  dir.create(
    region$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  boundary_sf <- st_read(
    region$boundary,
    quiet = TRUE
  )

  boundary_sf <- st_transform(
    boundary_sf,
    3978
  )

  boundary_vect <- vect(boundary_sf)

  model_files <- list.files(
    region$model_dir,
    pattern = "_rf\\.rds$",
    full.names = TRUE
  )

  if (length(model_files) == 0) next

  for (model_file in model_files) {
    model <- readRDS(model_file)

    response <- get_response_name(model_file)
    model_name <- get_model_name(
      model_file,
      response
    )

    predictors <- get_rf_predictors(model)

    predictor_rasters <- vector(
      "list",
      length(predictors)
    )

    for (i in seq_along(predictors)) {
      predictor <- predictors[i]

      r <- rast(
        get_predictor_path(predictor)
      )

      r <- crop(
        r,
        boundary_vect
      )

      r <- mask(
        r,
        boundary_vect
      )

      names(r) <- predictor
      predictor_rasters[[i]] <- r
    }

    predictor_stack <- rast(
      predictor_rasters
    )

    names(predictor_stack) <- predictors

    prediction_250m <- terra::predict(
      predictor_stack,
      model,
      fun = predict_rf_safe,
      na.rm = FALSE,
      cores = 1
    )

    names(prediction_250m) <- response

    prediction_250m <- mask(
      prediction_250m,
      boundary_vect
    )

    tif_250m <- file.path(
      region$output_dir,
      paste0(
        response,
        "_",
        model_name,
        "_prediction_250m.tif"
      )
    )

    writeRaster(
      prediction_250m,
      tif_250m,
      overwrite = TRUE,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES"
      )
    )

    # Aggregate after prediction so the model is still applied at 250 m
    prediction_1km <- aggregate(
      prediction_250m,
      fact = aggregation_factor,
      fun = mean,
      na.rm = TRUE
    )

    names(prediction_1km) <- response

    tif_1km <- file.path(
      region$output_dir,
      paste0(
        response,
        "_",
        model_name,
        "_prediction_1km.tif"
      )
    )

    writeRaster(
      prediction_1km,
      tif_1km,
      overwrite = TRUE,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES"
      )
    )

    prediction_records[[record_counter]] <- data.frame(
      region = region_name,
      region_label = region$label,
      response = response,
      model = model_name,
      tif_1km = tif_1km,
      boundary = region$boundary,
      stringsAsFactors = FALSE
    )

    record_counter <- record_counter + 1
  }
}

prediction_records <- do.call(
  rbind,
  prediction_records
)

rownames(prediction_records) <- NULL

# M1 and M2 use the same colour range within each response and domain
comparison_groups <- unique(
  prediction_records[
    ,
    c("region", "response"),
    drop = FALSE
  ]
)

for (group_index in seq_len(nrow(comparison_groups))) {
  region_name <- comparison_groups$region[group_index]
  response <- comparison_groups$response[group_index]

  response_records <- prediction_records[
    prediction_records$region == region_name &
      prediction_records$response == response,
    ,
    drop = FALSE
  ]

  common_min <- Inf
  common_max <- -Inf

  for (tif_path in response_records$tif_1km) {
    r <- rast(tif_path)

    raster_range <- global(
      r,
      fun = c("min", "max"),
      na.rm = TRUE
    )

    min_i <- raster_range[1, "min"]
    max_i <- raster_range[1, "max"]

    if (is.finite(min_i)) {
      common_min <- min(
        common_min,
        min_i
      )
    }

    if (is.finite(max_i)) {
      common_max <- max(
        common_max,
        max_i
      )
    }
  }

  boundary_sf <- st_read(
    response_records$boundary[1],
    quiet = TRUE
  )

  boundary_sf <- st_transform(
    boundary_sf,
    3978
  )

  for (i in seq_len(nrow(response_records))) {
    map_row <- response_records[
      i,
      ,
      drop = FALSE
    ]

    prediction <- rast(
      map_row$tif_1km
    )

    map_df <- as.data.frame(
      prediction,
      xy = TRUE,
      na.rm = FALSE
    )

    names(map_df)[3] <- "prediction"

    map_title <- paste0(
      get_response_label(response),
      " - ",
      map_row$model
    )

    p <- ggplot() +
      geom_raster(
        data = map_df,
        aes(
          x = x,
          y = y,
          fill = prediction
        )
      ) +
      geom_sf(
        data = boundary_sf,
        fill = NA,
        colour = "black",
        linewidth = 0.45
      ) +
      scale_fill_viridis_c(
        name = "Fe (%)",
        limits = c(
          common_min,
          common_max
        ),
        option = "C",
        na.value = "white"
      ) +
      coord_sf(
        crs = st_crs(3978),
        expand = FALSE
      ) +
      labs(
        title = map_title,
        subtitle = paste0(
          map_row$region_label,
          " | 250 m predictions aggregated to 1 km for display"
        ),
        x = NULL,
        y = NULL
      ) +
      theme_minimal(
        base_size = 11
      ) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_rect(
          fill = "white",
          colour = NA
        ),
        plot.background = element_rect(
          fill = "white",
          colour = NA
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        plot.subtitle = element_text(
          size = 9
        )
      )

    png_path <- file.path(
      dirname(map_row$tif_1km),
      paste0(
        response,
        "_",
        map_row$model,
        "_prediction_1km.png"
      )
    )

    ggsave(
      filename = png_path,
      plot = p,
      width = figure_width,
      height = figure_height,
      units = "in",
      dpi = figure_dpi,
      bg = "white"
    )
  }
}
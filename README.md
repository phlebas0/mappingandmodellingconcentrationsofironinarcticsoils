## Mapping Total and Reactive Iron in Glacial Till across the Slave Craton

### Overview

This repository contains the code and supporting data developed for my Master's in Data Science (Earth and Environment) thesis at Durham University.

The study integrated regional till geochemistry with geological, geophysical, topographic, climatic and permafrost-related datasets to:

1. Predict and map total and reactive iron concentrations across the study region.
2. Compare alternative representations of bedrock geology as predictors of iron.
3. Evaluate predictive performance and spatial generalisability under different validation strategies.

### Methodology

#### Data and feature engineering

Geochemical observations were integrated with environmental predictors, including bedrock geology, aeromagnetic data, terrain derivatives, climate, permafrost characteristics and soil conditions. Spatial datasets were processed and harmonised to a common 250 m prediction grid.

#### Random Forest modelling

Separate Random Forest regression models were developed for total and reactive Fe concentrations.

Two approaches to representing bedrock geology were compared:

- **M₁ - Categorical lithology:** Bedrock lithology classification alongside 21 environmental and geophysical predictors (22 predictors in total).
- **M₂ - Lithological proximity:** Euclidean distances to 47 bedrock lithologies alongside the same 21 predictors (68 predictors in total).

Final models used the R `randomForest` package with the following parameters:

- `ntree = 1001`
- `nodesize = 4`
- `mtry = floor(p / 3)`, where `p` is the number of predictors.

#### Model validation

Predictive performance and spatial generalisability were evaluated using:

- out-of-bag (OOB) validation,
- random 10-fold cross-validation,
- kNNDM spatial cross-validation (5-fold),
- spatial block cross-validation (5-fold; 50 km blocks), and
- residual spatial autocorrelation analysis using Moran's I.

Performance was assessed using R², RMSE, MAE and prediction bias. Comparing conventional and spatially-structured validation provided insight into the effects of spatial dependence, sampling heterogeneity, and geographical separation on predictive performance.

### Findings

The analysis identified several important considerations for spatial prediction:

- **Validation sensitivity:** Random cross-validation consistently produced higher performance estimates than spatially structured validation, suggesting optimistic estimates of generalisability to geographically separated observations. This is consistent with a growing body of literature highlighting the limitations of conventional cross-validation for spatially dependent environmental data (Roberts et al., 2017; Meyer & Pebesma, 2022).
- **Geological representation:** Lithological proximity predictors modestly improved conventional validation performance but did not consistently improve spatial generalisability. One possible explanation is that glacial transport and depositional processes weaken relationships between till geochemistry and geographically proximate bedrock (McClenaghan et al., 2023). More broadly, Cracknell & Reading (2014) demonstrated that the predictive contribution of spatial information in geological machine learning is sensitive to the distribution of training observations.
- **Response variability:** Reactive Fe generally exhibited stronger performance under conventional validation, although both responses remained challenging to predict at geographically separated locations. This may reflect differences in their underlying environmental controls: the reactive Fe proxy may capture spatial patterns associated with weathering and environmental conditions, whereas total Fe concentrations may be more strongly influenced by heterogeneous parent-material composition and glacial sediment provenance.
- **Spatial dependence:** Model performance was sensitive to the geographical distribution of observations, highlighting the importance of validation procedures that reflect the intended prediction task. However, spatial cross-validation does not necessarily provide an 'improved' or unbiased estimate of map-wide accuracy, particularly where observations are unevenly distributed or environmentally unrepresentative (e.g., Wadoux et al., 2021; Meyer & Pebesma, 2022).

### Repository contents

- `fit_random_forest_models.R` - Model training, evaluation and cross-validation.
- `predict_and_map_models.R` - Spatial prediction and mapping.
- `predictor_rasters/` - Environmental and geological predictor datasets.
- `response_var/` - Geochemical response data.

### References
- Cracknell, M. J. & Reading, A. M. (2014). Geological mapping using remote sensing data: A comparison of five machine learning algorithms, their response to variations in the spatial distribution of training data and the use of explicit spatial information. Computers & Geosciences, 63, 22–33. https://doi.org/10.1016/j.cageo.2013.10.008
- McClenaghan, M. B. et al. (2023). Review of till geochemistry and indicator mineral methods for mineral exploration in glaciated terrain. Geochemistry: Exploration, Environment, Analysis, 23(4), geochem2023-013. https://doi.org/10.1144/geochem2023-013
- Meyer, H. & Pebesma, E. (2022). Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications, 13, 2208. https://doi.org/10.1038/s41467-022-29838-9
- Roberts, D. R. et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography, 40(8), 913–929. https://doi.org/10.1111/ecog.02881
- Wadoux, A. M. J.-C., Heuvelink, G. B. M., de Bruin, S. & Brus, D. J. (2021). Spatial cross-validation is not the right way to evaluate map accuracy. Ecological Modelling, 457, 109692. https://doi.org/10.1016/j.ecolmodel.2021.109692

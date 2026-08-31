## Mapping Total and Reactive Iron in Glacial Till across the Slave Craton
Spatial Random Forest modelling of total and reactive iron in glacial till across the Slave Craton, northern Canada.

### Overview
This repository contains the code used for an MSc Environmental Data Science dissertation at Durham University.

The project uses **Random Forest** models to predict total and reactive iron concentrations in glacial till across the **Slave Craton, northern Canada**. Models combine till geochemistry with geological, geophysical, terrain, climate, permafrost and soil-moisture predictors.

### Models
Two representations of bedrock geology are compared:
- **M1 — Bedrock lithology:** categorical bedrock lithology + 21 environmental and geophysical predictors
- **M2 — Lithology distances:** 47 distance-to-lithology predictors + the same 21 additional predictors

Separate models are fitted for:
- **Total Fe** (`Fe_pct`)
- **Reactive Fe** (`Fe_react_pct`)

#### Validation
Model performance is evaluated using:
- Out-of-bag validation,
- random 10-fold cross-validation,
- 5-fold kNNDM spatial cross-validation,
- 5-fold 50 km block cross-validation, and
- residual Moran's I.

Final models use the `randomForest` package with:
```r
ntree = 1001
nodesize = 4
mtry = floor(p / 3)

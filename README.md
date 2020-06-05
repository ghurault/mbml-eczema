# Personalised prediction of daily eczema severity scores using a mechanistic machine learning model

This repository contains the code written for the article by **Hurault et al., "Personalised prediction of daily eczema severity scores  using a mechanistic machine learning model"**.
The code is written in the R language for statistical computing and the models using the probabilistic programming language [Stan](https://mc-stan.org/).

The `Models` folder contains:

- Our two models
  - [`BaseModel.stan`](Models/BaseModel.stan): the auto-regressive model based on an Exponentially Modified Gaussian distribution.
  - [`ExtendedModel.stan`](Models/ExtendedModel.stan): an extension of this model to include measurements present in SWET but not in the Flares dataset.
This model notably takes into account the quantity and potency of the treatment use in the determination of the different treatment responsiveness, as well as demographic factors.
- Two benchmarks (the other benchmarks, the uniform and historical forecasts, are computed in the [`validation.R`](validation.R) script)
  - [`Autoregression.stan`](Models/Autoregression.stan): an auto-regressive model (our model without Flares triggers).
  - [`RandomWalk.stan`](Models/RandomWalk.stan): a Gaussian random walk model.

These models, except the ExtendedModel, can be fitted to the Flares and SWET datasets in [`fitting.R`](fitting.R).
The ExtendedModel can be fitted in [`fitting_ext.R`](fitting_ext.R).

Forward chaining is implemented in [`validation.R`](validation.R) and these results analysed in [`results_validation.R`](results_validation.R).

Utility functions used within the scripts are available in [`functions.R`](functions.R).
In addition, we used functions from Guillem Hurault's personal package, [HuraultMisc](https://github.com/ghurault/HuraultMisc).

The Flares and SWET datasets are not available according to our data sharing agreement.
During the analysis, these datasets are loaded from a proprietary package `TanakaData` which includes the raw files as well as data processing functions.
Further processing is performed by functions written in [`functions_data.R`](functions_data.R).

Nonetheless, it is possible to generate fake data from the prior predictive distribution of the different proposed models.
We implement this, as well as prior predictive check and fake data check in [`prior_fake_check.R`](prior_fake_check.R).

Finally, [`plots.R`](plots.R) is used to produce several plots present in the paper.

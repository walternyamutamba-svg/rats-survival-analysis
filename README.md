# Rats Survival Analysis

This repository contains my project on **parametric survival analysis** using the **rats dataset** in R.

## Project Overview

The aim of this project is to analyze survival times in the rats dataset using standard survival analysis techniques and parametric survival models. The analysis includes:

- data exploration and summary statistics
- Kaplan–Meier survival analysis
- treatment group comparison
- fitting multiple parametric survival models
- model selection using AIC
- life function estimation for the best-fitting model
- interpretation of survival and hazard patterns

## Dataset

The project uses the **rats dataset** from the `survival` package in R.

Key variables include:

- `litter`: litter identifier
- `rx`: treatment group
- `time`: observed survival time
- `status`: event indicator
- `event`: recoded event variable used in the analysis
- `treatment`: treatment label derived from `rx`

## Files in this Repository

- `rats_survival_analysis.Rmd` — R Markdown source file containing the full analysis
- `rats_survival_presentation.pdf` — Beamer presentation slides in PDF format
- `README.md` — project description and repository guide

## Methods Used

The following methods were used in this project:

- Kaplan–Meier estimator
- Log-rank test
- Parametric survival models:
  - Exponential
  - Weibull
  - Log-normal
  - Log-logistic
  - Gamma
  - Gompertz
  - Generalized Gamma
- Model comparison using AIC

## Main Result

Among the fitted parametric models, the **Weibull model** provided the best fit based on the lowest AIC. The analysis suggests that the hazard of tumor occurrence increases over time.

## Software

This project was completed in **R** using packages including:

- `survival`
- `flexsurv`
- `ggplot2`
- `dplyr`
- `survminer`
- `patchwork`
- `knitr`
- `kableExtra`

## Author

Walter M Nyamutamba
Samuel Agyekum
Honorine
Belyse
Grace

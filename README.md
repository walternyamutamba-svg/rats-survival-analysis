# Rats Survival Analysis

Parametric survival analysis of the `rats` dataset using R.

## Live Report

Full rendered analysis available at:
https://rpubs.com/Grace_k/rats-survival-analysis

## Project Overview

This project analyses time-to-tumour onset in rats exposed to a
carcinogenic agent, using a full pipeline of survival analysis methods:

- Data exploration and descriptive summary statistics
- Kaplan–Meier estimation with log-rank test
- Cox proportional hazards models (unadjusted, adjusted, clustered SE)
- Proportional hazards assumption diagnostics (Schoenfeld residuals)
- Weibull Accelerated Failure Time (AFT) model
- 7 parametric survival models with AIC-based model selection
- Life function estimation — S(t), h(t), f(t), F(t) — overall and by group
- Shared gamma frailty model to account for within-litter clustering

## Dataset

The `rats` dataset from the `survival` package in R contains data from
a study where rats were monitored for tumour development after exposure
to a carcinogenic agent.

**Key variables:**

| Variable | Description |
|---|---|
| `time` | Observed survival time (days) |
| `status` / `event` | 1 = tumour occurred, 0 = censored |
| `treatment` | Control vs Treated (derived from `rx`) |
| `sex` | Sex of the rat (female / male) |
| `litter` | Litter identifier — clustering variable |

**Study design:** 2:1 allocation (200 Control, 100 Treated),
max follow-up = 104 days.

## Main Results

| Result | Value |
|---|---|
| Total observations | 300 |
| Events (tumours) | 42 (14%) |
| Log-rank p-value | 0.018 |
| Best parametric model | Weibull (AIC = 578.21) |
| Weibull shape k | 3.667 — increasing hazard confirmed |
| Cox HR (Treated vs Control) | 2.21 (p = 0.011) |
| Sex HR (Male vs Female) | 0.047 (p < 0.001) |
| Mean survival — Control | 173.2 days |
| Mean survival — Treated | 124.8 days |
| Frailty variance θ | 2.02 — substantial litter clustering |
| Frailty HR (Treated) | 2.069 (p = 0.022) |

**Key finding:** The treated group has approximately double the hazard
of tumour onset. Sex is the strongest predictor — male rats are 95%
less likely to develop tumours than females. The Weibull shape parameter
k > 1 confirms an increasing hazard over time.

## Files in this Repository

| File | Description |
|---|---|
| `Group4_SA_1.Rmd` | Full R Markdown source (complete analysis) |
| `Group4_SA_1.R` | Plain R script version |
| `Group4_SA_1.pdf` | Beamer presentation slides |
| `outputs/km_by_treatment.png` | KM curves by treatment group |
| `outputs/km_by_sex.png` | KM curves by sex |
| `outputs/aic_comparison.png` | AIC bar chart — all 7 models |
| `outputs/cox_adjusted_curves.png` | Adjusted Cox survival curves |
| `outputs/schoenfeld_residuals.png` | PH assumption diagnostics |
| `outputs/life_functions_overall.png` | 2x2 life function plots (overall) |
| `outputs/life_functions_by_group.png` | 2x2 life function plots (by group) |
| `outputs/pp_plot.png` | PP-plot goodness of fit |

## Methods Used

- Kaplan–Meier estimator
- Log-rank test
- Cox proportional hazards regression
- Schoenfeld residuals test (PH assumption)
- Weibull AFT model
- Parametric models: Exponential, Weibull, Log-Normal,
  Log-Logistic, Gamma, Gompertz
- Model comparison using AIC 
- Shared gamma frailty model for litter clustering

## Software

R packages used: `survival`, `flexsurv`, `survminer`, `ggplot2`,
`dplyr`, `knitr`, `kableExtra`, `broom`, `gridExtra`, `patchwork`

## Authors

- Walter M Nyamutamba
- Samuel Agyekum
- Honorine Mujawayezu
- Belyse Kaneza
- Grace Kitonyi

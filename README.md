# Rats Survival Analysis

> **Parametric, Non-Parametric, and Semi-Parametric Survival Analysis**  
> African Institute for Mathematical Sciences (AIMS), Rwanda — March 2026

--- 

## Live Reports

| Report | Link |
|---|---|
| Project 1 — Parametric Analysis | [View on RPubs](https://rpubs.com/Grace_k/rats-survival-analysis) |
| Project 2 — Extended Workflow | [View on RPubs](https://rpubs.com/Grace_k/1414829) |
| GitHub Repository | [walternyamutamba-svg/rats-survival-analysis](https://github.com/walternyamutamba-svg/rats-survival-analysis) |

---

## Project Overview

This repository contains a two-part survival analysis project applied to the
`rats` dataset from R's `survival` package. The dataset records tumour
incidence in rats from a carcinogenesis trial (Mantel et al., 1977).

**Project 1** focuses on parametric survival modelling — fitting and comparing
seven distributional families, selecting the best model via AIC, and deriving
the four life functions S(t), h(t), f(t), F(t).

**Project 2** extends the analysis with non-parametric and semi-parametric
methods — Kaplan–Meier estimation, Nelson–Aalen cumulative hazard, formal
hypothesis testing, full Cox regression with diagnostics, model refinement
for PH violations, and a comprehensive synthesis integrating all three
analytical frameworks.

---

## Dataset

The `rats` dataset from the `survival` package in R records data from a
carcinogenesis experiment by Mantel, Bohidar & Ciminera (1977).

| Variable | Type | Description |
|---|---|---|
| `time` | Numeric | Days to tumour onset or censoring |
| `status` / `event` | Binary | 1 = tumour occurred, 0 = censored |
| `treatment` | Factor | Control vs Treated (derived from `rx`) |
| `sex` | Factor | Female / Male |
| `litter` | Factor | Litter identifier — clustering variable (100 litters) |

**Study design:** 2:1 allocation — 200 Control, 100 Treated.
Max follow-up = 104 days.

---

## Main Results

### Project 1 — Parametric Survival Analysis

| Result | Value |
|---|---|
| Total observations | 300 |
| Events (tumours) | 42 (14%) |
| Best parametric model | **Weibull** (AIC = 578.21) |
| Weibull shape k | 3.667 — increasing hazard confirmed |
| Weibull scale λ | 160.6 days |
| Mean survival — Overall | ~150.6 days |
| Mean survival — Control | ~173.2 days |
| Mean survival — Treated | ~124.8 days |
| Difference in mean survival | 48.4 days |
| Frailty variance θ | 2.02 — substantial litter clustering |
| Frailty HR (Treated) | 2.069 (p = 0.022) |

### Project 2 — Extended Workflow

| Result | Value |
|---|---|
| Median survival | Not reached (event rate < 50%) |
| RMST — Control | 100.38 days |
| RMST — Treated | 98.55 days |
| Log-rank p (treatment) | 0.018 — significant |
| Log-rank p (sex) | < 0.001 — highly significant |
| Wilcoxon p (treatment) | ~0.023 — consistent with log-rank |
| Cox HR — Treated vs Control (adjusted) | 2.21 (p = 0.011) |
| Cox HR — Male vs Female (adjusted) | 0.047 (p < 0.001) |
| PH assumption — Treatment | Violated (p = 0.020) |
| PH assumption — Sex | Holds (p = 0.438) |
| Stratified Cox HR (treatment) | ~2.2 — consistent after stratification |
| Time-dependent coefficient | Significant — HR varies with log(t) |

**Key findings:**
- The carcinogenic drug doubles the hazard of tumour onset across all models
- Sex is the dominant predictor — male rats are 95% less likely to develop tumours than females
- The Weibull k > 1 and the Nelson–Aalen smoothed hazard both confirm an increasing hazard — the two estimates align closely, validating the parametric assumption
- The PH violation for treatment motivates parametric and AFT approaches over the standard Cox model
- Litter clustering is substantial (θ = 2.02) — the frailty model is essential for valid inference

---

## Methods Used

### Project 1
- Descriptive summary statistics
- Kaplan–Meier estimation (overall and by group)
- Log-rank test
- Cox proportional hazards models (unadjusted, adjusted, clustered SE)
- Schoenfeld residuals test (PH assumption check)
- Weibull AFT model with time ratios
- 7 parametric models: Exponential, Weibull, Log-Normal, Log-Logistic, Gamma, Gompertz, Generalized Gamma
- Model comparison using AIC and BIC
- Life function estimation — S(t), h(t), f(t), F(t) — overall and by treatment group
- PP-plot goodness of fit
- Shared gamma frailty model for litter clustering

### Project 2
- Kaplan–Meier estimator with stratification (treatment and sex)
- Median survival, RMST, survival probabilities at percentile time points
- Nelson–Aalen cumulative hazard estimator
- Non-parametric smoothed hazard estimation
- Overlay comparison of Nelson–Aalen and Weibull hazard
- Log-rank test and Wilcoxon (Breslow) test
- Univariable and multivariable Cox models
- Forest plot of hazard ratios
- Proportional hazards diagnostics (Schoenfeld residuals, dfbeta, Martingale, deviance)
- Stratified Cox model (stratified by sex)
- Time-dependent coefficient Cox model
- Full synthesis integrating parametric, non-parametric, and semi-parametric findings

---

## Files in this Repository

### Source Code
| File | Description |
|---|---|
| `rats_survival_analysis.Rmd` | Project 1 — Full R Markdown source |
| `rats_survival_analysis.R` | Project 1 — Plain R script |
| `assignment2_survival_analysis.Rmd` | Project 2 — Full R Markdown source |
| `rats_survival_presentation.pdf` | Beamer slides (13 slides) |

### Project 1 Outputs
| File | Description |
|---|---|
| `outputs/km_by_treatment.png` | KM curves by treatment group |
| `outputs/km_by_sex.png` | KM curves by sex |
| `outputs/aic_comparison.png` | AIC bar chart — all 7 parametric models |
| `outputs/cox_adjusted_curves.png` | Adjusted Cox survival curves |
| `outputs/schoenfeld_residuals.png` | Schoenfeld residuals — PH diagnostics |
| `outputs/life_functions_overall.png` | 2x2 life function plots (overall) |
| `outputs/life_functions_by_group.png` | 2x2 life function plots (by group) |
| `outputs/pp_plot.png` | PP-plot goodness of fit |

### Project 2 Outputs
| File | Description |
|---|---|
| `outputs/plots/km_overall.png` | Overall KM survival curve |
| `outputs/plots/km_by_treatment.png` | KM curves by treatment |
| `outputs/plots/km_by_sex.png` | KM curves by sex |
| `outputs/plots/nelson_aalen_overall.png` | Nelson–Aalen cumulative hazard |
| `outputs/plots/nelson_aalen_by_treatment.png` | Cumulative hazard by treatment |
| `outputs/plots/hazard_na_vs_weibull.png` | Smoothed hazard vs Weibull overlay |
| `outputs/plots/forest_plot.png` | Cox model forest plot |
| `outputs/plots/schoenfeld_residuals.png` | Schoenfeld residuals |
| `outputs/plots/dfbeta_residuals.png` | Influential observations |
| `outputs/plots/martingale_deviance.png` | Martingale and deviance residuals |
| `outputs/plots/km_vs_weibull.png` | KM vs Weibull parametric overlay |
| `outputs/tables/survival_at_timepoints.csv` | Survival probabilities table |
| `outputs/tables/hypothesis_tests.csv` | Log-rank and Wilcoxon results |
| `outputs/tables/cox_multivariable.csv` | Multivariable Cox model table |
| `outputs/tables/cox_stratified.csv` | Stratified Cox model table |
| `outputs/tables/cox_time_dependent.csv` | Time-dependent Cox model table |
| `outputs/tables/ph_schoenfeld_test.csv` | PH assumption test results |

---

## Software

**Language:** R

**Packages:** `survival`, `flexsurv`, `survminer`, `ggplot2`, `dplyr`,
`knitr`, `kableExtra`, `broom`, `gridExtra`, `patchwork`

---

## References

- Mantel, N., Bohidar, N. R., & Ciminera, J. L. (1977). Mantel-Haenszel analyses of litter-matched time-to-response data. *Cancer Research*, 37(11), 3863–3868.
- Therneau, T. M., & Grambsch, P. M. (2000). *Modeling Survival Data: Extending the Cox Model*. Springer.
- Jackson, C. (2016). flexsurv: A platform for parametric survival modeling in R. *Journal of Statistical Software*, 70(8), 1–33.

---

## Authors

| Name | Institution |
|---|---|
| Walter M Nyamutamba | AIMS Rwanda |
| Samuel Agyekum | AIMS Rwanda |
| Honorine Mujawayezu | AIMS Rwanda |
| Belyse Kaneza | AIMS Rwanda |
| Grace Kitonyi | AIMS Rwanda |

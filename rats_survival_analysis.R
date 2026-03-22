```{r}
# =============================================================
# Parametric Survival Analysis — Rats Dataset
# Authors: Walter M Nyamutamba, Samuel Agyekum, Honorine Mujawayezu,
#          Belyse Kaneza, Grace Kitonyi
# Date: March 2026
# =============================================================

# --- Packages ------------------------------------------------
library(survival)
library(flexsurv)
library(survminer)
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)
library(broom)
library(gridExtra)
library(patchwork)


# --- Data Preparation ----------------------------------------
data(rats)

rats <- rats %>%
  mutate(
    litter    = as.factor(litter),
    event     = as.numeric(status),
    treatment = factor(ifelse(rx == 1, "Treated", "Control"),
                       levels = c("Control", "Treated")),
    sex       = as.factor(sex)
  )

rats_clean <- rats %>%
  select(litter, time, event, treatment, sex) %>%
  na.omit()

cat("Sample size:", nrow(rats_clean), "\n")
cat("Events:", sum(rats_clean$event), "\n")
cat("Censored:", sum(rats_clean$event == 0), "\n")
cat("Event rate:", round(mean(rats_clean$event) * 100, 1), "%\n")


# --- Missing Values ------------------------------------------
cat("\nMissing values per variable:\n")
print(colSums(is.na(rats)))


# --- Descriptive Summary -------------------------------------
overall <- data.frame(
  Group          = "Overall",
  N              = nrow(rats_clean),
  Events         = sum(rats_clean$event),
  Censored       = sum(rats_clean$event == 0),
  Event_Rate_Pct = round(mean(rats_clean$event) * 100, 1),
  Min_Time       = min(rats_clean$time),
  Max_Time       = max(rats_clean$time),
  Mean_Time      = round(mean(rats_clean$time), 2),
  Median_Time    = median(rats_clean$time),
  SD_Time        = round(sd(rats_clean$time), 2)
)

by_group <- rats_clean %>%
  group_by(Group = treatment) %>%
  summarise(
    N              = n(),
    Events         = sum(event),
    Censored       = sum(event == 0),
    Event_Rate_Pct = round(mean(event) * 100, 1),
    Min_Time       = min(time),
    Max_Time       = max(time),
    Mean_Time      = round(mean(time), 2),
    Median_Time    = median(time),
    SD_Time        = round(sd(time), 2),
    .groups        = "drop"
  ) %>%
  mutate(Group = as.character(Group))

print(bind_rows(overall, by_group))


# --- EDA Plots -----------------------------------------------
p1 <- ggplot(rats_clean, aes(x = time, fill = treatment)) +
  geom_histogram(bins = 20, alpha = 0.75, color = "white",
                 position = "identity") +
  scale_fill_manual(values = c("Control" = "#2E86AB",
                               "Treated" = "#F18F01")) +
  labs(title = "Distribution of Survival Times",
       x = "Time (days)", y = "Count", fill = "") +
  theme_minimal()

p2 <- ggplot(rats_clean, aes(x = treatment, y = time,
                              fill = treatment)) +
  geom_boxplot(alpha = 0.85) +
  scale_fill_manual(values = c("#2E86AB", "#F18F01")) +
  labs(title = "Survival Time by Treatment",
       x = "", y = "Time (days)") +
  theme_minimal() + theme(legend.position = "none")

p3 <- ggplot(rats_clean, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_fill_manual(values = c("#A23B72", "#3B7080")) +
  labs(title = "Sex Distribution", x = "", y = "Count") +
  theme_minimal() + theme(legend.position = "none")

grid.arrange(p1, p2, p3, ncol = 3)


# --- Survival Object -----------------------------------------
surv_obj <- Surv(time = rats_clean$time, event = rats_clean$event)


# --- Kaplan-Meier --------------------------------------------
km_overall   <- survfit(surv_obj ~ 1, data = rats_clean)
km_treatment <- survfit(surv_obj ~ treatment, data = rats_clean)
km_sex       <- survfit(surv_obj ~ sex, data = rats_clean)

ggsurvplot(
  km_treatment, data = rats_clean,
  conf.int = TRUE, risk.table = TRUE,
  pval = TRUE, pval.method = TRUE,
  palette = c("#2E86AB", "#F18F01"),
  legend.labs = c("Control", "Treated"),
  legend.title = "Group",
  ggtheme = theme_minimal(),
  title = "Kaplan-Meier Curves by Treatment",
  subtitle = "Median survival not reached - event rate < 50%",
  xlab = "Time (days)", ylab = "S(t)",
  surv.median.line = "none",
  risk.table.title = "Number at Risk",
  tables.theme = theme_cleantable()
)

# Log-rank test
lr_test <- survdiff(surv_obj ~ treatment, data = rats_clean)
pval    <- 1 - pchisq(lr_test$chisq, df = 1)
cat(sprintf("\nLog-rank test: Chi-sq = %.4f, p = %.4f\n",
            lr_test$chisq, pval))

# RMST
rmst_groups <- summary(km_treatment)$table
cat(sprintf("Control RMST:  %.2f days\n",
            rmst_groups["treatment=Control", "rmean"]))
cat(sprintf("Treated RMST:  %.2f days\n",
            rmst_groups["treatment=Treated", "rmean"]))

# Survival at key time points
key_times <- c(25, 50, 75, 100)
sp        <- summary(km_treatment, times = key_times)
surv_table <- data.frame(
  Time      = sp$time,
  Group     = gsub("treatment=", "", sp$strata),
  Survival  = round(sp$surv,  4),
  Lower_95  = round(sp$lower, 4),
  Upper_95  = round(sp$upper, 4),
  N_at_risk = sp$n.risk
)
print(surv_table)


# --- Cox Models ----------------------------------------------
cox_simple <- coxph(Surv(time, event) ~ treatment,
                    data = rats_clean)
print(summary(cox_simple))

cox_multi <- coxph(Surv(time, event) ~ treatment + sex,
                   data = rats_clean)
print(summary(cox_multi))

cox_cluster <- coxph(
  Surv(time, event) ~ treatment + sex + cluster(litter),
  data = rats_clean
)
print(summary(cox_cluster))

# Tie-handling comparison
cox_b <- coxph(Surv(time, event) ~ treatment + sex,
               data = rats_clean, ties = "breslow")
cox_e <- coxph(Surv(time, event) ~ treatment + sex,
               data = rats_clean, ties = "efron")
cox_x <- coxph(Surv(time, event) ~ treatment + sex,
               data = rats_clean, ties = "exact")

ties_table <- data.frame(
  Method       = c("Breslow", "Efron", "Exact"),
  Treatment_HR = round(exp(c(
    coef(cox_b)["treatmentTreated"],
    coef(cox_e)["treatmentTreated"],
    coef(cox_x)["treatmentTreated"])), 3),
  Treatment_p  = round(c(
    summary(cox_b)$coefficients["treatmentTreated", "Pr(>|z|)"],
    summary(cox_e)$coefficients["treatmentTreated", "Pr(>|z|)"],
    summary(cox_x)$coefficients["treatmentTreated", "Pr(>|z|)"]),
    4)
)
print(ties_table)

# Adjusted survival curves
newdata_cox <- data.frame(
  treatment = factor(c("Control", "Treated"),
                     levels = c("Control", "Treated")),
  sex       = factor(c("f", "f"),
                     levels = levels(rats_clean$sex))
)
cox_surv <- survfit(cox_multi, newdata = newdata_cox)

ggsurvplot(
  cox_surv, data = newdata_cox, conf.int = TRUE,
  palette = c("#2E86AB", "#F18F01"),
  legend.labs = c("Control", "Treated"),
  ggtheme = theme_minimal(),
  title = "Adjusted Survival Curves from Cox Model (Female)",
  xlab = "Time (days)", ylab = "S(t)"
)


# --- PH Assumption Check -------------------------------------
ph_test <- cox.zph(cox_multi)
print(ph_test$table)

par(mfrow = c(1, 2))
plot(ph_test, var = 1,
     main = "Schoenfeld Residuals - Treatment",
     ylab = "Beta(t)", xlab = "Time (days)",
     col = "#F18F01", lwd = 2)
abline(h = coef(cox_multi)[1], col = "red", lty = 2, lwd = 2)
plot(ph_test, var = 2,
     main = "Schoenfeld Residuals - Sex",
     ylab = "Beta(t)", xlab = "Time (days)",
     col = "#A23B72", lwd = 2)
abline(h = coef(cox_multi)[2], col = "red", lty = 2, lwd = 2)
par(mfrow = c(1, 1))


# --- Weibull AFT Model ---------------------------------------
aft_weibull <- survreg(
  Surv(time, event) ~ treatment + sex,
  data = rats_clean, dist = "weibull"
)
print(summary(aft_weibull))

aft_coefs <- coef(aft_weibull)
aft_coefs <- aft_coefs[!names(aft_coefs) %in%
                         c("(Intercept)", "Log(scale)")]
cat("\nTime Ratios from Weibull AFT:\n")
print(round(exp(aft_coefs), 3))


# --- Parametric Model Fitting --------------------------------
dist_list <- c(
  "Exponential"       = "exp",
  "Weibull"           = "weibull",
  "Log-Normal"        = "lnorm",
  "Log-Logistic"      = "llogis",
  "Gamma"             = "gamma",
  "Gompertz"          = "gompertz",
  "Generalized Gamma" = "gengamma"
)

fits <- lapply(dist_list, function(d) {
  tryCatch(
    flexsurvreg(surv_obj ~ 1, data = rats_clean, dist = d),
    error = function(e) NULL
  )
})
names(fits) <- names(dist_list)
fits <- Filter(Negate(is.null), fits)

model_comparison <- data.frame(
  Model          = names(fits),
  Log_Likelihood = sapply(fits, function(m) round(logLik(m)[1], 3)),
  AIC            = sapply(fits, function(m) round(AIC(m), 3)),
  BIC            = sapply(fits, function(m) round(BIC(m), 3)),
  Num_Params     = sapply(fits, function(m) length(m$coefficients))
) %>%
  arrange(AIC)

rownames(model_comparison) <- NULL
print(model_comparison)

best_model_name <- model_comparison$Model[1]
best_model      <- fits[[best_model_name]]
cat(sprintf("\nBest model: %s  (AIC = %.3f)\n",
            best_model_name, model_comparison$AIC[1]))


# --- All Models vs KM ----------------------------------------
colors   <- c("#E53935","#1E88E5","#43A047","#FB8C00",
              "#8E24AA","#00ACC1","#6D4C41")
time_seq <- seq(0.1, max(rats_clean$time), length.out = 300)

plot(km_overall, col = "black", lwd = 2.5, conf.int = FALSE,
     xlab = "Time (days)", ylab = "S(t)",
     main = "All Parametric Models vs Kaplan-Meier",
     ylim = c(0, 1))
for (i in seq_along(fits)) {
  sp <- summary(fits[[i]], t = time_seq, type = "survival")[[1]]
  lines(sp$time, sp$est, col = colors[i], lwd = 1.8, lty = 2)
}
legend("bottomleft",
       legend = c("Kaplan-Meier", names(fits)),
       col    = c("black", colors[seq_along(fits)]),
       lwd = 2, lty = c(1, rep(2, length(fits))),
       cex = 0.75, bty = "n")


# --- MLE Parameter Estimates ---------------------------------
param_table <- as.data.frame(best_model$res)
param_table$Parameter <- rownames(param_table)
all_cols  <- colnames(param_table)
est_col   <- grep("^est$",          all_cols, value = TRUE)
lower_col <- grep("l.*%|lower|L95", all_cols, value = TRUE,
                  ignore.case = TRUE)[1]
upper_col <- grep("u.*%|upper|U95", all_cols, value = TRUE,
                  ignore.case = TRUE)[1]
se_col    <- grep("^se$",           all_cols, value = TRUE)
keep <- c("Parameter", est_col, lower_col, upper_col, se_col)
keep <- keep[keep %in% colnames(param_table)]
param_table <- param_table[, keep]
colnames(param_table) <- c("Parameter","Estimate",
                           "Lower_95CI","Upper_95CI",
                           "SE")[1:ncol(param_table)]
param_table[, -1] <- round(param_table[, -1], 4)
rownames(param_table) <- NULL
print(param_table)


# --- Life Functions Overall ----------------------------------
S_t <- summary(best_model, t = time_seq, type = "survival")[[1]]$est
h_t <- summary(best_model, t = time_seq, type = "hazard")[[1]]$est
f_t <- h_t * S_t
F_t <- 1 - S_t

if (best_model_name == "Weibull") {
  shape  <- best_model$res["shape", "est"]
  scale  <- best_model$res["scale", "est"]
  mean_T <- scale * gamma(1 + 1/shape)
  var_T  <- scale^2 * (gamma(1 + 2/shape) - (gamma(1 + 1/shape))^2)
  cat(sprintf("\nWeibull shape k:      %.4f\n", shape))
  cat(sprintf("Weibull scale lambda: %.4f\n", scale))
  cat(sprintf("Mean survival time:   %.2f days\n", mean_T))
  cat(sprintf("Variance:             %.2f days^2\n", var_T))
  cat(sprintf("SD:                   %.2f days\n", sqrt(var_T)))
}

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
plot(time_seq, S_t, type="l", col="#1E88E5", lwd=2.5,
     xlab="Time (days)", ylab="S(t)",
     main=paste("Survival Function -", best_model_name),
     ylim=c(0,1)); abline(h=0.5, lty=2, col="grey60"); grid()
plot(time_seq, h_t, type="l", col="#E53935", lwd=2.5,
     xlab="Time (days)", ylab="h(t)",
     main=paste("Hazard Function -", best_model_name)); grid()
plot(time_seq, f_t, type="l", col="#43A047", lwd=2.5,
     xlab="Time (days)", ylab="f(t)",
     main=paste("Density Function -", best_model_name)); grid()
plot(time_seq, F_t, type="l", col="#8E24AA", lwd=2.5,
     xlab="Time (days)", ylab="F(t)",
     main=paste("CDF -", best_model_name), ylim=c(0,1))
abline(h=0.5, lty=2, col="grey60"); grid()
par(mfrow = c(1, 1))


# --- Life Functions by Treatment Group -----------------------
rats_ctrl <- subset(rats_clean, treatment == "Control")
rats_trt  <- subset(rats_clean, treatment == "Treated")
best_dist <- dist_list[[best_model_name]]

fit_ctrl <- flexsurvreg(Surv(time, event) ~ 1,
                        data = rats_ctrl, dist = best_dist)
fit_trt  <- flexsurvreg(Surv(time, event) ~ 1,
                        data = rats_trt,  dist = best_dist)

time_seq2 <- seq(0.1, max(rats_clean$time), by = 1)

S_ctrl <- summary(fit_ctrl, t=time_seq2, type="survival")[[1]]$est
h_ctrl <- summary(fit_ctrl, t=time_seq2, type="hazard")[[1]]$est
f_ctrl <- h_ctrl * S_ctrl; F_ctrl <- 1 - S_ctrl

S_trt  <- summary(fit_trt, t=time_seq2, type="survival")[[1]]$est
h_trt  <- summary(fit_trt, t=time_seq2, type="hazard")[[1]]$est
f_trt  <- h_trt * S_trt; F_trt <- 1 - S_trt

if (best_model_name == "Weibull") {
  shape_ctrl <- fit_ctrl$res["shape","est"]
  scale_ctrl <- fit_ctrl$res["scale","est"]
  mean_ctrl  <- scale_ctrl * gamma(1 + 1/shape_ctrl)
  var_ctrl   <- scale_ctrl^2 * (gamma(1+2/shape_ctrl) -
                                  gamma(1+1/shape_ctrl)^2)
  shape_trt  <- fit_trt$res["shape","est"]
  scale_trt  <- fit_trt$res["scale","est"]
  mean_trt   <- scale_trt * gamma(1 + 1/shape_trt)
  var_trt    <- scale_trt^2 * (gamma(1+2/shape_trt) -
                                 gamma(1+1/shape_trt)^2)
  print(data.frame(
    Group   = c("Control","Treated"),
    Shape_k = round(c(shape_ctrl, shape_trt), 4),
    Scale   = round(c(scale_ctrl, scale_trt), 4),
    Mean    = round(c(mean_ctrl, mean_trt), 2),
    Var     = round(c(var_ctrl, var_trt), 2),
    SD      = round(c(sqrt(var_ctrl), sqrt(var_trt)), 2)
  ))
  cat(sprintf("Mean difference: %.1f days\n", mean_ctrl - mean_trt))
}

cols2 <- c("Control"="#2196F3","Treated"="#F44336")
par(mfrow=c(2,2), mar=c(4.5,4.5,3.5,1))
plot(time_seq2, S_ctrl, type="l", col=cols2[1], lwd=2.5,
     ylim=c(0,1), xlab="Time (days)", ylab="S(t)",
     main="Survival by Group")
lines(time_seq2, S_trt, col=cols2[2], lwd=2.5)
abline(h=0.5, lty=2, col="grey60")
legend("topright", legend=names(cols2), col=cols2, lwd=2, bty="n")
grid()
plot(time_seq2, h_ctrl, type="l", col=cols2[1], lwd=2.5,
     xlab="Time (days)", ylab="h(t)", main="Hazard by Group")
lines(time_seq2, h_trt, col=cols2[2], lwd=2.5)
legend("topleft", legend=names(cols2), col=cols2, lwd=2, bty="n")
grid()
plot(time_seq2, f_ctrl, type="l", col=cols2[1], lwd=2.5,
     xlab="Time (days)", ylab="f(t)", main="Density by Group")
lines(time_seq2, f_trt, col=cols2[2], lwd=2.5)
legend("topright", legend=names(cols2), col=cols2, lwd=2, bty="n")
grid()
plot(time_seq2, F_ctrl, type="l", col=cols2[1], lwd=2.5,
     ylim=c(0,1), xlab="Time (days)", ylab="F(t)",
     main="CDF by Group")
lines(time_seq2, F_trt, col=cols2[2], lwd=2.5)
abline(h=0.5, lty=2, col="grey60")
legend("bottomright", legend=names(cols2), col=cols2, lwd=2, bty="n")
grid()
par(mfrow=c(1,1))


# --- PP-Plot -------------------------------------------------
km_fit     <- survfit(surv_obj ~ 1, data = rats_clean)
km_times   <- km_fit$time
km_surv    <- km_fit$surv
param_surv <- summary(best_model, t=km_times,
                      type="survival")[[1]]$est

plot(1-km_surv, 1-param_surv, pch=16, col="#1E88E5", cex=0.8,
     xlab="Empirical CDF (KM)",
     ylab=paste("Fitted CDF (", best_model_name, ")"),
     main=paste("PP-Plot: KM vs", best_model_name),
     xlim=c(0,1), ylim=c(0,1))
abline(0, 1, col="red", lwd=2, lty=2)
legend("bottomright", legend="Perfect fit",
       col="red", lty=2, lwd=2, bty="n")
grid()


# --- Frailty Model -------------------------------------------
frailty_model <- coxph(
  Surv(time, event) ~ treatment +
    frailty(litter, distribution = "gamma"),
  data = rats_clean
)
print(summary(frailty_model))

hr <- exp(coef(frailty_model)["treatmentTreated"])
cat(sprintf("\nFrailty HR (Treated vs Control): %.4f\n", hr))

tryCatch({
  theta <- frailty_model$history[[1]]$theta
  cat(sprintf("Frailty variance theta: %.4f\n", theta))
  if (theta > 0.1) {
    cat("Meaningful within-litter clustering detected.\n")
  }
}, error = function(e) {
  cat("Check frailty_model$history manually.\n")
})


# --- Save Plots to outputs/ ----------------------------------
if (!dir.exists("outputs")) dir.create("outputs")

png("outputs/km_by_treatment.png", width=900, height=700, res=120)
print(ggsurvplot(km_treatment, data=rats_clean, conf.int=TRUE,
  risk.table=TRUE, pval=TRUE, palette=c("#2E86AB","#F18F01"),
  legend.labs=c("Control","Treated"), ggtheme=theme_minimal(),
  title="KM Curves by Treatment", xlab="Time (days)",
  surv.median.line="none"))
dev.off()

png("outputs/schoenfeld_residuals.png",
    width=900, height=500, res=120)
par(mfrow=c(1,2))
plot(ph_test, var=1, main="Schoenfeld - Treatment",
     ylab="Beta(t)", xlab="Time (days)", col="#F18F01", lwd=2)
abline(h=coef(cox_multi)[1], col="red", lty=2, lwd=2)
plot(ph_test, var=2, main="Schoenfeld - Sex",
     ylab="Beta(t)", xlab="Time (days)", col="#A23B72", lwd=2)
abline(h=coef(cox_multi)[2], col="red", lty=2, lwd=2)
par(mfrow=c(1,1))
dev.off()

png("outputs/life_functions_overall.png",
    width=900, height=700, res=120)
par(mfrow=c(2,2), mar=c(4.5,4.5,3,1))
plot(time_seq, S_t, type="l", col="#1E88E5", lwd=2.5,
     xlab="Time (days)", ylab="S(t)",
     main="Survival Function", ylim=c(0,1)); grid()
plot(time_seq, h_t, type="l", col="#E53935", lwd=2.5,
     xlab="Time (days)", ylab="h(t)",
     main="Hazard Function"); grid()
plot(time_seq, f_t, type="l", col="#43A047", lwd=2.5,
     xlab="Time (days)", ylab="f(t)",
     main="Density Function"); grid()
plot(time_seq, F_t, type="l", col="#8E24AA", lwd=2.5,
     xlab="Time (days)", ylab="F(t)",
     main="CDF", ylim=c(0,1)); grid()
par(mfrow=c(1,1))
dev.off()

png("outputs/life_functions_by_group.png",
    width=900, height=700, res=120)
par(mfrow=c(2,2), mar=c(4.5,4.5,3.5,1))
plot(time_seq2, S_ctrl, type="l", col="#2196F3", lwd=2.5,
     ylim=c(0,1), xlab="Time (days)", ylab="S(t)",
     main="Survival by Group")
lines(time_seq2, S_trt, col="#F44336", lwd=2.5)
legend("topright", legend=c("Control","Treated"),
       col=c("#2196F3","#F44336"), lwd=2, bty="n"); grid()
plot(time_seq2, h_ctrl, type="l", col="#2196F3", lwd=2.5,
     xlab="Time (days)", ylab="h(t)", main="Hazard by Group")
lines(time_seq2, h_trt, col="#F44336", lwd=2.5)
legend("topleft", legend=c("Control","Treated"),
       col=c("#2196F3","#F44336"), lwd=2, bty="n"); grid()
plot(time_seq2, f_ctrl, type="l", col="#2196F3", lwd=2.5,
     xlab="Time (days)", ylab="f(t)", main="Density by Group")
lines(time_seq2, f_trt, col="#F44336", lwd=2.5)
legend("topright", legend=c("Control","Treated"),
       col=c("#2196F3","#F44336"), lwd=2, bty="n"); grid()
plot(time_seq2, F_ctrl, type="l", col="#2196F3", lwd=2.5,
     ylim=c(0,1), xlab="Time (days)", ylab="F(t)",
     main="CDF by Group")
lines(time_seq2, F_trt, col="#F44336", lwd=2.5)
legend("bottomright", legend=c("Control","Treated"),
       col=c("#2196F3","#F44336"), lwd=2, bty="n"); grid()
par(mfrow=c(1,1))
dev.off()

cat("\nAll plots saved to outputs/\n")
cat("Analysis complete!\n")
```


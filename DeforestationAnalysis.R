#======================================================================
# Assessing Socioeconomic and Fiscal Drivers of Deforestation in Indonesia Using Principal Component Analysis and Bayesian Modeling  
# Alfiyyah Hasanah
# R script for EDA, PCA, Bayesian models, and sensitivity analysis
#======================================================================

#---------------------------
# SETUP & LIBRARIES
#---------------------------
getwd()

library(dplyr)
library(ggplot2)
library(broom)
library(car)
library(bayesplot)
library(rstanarm)
library(tidybayes)
library(broom.mixed)
library(loo)
library(patchwork)
library(corrplot)
library(GGally)
library(e1071)

#---------------------------
# DATA IMPORT & SIMPLE PREVIEW
#---------------------------
df <- read.csv("DATA/Indonesia_Deforestation_MergedData.csv", header = TRUE, sep = ";")

# Quick checks
str(df)
sapply(df, function(x) sum(is.na(x)))
summary(df$TreeCoverLoss)
table(df$Year)

#---------------------------
# EXPLORATORY DATA ANALYSIS (EDA)
#---------------------------
# Basic distribution and quantiles for tree cover loss
hist(df$TreeCoverLoss, breaks = 20,
     main = "Histogram: Tree Cover Loss (ha)",
     xlab = "TreeCoverLoss")
quantile(df$TreeCoverLoss, probs = c(0, .25, .5, .75, .9, .99))

# Remove Year column for distribution plots
df_no_year <- df %>% select(-Year)

# Create list of ggplot histograms for each numeric variable
plot_list <- lapply(names(df_no_year), function(v) {
  ggplot(df, aes(x = .data[[v]])) +
    geom_histogram(fill = "skyblue", color = "white", bins = 30) +
    labs(title = v, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8))
})

# Combine all distribution plots
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 3)
combined_plot + plot_annotation(title = "Distributions of Key Variables")  # save or export as Figure X

# Correlation heatmap (exclude Year column)
num_vars <- df
if ("Year" %in% names(num_vars)) num_vars_corr <- num_vars %>% select(-Year) else num_vars_corr <- num_vars
cor_mat <- cor(num_vars_corr, use = "pairwise.complete.obs")
corrplot(cor_mat,
         method = "color",
         type = "upper",
         tl.cex = 0.6,
         tl.col = "black",
         tl.srt = 45)

# Pairwise scatter matrix for key predictors
pairs(~ TreeCoverLoss + PopulationTotal + ForestArea + ProtectedArea + UrbanGrowth +
        EnergyUsePerCapita + GDPGrowth + AgriValue,
      data = df, main = "Scatterplot Matrix of Key Variables")

# Time series of tree cover loss
ggplot(df, aes(x = Year, y = TreeCoverLoss)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "orange") +
  labs(title = "Tree Cover Loss Over Time", x = "Year", y = "Tree Cover Loss (ha)") +
  theme_minimal()

#---------------------------
# CHECK LINEAR MODEL RESIDUALS (ASSUMPTION DIAGNOSTICS)
#---------------------------
# Fit a frequentist linear model for residual diagnostics before Bayesian modeling
lm_model <- lm(log(TreeCoverLoss) ~ PopulationTotal + ForestArea +
                 ProtectedArea + UrbanGrowth + EnergyUsePerCapita +
                 GDPGrowth + AgriValue + FDIInflow +
                 ODAUSD + Flood + Fire + Drought, data = df)

res <- residuals(lm_model)

# Residual plots and normality checks
hist(res, main = "Histogram of Residuals", xlab = "Residuals", col = "lightblue")
qqnorm(res); qqline(res, col = "red")
shapiro.test(res)   # Shapiro-Wilk normality test for residuals
plot(fitted(lm_model), res, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values"); abline(h = 0, col = "red")
plot(res, type = "o", main = "Residuals over Time")
acf(res, main = "ACF of Residuals")

#---------------------------
# LOG TRANSFORMATION (response)
#---------------------------
df <- df %>% mutate(log_loss = log(TreeCoverLoss))

p1 <- ggplot(df, aes(x = TreeCoverLoss)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "white") +
  ggtitle("Tree Cover Loss (Raw)") + xlab("Loss (ha)")

p2 <- ggplot(df, aes(x = log_loss)) +
  geom_histogram(bins = 20, fill = "lightcoral", color = "white") +
  ggtitle("Tree Cover Loss (Log Transformed)") + xlab("log(Loss)")

p1 + p2  # place side-by-side (placeholder: Figure X)

# Compute skewness to show effect of transform
skewness(df$TreeCoverLoss, na.rm = TRUE)
skewness(log(df$TreeCoverLoss), na.rm = TRUE)

#---------------------------
# STANDARDIZATION (z-scores) FOR PREDICTORS
#---------------------------
pred_vars <- c("PopulationTotal","ForestArea","ProtectedArea","UrbanGrowth",
               "EnergyUsePerCapita","GDPGrowth","AgriValue",
               "FDIInflow","ODAUSD","Flood","Fire","Drought")

pred_vars <- pred_vars[pred_vars %in% names(df)]

df_std <- df %>%
  mutate(across(all_of(pred_vars),
                ~ (.-mean(., na.rm = TRUE))/sd(., na.rm = TRUE),
                .names = "z_{col}"))

z_preds <- paste0("z_", pred_vars)

#---------------------------
# MULTICOLLINEARITY CHECK (VIF and correlation)
#---------------------------
lm_check <- lm(log_loss ~ z_PopulationTotal + z_ForestArea + z_ProtectedArea +
                 z_UrbanGrowth + z_EnergyUsePerCapita + z_GDPGrowth +
                 z_AgriValue + z_FDIInflow + z_ODAUSD +
                 z_Flood + z_Fire + z_Drought,
               data = df_std)

# Variance Inflation Factors
vif_vals <- car::vif(lm_check)
round(vif_vals, 2)

# Correlation matrix of standardized predictors
round(cor(df_std[, z_preds], use = "pairwise.complete.obs"), 2)

#---------------------------
# PCA FOR HIGHLY CORRELATED BLOCK
#---------------------------
pca_vars <- z_preds
pca_data <- df_std %>% select(all_of(pca_vars)) %>% na.omit()
pca_block <- prcomp(pca_data, scale. = FALSE)

# PCA summaries and loadings
summary(pca_block)
round(pca_block$rotation[, 1:3], 2)

# Scree plot with pink line and labeled PCs
qplot(1:length(pca_block$sdev),
      pca_block$sdev^2 / sum(pca_block$sdev^2)) +
  geom_line(color = "#FF1493", size = 1) +
  geom_point(color = "#FF1493", size = 3) +
  scale_x_continuous(breaks = 1:length(pca_block$sdev),
                     labels = paste0("PC", 1:length(pca_block$sdev))) +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Scree Plot") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        panel.grid.minor = element_blank())

# Save scores back to df_std (PC1, PC2, PC3)
scores <- as.data.frame(pca_block$x)
df_std$PC1 <- NA; df_std$PC2 <- NA; df_std$PC3 <- NA
df_std$PC1[rownames(df_std) %in% rownames(scores)] <- scores$PC1
df_std$PC2[rownames(df_std) %in% rownames(scores)] <- scores$PC2
df_std$PC3[rownames(df_std) %in% rownames(scores)] <- scores$PC3

#---------------------------
# LAGGED VARIABLES
#---------------------------
df_std <- df_std %>%
  arrange(Year) %>%
  mutate(loss_lag1 = lag(TreeCoverLoss, 1),
         log_loss_lag1 = lag(log_loss, 1))

plot(df_std$Year, df_std$loss_lag1, type = "b", main = "Lag-1 Loss (ha)")

#---------------------------
# BAYESIAN MODELS (rstanarm::stan_glm)
#---------------------------

# Model A: log_loss ~ all standardized predictors + lag
formula_A <- log_loss ~ z_PopulationTotal + z_ForestArea + z_ProtectedArea +
  z_UrbanGrowth + z_EnergyUsePerCapita + z_GDPGrowth +
  z_AgriValue + z_FDIInflow + z_ODAUSD +
  z_Flood + z_Fire + z_Drought + log_loss_lag1

prior_model_A <- stan_glm(formula_A, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_A <- update(prior_model_A, prior_PD = FALSE)
tidy(model_A, conf.int = TRUE, conf.level = 0.95)
pp_check(model_A)
loo_A <- loo(model_A)

# Model B: log_loss ~ PC1
formula_B <- log_loss ~ PC1

prior_model_B <- stan_glm(formula_B, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_B <- update(prior_model_B, prior_PD = FALSE)
tidy(model_B, conf.int = TRUE, conf.level = 0.95)
pp_check(model_B)
loo_B <- loo(model_B)

# Model C: log_loss ~ PC1 + PC2 + PC3
formula_C <- log_loss ~ PC1 + PC2 + PC3

prior_model_C <- stan_glm(formula_C, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_C <- update(prior_model_C, prior_PD = FALSE)
tidy(model_C, conf.int = TRUE, conf.level = 0.95)
pp_check(model_C)
loo_C <- loo(model_C)

# Model D: log_loss ~ PC1 + PC2
formula_D <- log_loss ~ PC1 + PC2

prior_model_D <- stan_glm(formula_D, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_D <- update(prior_model_D, prior_PD = FALSE)
tidy(model_D, conf.int = TRUE, conf.level = 0.95)
pp_check(model_D)
loo_D <- loo(model_D)

# Model E: log_loss ~ PC1 + PC3
formula_E <- log_loss ~ PC1 + PC3

prior_model_E <- stan_glm(formula_E, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_E <- update(prior_model_E, prior_PD = FALSE)
tidy(model_E, conf.int = TRUE, conf.level = 0.95)
pp_check(model_E)
loo_E <- loo(model_E)

# Model F: log_loss ~ PC2 + PC3
formula_F <- log_loss ~ PC2 + PC3

prior_model_F <- stan_glm(formula_F, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_F <- update(prior_model_F, prior_PD = FALSE)
tidy(model_F, conf.int = TRUE, conf.level = 0.95)
pp_check(model_F)
loo_F <- loo(model_F)

# Model G: log_loss ~ PC1 + lag
formula_G <- log_loss ~ PC1 + log_loss_lag1

prior_model_G <- stan_glm(formula_G, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_G <- update(prior_model_G, prior_PD = FALSE)
tidy(model_G, conf.int = TRUE, conf.level = 0.95)
pp_check(model_G)
loo_G <- loo(model_G)

# Model H: log_loss ~ PC1 + Budget (load Budget externally)
budget <- read.csv("DATA/APBLingkungan20082024.csv", header = TRUE, sep = ";")
df_std <- df_std %>%
  left_join(budget %>% select(Year, Budget), by = "Year") %>%
  mutate(z_Budget = as.numeric(scale(Budget)))

formula_H <- log_loss ~ PC1 + z_Budget

prior_model_H <- stan_glm(formula_H, data = df_std, family = gaussian(),
                          prior = normal(0, 10), prior_intercept = normal(0, 5),
                          prior_aux = exponential(1), chains = 4, iter = 10000, seed = 2025)

model_H <- update(prior_model_H, prior_PD = FALSE)
tidy(model_H, conf.int = TRUE, conf.level = 0.95)
pp_check(model_H)
loo_H <- loo(model_H)


#---------------------------
# MODEL COMPARISONS (LOO)
#---------------------------

# All models
tidy(model_A, conf.int = TRUE, conf.level = 0.95); loo_A  # all standardized predictors + lag
tidy(model_B, conf.int = TRUE, conf.level = 0.95); loo_B  # PC1 only
tidy(model_C, conf.int = TRUE, conf.level = 0.95); loo_C  # PC1 + PC2 + PC3
tidy(model_D, conf.int = TRUE, conf.level = 0.95); loo_D  # PC1 + PC2
tidy(model_E, conf.int = TRUE, conf.level = 0.95); loo_E  # PC1 + PC3
tidy(model_F, conf.int = TRUE, conf.level = 0.95); loo_F  # PC2 + PC3
tidy(model_G, conf.int = TRUE, conf.level = 0.95); loo_G  # PC1 + lag
tidy(model_H, conf.int = TRUE, conf.level = 0.95); loo_H  # PC1 + Budget

#---------------------------
# PP_CHECK PLOTS
#---------------------------
pp_A <- pp_check(model_A) + ggtitle("Model A")
pp_B <- pp_check(model_B) + ggtitle("Model B")
pp_C <- pp_check(model_C) + ggtitle("Model C")
pp_D <- pp_check(model_D) + ggtitle("Model D")
pp_E <- pp_check(model_E) + ggtitle("Model E")
pp_F <- pp_check(model_F) + ggtitle("Model F")
pp_G <- pp_check(model_G) + ggtitle("Model G")
pp_H <- pp_check(model_H) + ggtitle("Model H")
(pp_A | pp_B | pp_C | pp_D) / (pp_E | pp_F | pp_G | pp_H)  # preview grid (Figure placeholder)

#---------------------------
# SENSITIVITY ANALYSIS: Model B under alternative priors
#---------------------------
# Define prior sets (as in your original script)
priors_weak <- list(prior = normal(0, 10), prior_intercept = normal(0, 5), prior_aux = exponential(1))
priors_strong <- list(prior = normal(0, 2), prior_intercept = normal(0, 2), prior_aux = exponential(1))
priors_flat <- list(prior = normal(0, 50), prior_intercept = normal(0, 10), prior_aux = exponential(0.1))
priors_student <- list(prior = student_t(df = 3, location = 0, scale = 5),
                       prior_intercept = student_t(df = 3, location = 0, scale = 2),
                       prior_aux = exponential(1))

# Fit model under each prior
# Model B: log_loss ~ PC1
prior_model_weak <- stan_glm(formula_B, data = df_std, family = gaussian(),
                          prior = priors_weak$prior, prior_intercept = priors_weak$prior_intercept,
                          prior_aux = priors_weak$prior_aux, chains = 4, iter = 10000, seed = 2025)

model_weak <- update(prior_model_weak, prior_PD = FALSE)

prior_model_strong <- stan_glm(formula_B, data = df_std, family = gaussian(),
                             prior = priors_strong$prior, prior_intercept = priors_strong$prior_intercept,
                             prior_aux = priors_strong$prior_aux, chains = 4, iter = 10000, seed = 2025)

model_strong <- update(prior_model_strong, prior_PD = FALSE)

prior_model_flat <- stan_glm(formula_B, data = df_std, family = gaussian(),
                               prior = priors_flat$prior, prior_intercept = priors_flat$prior_intercept,
                               prior_aux = priors_flat$prior_aux, chains = 4, iter = 10000, seed = 2025)

model_flat <- update(prior_model_flat, prior_PD = FALSE)

prior_model_student <- stan_glm(formula_B, data = df_std, family = gaussian(),
                             prior = priors_student$prior, prior_intercept = priors_student$prior_intercept,
                             prior_aux = priors_student$prior_aux, chains = 4, iter = 10000, seed = 2025)

model_student <- update(prior_model_student, prior_PD = FALSE)

# Summarize and compare coefficient (PC1) across priors
results <- bind_rows(
  tidy(model_weak, conf.int = TRUE, conf.level = 0.95) %>% mutate(prior = "Weak (N(0,10))"),
  tidy(model_strong, conf.int = TRUE, conf.level = 0.95) %>% mutate(prior = "Strong (N(0,2))"),
  tidy(model_flat, conf.int = TRUE, conf.level = 0.95) %>% mutate(prior = "Flat (N(0,50))"),
  tidy(model_student, conf.int = TRUE, conf.level = 0.95) %>% mutate(prior = "Student-t(3,0,5)")
)

results %>% filter(term == "PC1") %>% select(prior, estimate, conf.low, conf.high)

loo_results <- list(Weak = loo(model_weak), Strong = loo(model_strong),
                    Flat = loo(model_flat), Student = loo(model_student))
loo_compare(loo_results$Weak, loo_results$Strong, loo_results$Flat, loo_results$Student)

# Visualize posterior sensitivity for PC1
mcmc_areas(
  list(
    "Weak (N(0,10))" = as.matrix(model_weak, pars = "PC1"),
    "Strong (N(0,2))" = as.matrix(model_strong, pars = "PC1"),
    "Flat (N(0,50))" = as.matrix(model_flat, pars = "PC1"),
    "Student-t(3,0,5)" = as.matrix(model_student, pars = "PC1")
  ),
  prob = 0.95
) + ggtitle("Posterior Sensitivity of PC1 to Different Priors")

# Posterior predictive checks for each prior model (grid)
pp_weak <- pp_check(model_weak) + ggtitle("Model B with Weakly Informative Prior")
pp_strong <- pp_check(model_strong) + ggtitle("Model B with Strongly Informative Prior")
pp_flat <- pp_check(model_flat) + ggtitle("Model B with Flat Prior")
pp_student <- pp_check(model_student) + ggtitle("Model B with Student-t Prior")
(pp_weak | pp_strong) / (pp_flat | pp_student)

#======================================================================
# End of script
#======================================================================

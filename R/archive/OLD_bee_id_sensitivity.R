# Follow-up checks for the volunteer ID model
#   1. Residual diagnostics (DHARMa)
#   2. Sensitivity analyses on the training effect (glmmTMB)
#   3. GEE cross-check (geepack)
#
# Run bee_id_cleaning.R first: it writes ../output/bee_id_clean.csv.
# Paths are relative to the R/ folder (where renv lives).
#
# Setup (run once, then renv::snapshot()):
#   install.packages(c("geepack", "DHARMa"))

library(readr)
library(dplyr)
library(stringr)
library(glmmTMB)
library(DHARMa)
library(geepack)

# ---------------------------------------------------------------- data
dat <- read_csv("/Users/aubreecurtis/Research/species-volunteer-glmm/bee_id_clean.csv", show_col_types = FALSE) |>
  mutate(training     = factor(training, levels = c("1", "2", "3+")),
         year         = factor(year),
         vid          = factor(vid),
         true_species = factor(true_species))

cat("Rows:", nrow(dat), " Volunteers:", n_distinct(dat$vid), "\n")

# ================================================================ 1. DIAGNOSTICS
m_main <- glmmTMB(
  y ~ training + master_gardener + year + (1 | vid) + (1 | true_species),
  data = dat, family = binomial
)

sim <- simulateResiduals(m_main, n = 1000)

pdf("output/diagnostics.pdf", width = 9, height = 6)
plot(sim)                            # QQ plot + residuals vs. predicted
plotResiduals(sim, dat$training)     # residuals by training level
plotResiduals(sim, dat$true_species) # residuals by species
dev.off()

cat("\n--- Diagnostics (p < 0.05 suggests a problem) ---\n")
print(testUniformity(sim, plot = FALSE))
print(testDispersion(sim, plot = FALSE))
print(testCategorical(sim, dat$true_species, plot = FALSE))
print(testCategorical(sim, dat$training, plot = FALSE))

# ================================================================ 2. SENSITIVITY
# Fits the model and its no-training version, and returns the training odds
# ratios (relative to level 1) with Wald 95% CIs and the likelihood-ratio p.
fit_scenario <- function(d, label,
                         fixed  = "training + master_gardener + year",
                         random = "(1 | vid) + (1 | true_species)") {
  d <- droplevels(d)
  f_full <- as.formula(paste("y ~", fixed, "+", random))
  f_null <- as.formula(paste("y ~", sub("^training \\+ ", "", fixed), "+", random))

  m  <- glmmTMB(f_full, data = d, family = binomial)
  m0 <- glmmTMB(f_null, data = d, family = binomial)
  cf <- summary(m)$coefficients$cond
  p_lrt <- anova(m0, m)[["Pr(>Chisq)"]][2]

  or_ci <- function(term) {
    est <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]
    c(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))
  }
  a <- or_ci("training2"); b <- or_ci("training3+")

  tibble(scenario = label, n_obs = nrow(d), n_vol = n_distinct(d$vid),
         or_2 = a[1], lo_2 = a[2], hi_2 = a[3],
         or_3plus = b[1], lo_3plus = b[2], hi_3plus = b[3],
         p_training = p_lrt)
}

# --- scenario inputs
top_vid <- names(which.max(table(dat$vid)))          # volunteer with the most rows
flip_vids <- dat |>                                   # training level differs within a year
  group_by(vid, year) |>
  summarise(k = n_distinct(training), .groups = "drop") |>
  filter(k > 1) |>
  pull(vid) |> as.character() |> unique()
bee_only <- dat |>                                    # drop fly and wasp true IDs
  filter(!str_detect(str_to_lower(correctid), "fly|flies|wasp"))

cat("\nTop volunteer:", top_vid, "(", sum(dat$vid == top_vid), "rows )\n")
cat("Training-flip volunteers:", paste(flip_vids, collapse = ", "),
    "(", sum(dat$vid %in% flip_vids), "rows )\n")
cat("Rows dropped for bees-only:", nrow(dat) - nrow(bee_only), "\n")

results <- bind_rows(
  fit_scenario(dat, "Main model"),
  fit_scenario(filter(dat, vid != top_vid),
               paste0("Without top volunteer (", top_vid, ")")),
  fit_scenario(filter(dat, !vid %in% flip_vids),
               "Without training-flip volunteers"),
  fit_scenario(bee_only, "Bees only (no fly/wasp true IDs)"),
  fit_scenario(dat, "Species as fixed effect",
               fixed  = "training + master_gardener + year + true_species",
               random = "(1 | vid)")
) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n--- Training effect across scenarios (odds ratios vs. training level 1) ---\n")
print(as.data.frame(results), row.names = FALSE)
write_csv(results, "/Users/aubreecurtis/Research/species-volunteer-glmm/output/sensitivity_summary.csv")

# ================================================================ 3. GEE CROSS-CHECK
# Population-averaged logistic regression with cluster-robust SEs.
# GEE has no random effects, so species enters as a fixed effect.
d_gee <- dat |> arrange(vid)     # geeglm needs rows sorted by cluster

gee_full <- geeglm(y ~ training + master_gardener + year + true_species,
                   id = vid, data = d_gee, family = binomial,
                   corstr = "exchangeable")
gee_null <- geeglm(y ~ master_gardener + year + true_species,
                   id = vid, data = d_gee, family = binomial,
                   corstr = "exchangeable")

cat("\n--- GEE (population-averaged) ---\n")
cf <- summary(gee_full)$coefficients
gee_tab <- data.frame(term = rownames(cf),
                      OR = exp(cf$Estimate),
                      lo = exp(cf$Estimate - 1.96 * cf$Std.err),
                      hi = exp(cf$Estimate + 1.96 * cf$Std.err),
                      p  = cf[["Pr(>|W|)"]])
gee_tab[, 2:5] <- round(gee_tab[, 2:5], 3)
print(gee_tab, row.names = FALSE)

cat("\nWald test for training (both levels together):\n")
print(anova(gee_null, gee_full))

write_csv(gee_tab, "output/gee_summary.csv")

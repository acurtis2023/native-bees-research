# Volunteer species-ID accuracy ~ training
# Cleaning + exclusions + logistic GLMM
#
# Columns expected in the workbook:
#   obs, year, vid, vtrain, mg, surveyid, correctid, morphacc
#
# Setup (run once):
#   install.packages(c("readxl", "dplyr", "stringr", "purrr", "readr",
#                      "glmmTMB", "emmeans", "DHARMa"))

library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

path <- "/Users/aubreecurtis/Research/species-volunteer-glmm/data/FinalStatistics.xlsx"
MIN_SPECIES_N <- 30              # rarer true species are pooled into "Other"

# ---------------------------------------------------------------- read
# Read everything as text: correctid and morphacc mix numbers (9) and words.
raw <- read_excel(path, col_types = "text")

df <- raw |>
  mutate(across(c(obs, year, vid, vtrain), as.integer),
         across(c(mg, surveyid, correctid, morphacc), str_trim))

# ---------------------------------------------------------------- labels
# Standardize spellings. Keys are lower-case with spaces removed.
canon <- c(
  bumblebees = "BumbleBees",   bumblebee = "BumbleBees",
  honeybees = "HoneyBees",
  stripedsweatbees = "StripedSweatBees",
  tinydarkbees = "TinyDarkBees", tinydarkbeess = "TinyDarkBees",
  hairybellybees = "HairyBellyBees",
  hairylegbees = "HairyLegBees",
  greenmetallicbees = "GreenMetallicBees",
  syrphidflies = "SyrphidFlies",
  beeflies = "BeeFlies",
  wasps = "Wasps",             wasp = "Wasps",
  cuckoobees = "CuckooBees",
  flies = "Flies"              # a fly that isn't a named fly category
)

# Turn "HoneyBees, BumbleBees" into the set c("BumbleBees", "HoneyBees").
# Anything not in `canon` (e.g. "HairyBellyBees?") becomes "UNRECOGNIZED".
parse_ids <- function(x) {
  key <- str_to_lower(str_remove_all(x, "\\s"))
  if (key == "honeybeestinydarkbees") {          # missing comma in the source
    return(c("HoneyBees", "TinyDarkBees"))
  }
  out <- unname(canon[str_split(key, ",")[[1]]])
  out[is.na(out)] <- "UNRECOGNIZED"
  sort(unique(out))
}

survey_sets <- map(df$surveyid, parse_ids)
true_sets   <- map(df$correctid, parse_ids)

df <- df |>
  mutate(
    n_survey     = lengths(survey_sets),
    n_true       = lengths(true_sets),
    same_set     = map2_lgl(survey_sets, true_sets, setequal),
    unrecognized = map_lgl(survey_sets, ~ "UNRECOGNIZED" %in% .x) |
                   map_lgl(true_sets,   ~ "UNRECOGNIZED" %in% .x),
    true_species = if_else(n_true == 1, map_chr(true_sets, 1), "Multiple")
  )

# ---------------------------------------------------------------- mg check
# A single change within a year (e.g. No -> Yes) can be real: someone became a
# master gardener mid-year. Status that flips back and forth (Yes, No, Yes ...)
# cannot be real, so in those volunteer-years the minority-value rows are
# treated as errors. Rows are ordered by obs, assumed to follow time order.
df <- df |>
  arrange(obs) |>
  group_by(vid, year) |>
  mutate(mg_changes  = sum(mg != lag(mg), na.rm = TRUE),
         mg_modal    = names(which.max(table(mg))),
         mg_conflict = mg_changes > 1 & mg != mg_modal) |>
  ungroup()

# ---------------------------------------------------------------- exclusions
# Each row gets at most one reason (first match wins). The 9s are excluded by
# team decision; the remaining rules drop rows whose scoring is uncertain.
df <- df |>
  mutate(exclude_reason = case_when(
    correctid == "9" | morphacc == "9"        ~ "Unknown true ID (9)",
    unrecognized                               ~ "Unrecognized label",
    morphacc == "Yes" & !same_set              ~ "Scored Yes but species differ",
    morphacc == "No"  & same_set               ~ "Scored No but IDs match",
    n_true > 1 & n_survey == 1                 ~ "Multi-species true ID, single-species answer",
    mg_conflict                                ~ "Master gardener status flips back and forth within year",
    TRUE                                       ~ NA_character_
  ))

excluded <- df |>
  filter(!is.na(exclude_reason)) |>
  select(obs, year, vid, vtrain, mg, surveyid, correctid, morphacc, exclude_reason)
write_csv(excluded, "/Users/aubreecurtis/Research/species-volunteer-glmm/data/excluded_rows.csv")   # for the methods section / the team

cat("\nRows read:", nrow(df), "\n")
cat("\nExclusions by reason:\n"); print(count(excluded, exclude_reason, sort = TRUE))
cat("\nShare of rows excluded, by training level:\n")
print(df |> group_by(vtrain) |>
        summarise(n = n(), excluded = sum(!is.na(exclude_reason)),
                  pct_excluded = round(100 * excluded / n, 1)))

# ---------------------------------------------------------------- analysis set
dat <- df |>
  filter(is.na(exclude_reason)) |>
  mutate(
    y               = as.integer(morphacc == "Yes"),
    training        = factor(vtrain, levels = 1:3, labels = c("1", "2", "3+")),
    year            = factor(year),
    master_gardener = as.integer(mg == "Yes"),
    vid             = factor(vid)
  )

rare <- dat |> count(true_species) |>
  filter(n < MIN_SPECIES_N, true_species != "Multiple") |> pull(true_species)
dat <- dat |>
  mutate(true_species = factor(if_else(true_species %in% rare, "Other", true_species)))

cat("\nAnalysis set:", nrow(dat), "rows,", n_distinct(dat$vid), "volunteers\n")
print(dat |> group_by(training) |>
        summarise(n = n(), volunteers = n_distinct(vid), accuracy = round(mean(y), 3)))
print(dat |> group_by(true_species) |>
        summarise(n = n(), accuracy = round(mean(y), 3)) |> arrange(desc(n)))

write_csv(dat, "bee_id_clean.csv")

# ---------------------------------------------------------------- model
library(glmmTMB)
library(emmeans)

m <- glmmTMB(
  y ~ training + master_gardener + year + (1 | vid) + (1 | true_species),
  data = dat, family = binomial
)
summary(m)

# Overall test of training (both coefficients together)
m0 <- update(m, . ~ . - training)
anova(m0, m)

# Predicted accuracy by training level, with confidence intervals
emmeans(m, "training", type = "response")
pairs(emmeans(m, "training"), type = "response")   # odds ratios between levels
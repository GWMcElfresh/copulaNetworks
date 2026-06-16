# Example: ImpacTB stratified copula workflow
#
# Prerequisites: run data prep in Stratified_Copula_Networks.Rmd through the MICE
# imputation chunk so `bn_data_imp` exists with disease posteriors merged.
#
# This script shows how to port the analysis using copulaNetworks function arguments
# (no YAML config).

library(copulaNetworks)
library(dplyr)

# ---------------------------------------------------------------------------
# ImpacTB-specific node group styling (optional)
# ---------------------------------------------------------------------------
impactb_node_groups <- function(vars) {
  dplyr::case_when(
    grepl("Protect|Severity", vars) ~ "Disease Outcome",
    grepl("M1|Alv|Myeloid|Total_Myeloid", vars) ~ "Myeloid",
    grepl("CD4|CD8|Activated", vars) ~ "T Cell",
    grepl("Bcell|CD40|Plasma|GC_B|Total_Bcell", vars) ~ "B Cell",
    grepl("MCore|BALT|Lymph|EarlyGran|Necrotic|Intermediate", vars) ~ "Spatial",
    TRUE ~ "Other"
  )
}

# ---------------------------------------------------------------------------
# Step 0: prepare data and strata
# ---------------------------------------------------------------------------
# Assumes bn_data_imp is in the global environment from the Rmd pipeline.

if (!exists("bn_data_imp")) {
  stop(
    "bn_data_imp not found. Source the ImpacTB Rmd through imputation first, ",
    "or load a saved checkpoint."
  )
}

prep <- prepare_copula_data(
  data = bn_data_imp,
  id_cols = c("SubjectId", "cDNA_ID"),
  strata_cols = c("Vaccine", "Timepoint", "Tissue", "Challenge"),
  exclude_cols = c("CFU_Homogenate", "LungPathScore", "Protection"),
  strata_specs = list(
    cohort = list(
      mutate = quote(Cohort = dplyr::case_when(
        grepl("BCG", Vaccine) ~ "BCG",
        grepl("RhCMV-TB", Vaccine) ~ "RhCMV-TB",
        Vaccine %in% c("Unvaccinated", "RhCMV/Gag") ~ "Control",
        TRUE ~ "Other"
      )),
      group_by = "Cohort",
      min_n = 10
    ),
    late = list(
      filter = quote(Timepoint %in% c("Day 21", "Day 28")),
      group_by = c("Cohort", "Timepoint"),
      name_sep = " | ",
      min_n = 10
    )
  ),
  out_dir = "checkpoints/copula_impactb"
)

# ---------------------------------------------------------------------------
# Step 1: fit copula models
# ---------------------------------------------------------------------------
fits <- fit_copula_strata(
  prep,
  method = "stars",
  nlambda = 40,
  stars_thresh = 0.1,
  min_n = 10,
  include_full = TRUE,
  out_dir = "checkpoints/copula_impactb"
)

# ---------------------------------------------------------------------------
# Step 2: single-stratum diagnostics
# ---------------------------------------------------------------------------
plot_all_strata(
  fits,
  out_dir = "checkpoints/copula_impactb/figures",
  node_groups = impactb_node_groups,
  min_pcor = 0.01
)

# ---------------------------------------------------------------------------
# Step 3: structure comparison across strata
# ---------------------------------------------------------------------------
if (all(c("cohort::BCG", "cohort::RhCMV-TB") %in% names(fits$fits))) {
  cmp_bcg_rhcmv <- compare_two_strata(
    fits$fits[["cohort::BCG"]],
    fits$fits[["cohort::RhCMV-TB"]],
    label_a = "BCG",
    label_b = "RhCMV-TB",
    matrices = c("pcor", "copula_cor")
  )

  plot_stratum_comparison(
    cmp_bcg_rhcmv,
    out_dir = "checkpoints/copula_impactb/comparisons/BCG_vs_RhCMV-TB"
  )
}

# ---------------------------------------------------------------------------
# One-call alternative
# ---------------------------------------------------------------------------
# result <- run_copula_pipeline(
#   data = bn_data_imp,
#   id_cols = c("SubjectId", "cDNA_ID"),
#   strata_cols = c("Vaccine", "Timepoint", "Tissue", "Challenge"),
#   exclude_cols = c("CFU_Homogenate", "LungPathScore", "Protection"),
#   strata_specs = list(
#     cohort = list(
#       mutate = quote(Cohort = dplyr::case_when(
#         grepl("BCG", Vaccine) ~ "BCG",
#         grepl("RhCMV-TB", Vaccine) ~ "RhCMV-TB",
#         Vaccine %in% c("Unvaccinated", "RhCMV/Gag") ~ "Control",
#         TRUE ~ "Other"
#       )),
#       group_by = "Cohort",
#       min_n = 10
#     )
#   ),
#   compare_pairs = list(c("cohort::BCG", "cohort::RhCMV-TB")),
#   out_dir = "checkpoints/copula_impactb_pipeline",
#   node_groups = impactb_node_groups
# )

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

prep <- PrepareCopulaData(
  data = bn_data_imp,
  idCols = c("SubjectId", "cDNA_ID"),
  strataCols = c("Vaccine", "Timepoint", "Tissue", "Challenge"),
  excludeCols = c("CFU_Homogenate", "LungPathScore", "Protection"),
  strataSpecs = list(
    cohort = list(
      mutate = quote(Cohort = dplyr::case_when(
        grepl("BCG", Vaccine) ~ "BCG",
        grepl("RhCMV-TB", Vaccine) ~ "RhCMV-TB",
        Vaccine %in% c("Unvaccinated", "RhCMV/Gag") ~ "Control",
        TRUE ~ "Other"
      )),
      group_by = "Cohort",
      minN = 10
    ),
    late = list(
      filter = quote(Timepoint %in% c("Day 21", "Day 28")),
      group_by = c("Cohort", "Timepoint"),
      nameSep = " | ",
      minN = 10
    )
  ),
  outDir = "checkpoints/copula_impactb"
)

# ---------------------------------------------------------------------------
# Step 1: fit copula models
# ---------------------------------------------------------------------------
fits <- FitCopulaStrata(
  prep,
  method = "stars",
  nlambda = 40,
  starsThresh = 0.1,
  minN = 10,
  includeFull = TRUE,
  outDir = "checkpoints/copula_impactb"
)

# ---------------------------------------------------------------------------
# Step 2: single-stratum diagnostics
# ---------------------------------------------------------------------------
PlotAllStrata(
  fits,
  outDir = "checkpoints/copula_impactb/figures",
  nodeGroups = impactb_node_groups,
  minPcor = 0.01
)

# ---------------------------------------------------------------------------
# Step 3: structure comparison across strata
# ---------------------------------------------------------------------------
if (all(c("cohort::BCG", "cohort::RhCMV-TB") %in% names(fits$fits))) {
  cmp_bcg_rhcmv <- CompareTwoStrata(
    fits$fits[["cohort::BCG"]],
    fits$fits[["cohort::RhCMV-TB"]],
    labelA = "BCG",
    labelB = "RhCMV-TB",
    matrices = c("pcor", "copulaCor")
  )

  PlotStratumComparison(
    cmp_bcg_rhcmv,
    outDir = "checkpoints/copula_impactb/comparisons/BCG_vs_RhCMV-TB"
  )
}

# ---------------------------------------------------------------------------
# One-call alternative
# ---------------------------------------------------------------------------
# result <- RunCopulaPipeline(
#   data = bn_data_imp,
#   idCols = c("SubjectId", "cDNA_ID"),
#   strataCols = c("Vaccine", "Timepoint", "Tissue", "Challenge"),
#   excludeCols = c("CFU_Homogenate", "LungPathScore", "Protection"),
#   strataSpecs = list(
#     cohort = list(
#       mutate = quote(Cohort = dplyr::case_when(
#         grepl("BCG", Vaccine) ~ "BCG",
#         grepl("RhCMV-TB", Vaccine) ~ "RhCMV-TB",
#         Vaccine %in% c("Unvaccinated", "RhCMV/Gag") ~ "Control",
#         TRUE ~ "Other"
#       )),
#       group_by = "Cohort",
#       minN = 10
#     )
#   ),
#   comparePairs = list(c("cohort::BCG", "cohort::RhCMV-TB")),
#   outDir = "checkpoints/copula_impactb_pipeline",
#   nodeGroups = impactb_node_groups
# )

#' Run the full two-phase factor-vine copula pipeline
#'
#' Phase 1 fits a factor copula prior on a large reference cohort; Phase 2
#' updates with graphical and/or vine models on a small sample; optional
#' consistency testing via simulation and/or Bayesian updating.
#'
#' @param priorData Large reference data frame (N x variables).
#' @param updateData Small update data frame (n x variables).
#' @param nodeCols Character vector of node column names.
#' @param nFactors Number of latent factors for Phase 1 (1 or 2).
#' @param linkingCopula Linking copula family per variable (Phase 1).
#' @param nQuad Quadrature points for FactorCopula integration.
#' @param phase2Method One of `"graphical"`, `"vine"`, or `"both"`.
#' @param testMethod One of `"simulation"`, `"bayes"`, `"both"`, or `"none"`.
#' @param nRep Simulation replicates for [TestPriorConsistency()].
#' @param nlambda Glasso path length for graphical update.
#' @param glassoMethod Lambda selection for graphical update.
#' @param starsThresh StARS threshold.
#' @param outDir Optional directory for RDS checkpoints.
#' @param seed Optional random seed for simulation test.
#' @return List with `priorFit`, `updateFit`, `comparison`, `consistencyTest`,
#'   and `meta`.
#' @export
RunFactorVinePipeline <- function(priorData,
                                  updateData,
                                  nodeCols,
                                  nFactors = 1L,
                                  linkingCopula = "bvn",
                                  nQuad = 25L,
                                  phase2Method = c("graphical", "vine", "both"),
                                  testMethod = c("simulation", "bayes", "both", "none"),
                                  nRep = 500L,
                                  nlambda = 40,
                                  glassoMethod = c("stars", "ebic"),
                                  starsThresh = 0.1,
                                  outDir = NULL,
                                  seed = NULL) {
  phase2_method <- match.arg(phase2Method)
  test_method <- match.arg(testMethod)
  glasso_method <- match.arg(glassoMethod)
  t_start <- Sys.time()
  warnings_captured <- character(0)

  withCallingHandlers(
    {
      prior_fit <- FitFactorCopulaPrior(
        priorData,
        nodeCols = nodeCols,
        nFactors = nFactors,
        linkingCopula = linkingCopula,
        nQuad = nQuad
      )

      factor_check <- CheckFactorStructure(prior_fit$impliedCor, nFactors = nFactors)

      update_fit <- FitCopulaUpdate(
        updateData,
        priorFit = prior_fit,
        method = phase2_method,
        nlambda = nlambda,
        glassoMethod = glasso_method,
        starsThresh = starsThresh
      )

      comparison <- NULL
      if (!is.null(update_fit$graphical)) {
        comparison <- ComparePriorToUpdate(prior_fit, update_fit)
      }

      consistency_test <- list()
      if (test_method %in% c("simulation", "both")) {
        consistency_test$simulation <- TestPriorConsistency(
          prior_fit,
          updateData,
          nRep = nRep,
          seed = seed
        )
      }
      if (test_method %in% c("bayes", "both") && requireNamespace("cmdstanr", quietly = TRUE)) {
        consistency_test$bayes <- FitBayesianFactorUpdate(
          prior_fit,
          updateData,
          computeBayesFactor = requireNamespace("bridgesampling", quietly = TRUE),
          seed = seed
        )
        if (!is.null(consistency_test$bayes$bayesFactor)) {
          consistency_test$bayesFactor <- consistency_test$bayes$bayesFactor$logMargLik
        }
      } else if (test_method %in% c("bayes", "both")) {
        warning("cmdstanr not installed - skipping Bayesian consistency test.", call. = FALSE)
      }

      if (test_method == "simulation" && !is.null(consistency_test$simulation)) {
        consistency_test$pValue <- consistency_test$simulation$pValue
      }

      meta <- list(
        elapsedSec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
        factorStructure = factor_check,
        phase2Method = phase2_method,
        testMethod = test_method,
        packageVersion = utils::packageVersion("copulaNetworks"),
        warnings = warnings_captured
      )

      result <- list(
        priorFit = prior_fit,
        updateFit = update_fit,
        comparison = comparison,
        consistencyTest = consistency_test,
        meta = meta
      )

      if (!is.null(outDir)) {
        SaveCheckpoint(result, outDir, filename = "factor_vine_result.rds")
      }

      result
    },
    warning = function(w) {
      warnings_captured <<- c(warnings_captured, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

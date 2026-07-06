.factor_update_stan_cache <- new.env(parent = emptyenv())

#' Compile or retrieve cached Stan model for factor update
#'
#' @return A `CmdStanModel` object.
#' @keywords internal
get_factor_update_stan_model <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required for FitBayesianFactorUpdate(). ",
      "Install with install.packages('cmdstanr'), then run ",
      "cmdstanr::install_cmdstan() to install CmdStan.",
      call. = FALSE
    )
  }
  if (!is.null(.factor_update_stan_cache$mod)) {
    return(.factor_update_stan_cache$mod)
  }
  stan_file <- system.file("stan", "factor_update.stan", package = "copulaNetworks")
  if (!file.exists(stan_file)) {
    stop("Stan model file not found in inst/stan/.", call. = FALSE)
  }
  .factor_update_stan_cache$mod <- cmdstanr::cmdstan_model(stan_file)
  .factor_update_stan_cache$mod
}

#' Bayesian Gaussian copula update with priors centered at Phase 1
#'
#' Requires optional **cmdstanr** (and **bridgesampling** for Bayes factors).
#' ponytail: single-factor Gaussian copula via multivariate normal on scores.
#'
#' @param priorFit Output of [FitFactorCopulaPrior()].
#' @param updateData Data frame with node columns for the update cohort.
#' @param chains Number of MCMC chains.
#' @param iter Total iterations per chain (warmup + sampling; split evenly).
#' @param computeBayesFactor If `TRUE` and bridgesampling is available, compute
#'   a marginal-likelihood estimate via bridge sampling on the fitted model.
#' @param seed Optional random seed.
#' @return List with `fit` (`CmdStanMCMC`), `summary`, and optional `bayesFactor`.
#' @export
FitBayesianFactorUpdate <- function(priorFit,
                                    updateData,
                                    chains = 2L,
                                    iter = 1000L,
                                    computeBayesFactor = FALSE,
                                    seed = NULL) {
  node_cols <- priorFit$nodeCols
  uniform_matrix <- ApplyMarginalSpec(updateData, priorFit$marginalSpec, nodeCols = node_cols)
  z_matrix <- qnorm(uniform_matrix)
  z_matrix <- z_matrix[, node_cols, drop = FALSE]

  mu0 <- colMeans(z_matrix, na.rm = TRUE)
  r0 <- priorFit$impliedCor[node_cols, node_cols, drop = FALSE]
  diag(r0) <- 1

  stan_data <- list(
    N = nrow(z_matrix),
    D = ncol(z_matrix),
    Z = z_matrix,
    mu0 = mu0,
    R0 = r0
  )

  iter_warmup <- iter %/% 2L
  iter_sampling <- iter - iter_warmup

  stan_mod <- get_factor_update_stan_model()
  sample_args <- list(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  if (!is.null(seed)) {
    sample_args$seed <- seed
  }
  fit <- do.call(stan_mod$sample, sample_args)

  fit_summary <- fit$summary()

  result <- list(
    fit = fit,
    summary = fit_summary,
    stanData = stan_data
  )

  if (isTRUE(computeBayesFactor) && requireNamespace("bridgesampling", quietly = TRUE)) {
    result$bayesFactor <- tryCatch(
      {
        bs <- bridgesampling::bridge_sampler(fit, silent = TRUE)
        list(logMargLik = bs$logml, method = "bridge_sampler")
      },
      error = function(e) {
        warning("Bayes factor computation failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  } else if (isTRUE(computeBayesFactor)) {
    warning(
      "Package 'bridgesampling' not installed - skipping Bayes factor.",
      call. = FALSE
    )
  }

  result
}

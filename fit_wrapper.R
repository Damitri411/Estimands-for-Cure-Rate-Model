## fit_wrapper.R
##
## Thin wrappers around:
##   (a) survival::coxph()  -> the conventional, conflated population-level HR
##   (b) smcure::smcure()   -> the PH mixture-cure (PHMC) model used throughout
##                             the manuscript (Section 2.5.3: "Maximum likelihood
##                             estimation is implemented using standard software
##                             (e.g., R package smcure)")
##
## These wrappers translate smcure's internal sign convention into the
## manuscript's notation.
##
## SIGN CONVENTION NOTE (important):
## smcure's `cureform` model is fit internally as glm(w ~ Z, family=binomial),
## where w = E-step latent weight = P(SUSCEPTIBLE | data). So smcure's reported
## coefficients fit$b are for logit[P(susceptible)] = logit[1 - pi(z,x)].
## The manuscript defines incidence coefficients alpha via
## logit[pi(z,x)] = alpha0 + alpha1 z + alpha2'x, where pi = P(NON-susceptible).
## Since logit[1-pi] = -logit[pi], we have alpha = -fit$b.
## All functions below return coefficients already converted to the
## manuscript's alpha (pi = non-susceptible / "cure") convention.

suppressMessages({
  library(survival)
  library(smcure)
})

#' Fit the conventional (conflated) Cox proportional hazards model
#'
#' @return list(hr, log_hr, se_log_hr, ci_lower, ci_upper, fit)
fit_cox_conflated <- function(time, event, trt, covariates = NULL, data) {
  covs <- if (is.null(covariates)) "" else paste("+", paste(covariates, collapse = " + "))
  form <- as.formula(paste("Surv(time, event) ~ trt", covs))
  fit <- coxph(form, data = data)
  s <- summary(fit)
  list(
    hr        = unname(s$coefficients["trt", "exp(coef)"]),
    log_hr    = unname(s$coefficients["trt", "coef"]),
    se_log_hr = unname(s$coefficients["trt", "se(coef)"]),
    ci_lower  = unname(s$conf.int["trt", "lower .95"]),
    ci_upper  = unname(s$conf.int["trt", "upper .95"]),
    fit       = fit
  )
}

#' Fit the PH mixture-cure (PHMC) model via smcure and return manuscript-notation results
#'
#' Point-estimate-only fit (fast, ~0.1-0.3s). For CIs/coverage use
#' `phmc_bootstrap_ci()` below, which calls this function repeatedly on
#' resampled data (this is more flexible than smcure's own Var=TRUE path,
#' since it lets us obtain CIs for the risk difference and odds ratio, not
#' just the latency HR).
#'
#' @return list with: incidence coefficients (alpha, manuscript convention:
#'   logit[pi(z,x)], pi = P(non-susceptible)), latency Cox coefficients (beta),
#'   arm-specific cure rates at the covariate mean (pi1_hat, pi0_hat), the
#'   incidence risk difference (rd_hat = pi1_hat - pi0_hat), incidence odds
#'   ratio (or_hat), and latency hazard ratio (hr_hat).
fit_phmc_smcure <- function(time, event, trt, covariates = NULL, data,
                             emmax = 100, eps = 1e-7) {

  covs <- if (is.null(covariates)) "" else paste("+", paste(covariates, collapse = " + "))
  lat_form <- as.formula(paste("Surv(time, event) ~ trt", covs))
  inc_form <- as.formula(paste("~ trt", covs))

  d <- data
  d$time <- time
  d$event <- event
  d$trt <- trt

  invisible(capture.output(
    out <- smcure(lat_form, cureform = inc_form, data = d, model = "ph",
                  Var = FALSE, emmax = emmax, eps = eps)
  ))

  ## NOTE: smcure returns names(out$b) as R's auto-generated "Z[, -1]trt" etc.,
  ## but out$bnm / out$betanm hold the clean covariate names -- use those.
  alpha <- -out$b                      # convert to manuscript's pi (non-susceptible) convention
  names(alpha) <- out$bnm
  beta <- out$beta
  names(beta) <- out$betanm

  ## per-arm cure rate: pi(z, xbar) at the sample mean of the covariates,
  ## trt fixed at 0/1 (matches how Table 5 reports arm-specific cure rates)
  Z_names <- out$bnm[-1]  # drop intercept
  xbar <- sapply(Z_names, function(nm) if (nm == "trt") NA else mean(d[[nm]]))
  pi_at <- function(trt_val) {
    zz <- xbar; zz["trt"] <- trt_val
    zz[is.na(zz)] <- 0
    lp <- alpha["(Intercept)"] + sum(alpha[Z_names] * zz)
    plogis(lp)
  }
  pi1_hat <- pi_at(1)
  pi0_hat <- pi_at(0)
  or_hat  <- (pi1_hat / (1 - pi1_hat)) / (pi0_hat / (1 - pi0_hat))
  rd_hat  <- pi1_hat - pi0_hat
  hr_hat  <- exp(beta["trt"])

  list(
    alpha = alpha, beta = beta,
    pi1_hat = unname(pi1_hat), pi0_hat = unname(pi0_hat),
    rd_hat = unname(rd_hat), or_hat = unname(or_hat), hr_hat = unname(hr_hat),
    fit = out
  )
}

#' Percentile-bootstrap 95% CIs for the PHMC incidence risk difference,
#' incidence odds ratio, and latency hazard ratio (Section 2.5.4: "Standard
#' errors derive from the observed information matrix or bootstrap
#' resampling"). Stratified by event status, as in smcure's own bootstrap.
#'
#' @return list(point = fit_phmc_smcure(...) on original data,
#'              rd_ci, or_ci, hr_ci  [each a length-2 vector])
phmc_bootstrap_ci <- function(time, event, trt, covariates = NULL, data,
                               nboot = 50, emmax = 100, eps = 1e-7) {
  d <- data; d$time <- time; d$event <- event; d$trt <- trt
  point <- fit_phmc_smcure(time, event, trt, covariates, d, emmax, eps)

  idx1 <- which(d$event == 1); idx0 <- which(d$event == 0)
  n1 <- length(idx1); n0 <- length(idx0)

  rd_b <- or_b <- hr_b <- rep(NA_real_, nboot)
  for (b in seq_len(nboot)) {
    bs_idx <- c(sample(idx1, n1, replace = TRUE), sample(idx0, n0, replace = TRUE))
    db <- d[bs_idx, ]
    out_b <- tryCatch(
      fit_phmc_smcure(db$time, db$event, db$trt, covariates, db, emmax, eps),
      error = function(e) NULL
    )
    if (!is.null(out_b)) {
      rd_b[b] <- out_b$rd_hat
      or_b[b] <- out_b$or_hat
      hr_b[b] <- out_b$hr_hat
    }
  }
  qci <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))

  list(point = point,
       rd_ci = qci(rd_b), or_ci = qci(or_b), hr_ci = qci(hr_b),
       n_boot_ok = sum(!is.na(hr_b)))
}

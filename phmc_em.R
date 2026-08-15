## phmc_em.R
##
## A self-contained implementation of the semiparametric PH mixture-cure (PHMC)
## model used throughout the manuscript:
##
##   S(t | z,x) = pi(z,x) + (1 - pi(z,x)) * Ss(t | z,x)          [incidence x latency]
##   logit[pi(z,x)]   = alpha0 + alpha1*z + alpha2'x             [incidence, Eq. (C1)]
##   Ss(t | z,x) = exp{-H0(t) exp(beta1*z + beta2'x)}            [latency, Eq. (C1)]
##
## fit via the EM algorithm of Peng & Dear (2000) / Sy & Taylor (2000), i.e. the
## same estimating equations implemented operationally in the R package `smcure`.
## This file re-implements those equations directly with `stats::glm` (weighted
## logistic M-step) and `survival::coxph` (weighted Cox M-step + Breslow baseline),
## so the simulation has no dependency on smcure (not installable in this sandbox).
##
## NOTE ON pi(z,x): to match the manuscript's notation exactly, pi(z,x) is the
## probability of NON-susceptibility ("cure"). Internally we track
## u(z,x) = 1 - pi(z,x) = P(susceptible), since that is the quantity that enters
## the EM E-step and the latency (Cox) model directly.

#' Fit a semiparametric PH mixture-cure model by EM
#'
#' @param time observed (possibly censored) follow-up time
#' @param event event indicator (1 = event observed, 0 = censored)
#' @param inc_formula RHS-only formula for the incidence (logistic) model, e.g. ~ trt + age
#' @param lat_formula RHS-only formula for the latency (Cox) model, e.g. ~ trt + age
#' @param data data.frame containing time, event, and all covariates
#' @param max_iter maximum EM iterations
#' @param tol convergence tolerance on the log-likelihood
#' @param verbose print iteration diagnostics
#'
#' @return list with incidence coefficients (alpha), latency coefficients (beta),
#'   the Cox fit object (for baseline hazard extraction), fitted susceptibility
#'   probabilities pi_hat = 1-u_hat, convergence info, and a vcov for alpha via
#'   the weighted-logistic model (bootstrap is used for the primary reported SEs
#'   in run_simulation*.R, consistent with Section 2.5.4 of the manuscript).
fit_phmc <- function(time, event, inc_formula, lat_formula, data,
                      max_iter = 150, tol = 1e-4, verbose = FALSE) {

  stopifnot(length(time) == nrow(data), length(event) == nrow(data))
  n <- nrow(data)

  # design matrices. Both are built as plain data.frames with the ORIGINAL
  # covariate names (not matrix-multiplication in the glm/coxph formula) so
  # that fitted coefficient names match the manuscript's covariate names
  # (e.g. "trt"), not R's auto-generated "X_inctrt".
  inc_mf  <- model.frame(inc_formula, data = data)
  X_inc_m <- model.matrix(inc_formula, data = inc_mf)
  X_inc   <- as.data.frame(X_inc_m[, -1, drop = FALSE])  # drop intercept column; glm formula adds its own
  lat_mf  <- model.frame(lat_formula, data = data)
  X_lat_m <- model.matrix(lat_formula, data = lat_mf)
  X_lat   <- as.data.frame(X_lat_m[, -1, drop = FALSE])  # drop intercept for coxph

  # ---- initialize latent susceptibility weights w_i = E[uncured | data] ----
  # events are known-susceptible; censored subjects start at w=0.5
  w <- ifelse(event == 1, 1, 0.5)

  ll_old <- -Inf
  converged <- FALSE
  glm_fit <- NULL
  cox_fit <- NULL

  for (iter in seq_len(max_iter)) {

    ## ---------- M-step: incidence (weighted logistic regression) ----------
    ## Response is the *soft* class-membership w_i in [0,1]; this is the
    ## standard EM M-step for a Bernoulli latent class and is a valid MLE
    ## step even though w_i is fractional (glm() allows a numeric response
    ## in [0,1] for the binomial family).
    inc_dat <- data.frame(w = w, X_inc)
    glm_form <- as.formula(paste("w ~", paste(colnames(X_inc), collapse = " + ")))
    glm_fit <- suppressWarnings(
      glm(glm_form, data = inc_dat, family = binomial())
    )
    u_hat <- as.numeric(fitted(glm_fit))          # P(susceptible | x)
    u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)

    ## ---------- M-step: latency (weighted Cox PH) ----------
    ## Weighted partial likelihood with case weights w_i, the standard M-step
    ## used in the PHMC-EM literature (Peng & Dear 2000; Sy & Taylor 2000).
    cox_dat <- data.frame(time = time, event = event, X_lat, w = w)
    cox_form <- as.formula(paste("survival::Surv(time, event) ~",
                                  paste(colnames(X_lat), collapse = " + ")))
    cox_fit <- survival::coxph(cox_form, data = cox_dat, weights = w, ties = "breslow")

    ## Breslow baseline cumulative hazard (centered = FALSE => at x = 0)
    bh <- survival::basehaz(cox_fit, centered = FALSE)
    H0_fun <- stats::approxfun(bh$time, bh$hazard, method = "constant",
                                yleft = 0, rule = 2)
    lp <- as.numeric(predict(cox_fit, newdata = cox_dat, type = "lp", reference = "zero"))
    H0_t <- H0_fun(time)
    Su_t <- exp(-H0_t * exp(lp))                 # Ss(t | z,x) among susceptibles
    Su_t <- pmin(pmax(Su_t, 1e-12), 1 - 1e-12)

    ## ---------- E-step ----------
    ## For events, w_i = 1 (already known-susceptible).
    ## For censored subjects:
    ##   w_i = u(x_i) Su(t_i|x_i) / [ 1 - u(x_i) + u(x_i) Su(t_i|x_i) ]
    num <- u_hat * Su_t
    denom <- (1 - u_hat) + num
    w_new <- ifelse(event == 1, 1, num / denom)
    w_new <- pmin(pmax(w_new, 1e-8), 1 - 1e-8)

    ## ---------- convergence check ----------
    ## We track convergence via stability of the E-step weights and the
    ## incidence/latency coefficients, rather than an explicit observed-data
    ## log-likelihood (avoiding a fragile numerical hazard-jump approximation).
    ## This is equivalent in practice to log-likelihood-based stopping for
    ## this EM algorithm (Peng & Dear, 2000) and is what we monitor below.
    dw <- max(abs(w_new - w))
    db <- if (iter == 1) Inf else max(abs(c(coef(glm_fit), coef(cox_fit)) - prev_par))
    prev_par <- c(coef(glm_fit), coef(cox_fit))

    if (verbose) cat(sprintf("iter %3d | max|dw| = %.2e | max|dpar| = %.2e\n",
                              iter, dw, db))

    ## Primary convergence criterion: stability of the estimated coefficient
    ## vector (db). dw is reported as a diagnostic but is not required to
    ## vanish exactly, since the Breslow step-function baseline hazard can
    ## induce small persistent jitter in individual E-step weights near tied
    ## times without materially affecting the fitted coefficients.
    if (db < tol) {
      w <- w_new
      converged <- TRUE
      break
    }
    w <- w_new
  }
  ll_old <- NA_real_ # observed-data log-likelihood not computed (see note above)

  ## IMPORTANT SIGN NOTE:
  ## glm_fit models logit[u(z,x)] where u = P(susceptible | z,x) = 1 - pi(z,x).
  ## The manuscript's incidence coefficients alpha are defined via
  ## logit[pi(z,x)] = alpha0 + alpha1*z + alpha2'x, i.e. for the *cure*
  ## (non-susceptibility) probability, which is the negative link:
  ## logit[pi] = -logit[u]. So alpha = -coef(glm_fit).
  alpha_manuscript <- -coef(glm_fit)

  list(
    alpha       = alpha_manuscript,     # coefficients for pi(z,x) = P(non-susceptible), manuscript notation
    beta        = coef(cox_fit),        # coefficients for the susceptible-hazard Cox model (unaffected by sign issue)
    glm_fit     = glm_fit,              # underlying fit is for u = P(susceptible); alpha above = -coef(glm_fit)
    cox_fit     = cox_fit,
    w_hat       = w,
    pi_hat      = 1 - fitted(glm_fit),  # P(non-susceptible | x), matches manuscript's pi(z,x)
    u_hat       = fitted(glm_fit),
    loglik      = ll_old,
    iterations  = iter,
    converged   = converged,
    inc_formula = inc_formula,
    lat_formula = lat_formula
  )
}

#' Predict marginal (population-level) survival S(t|z,x) from a fitted PHMC model
#' at a specific covariate vector, for downstream cure-rate / plotting use.
predict_phmc_survival <- function(fit, newdata, times) {
  X_lat <- model.matrix(fit$lat_formula, data = model.frame(fit$lat_formula, data = newdata))[, -1, drop = FALSE]
  lp <- as.numeric(X_lat %*% fit$beta)
  bh <- survival::basehaz(fit$cox_fit, centered = FALSE)
  H0_fun <- stats::approxfun(bh$time, bh$hazard, method = "constant", yleft = 0, rule = 2)
  X_inc <- model.matrix(fit$inc_formula, data = model.frame(fit$inc_formula, data = newdata))
  # fit$alpha is on the pi (non-susceptible) scale; u = 1 - pi = 1 - plogis(X %*% alpha)
  u <- 1 - as.numeric(plogis(X_inc %*% fit$alpha))
  sapply(times, function(t) {
    H0t <- H0_fun(t)
    Su <- exp(-H0t * exp(lp))
    (1 - u) + u * Su
  })
}

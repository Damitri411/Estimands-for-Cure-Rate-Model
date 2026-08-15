## monotonicity_sensitivity.R
##
## Reproduces the sensitivity analysis of Supplement Section 1: bounds on the
## latency hazard ratio under departures from monotonicity
## P(Ci(0)=1, Ci(1)=0) = 0, governed by the violator fraction epsilon and the
## violators'-to-always-susceptible hazard ratio kappa (Eq. 1-2 of the
## Supplement). This is the script referenced in the manuscript and
## Supplement as "monotonicity_sensitivity.R" (reproduced here verbatim from
## the abridged version given in Supplement Section 2, with MAKE_PLOT support
## added).
##
## Set MAKE_PLOT=1 in the environment, or run interactively, to produce the
## band plot (Figure: latency HR bounds vs epsilon).

pert_factor <- function(eps, kappa) 1 / ((1 - eps) + eps * kappa)
latency_bound <- function(HR_obs, eps, kappa) HR_obs * pert_factor(eps, kappa)

trials <- list(
  "SURPASS-CVOT" = list(HR_obs = 1.003, pi_treated = 0.825),
  "REWIND"       = list(HR_obs = 1.126, pi_treated = 0.754)
)

eps_grid <- seq(0, 0.10, by = 0.005)
kappa_lo <- 0.5; kappa_hi <- 2.0; kappa_mid <- 1.5

## ---- Table 1 reproduction: bounds at eps = 0.10 ----
cat("Latency-HR bounds at eps = 0.10, kappa in [0.5, 2.0]\n")
cat(sprintf("%-13s %10s %10s %10s\n", "Trial", "HR_obs", "Lower", "Upper"))
for (nm in names(trials)) {
  HR <- trials[[nm]]$HR_obs
  lo <- latency_bound(HR, 0.10, kappa_hi)  # kappa>1 -> factor<1
  hi <- latency_bound(HR, 0.10, kappa_lo)  # kappa<1 -> factor>1
  cat(sprintf("%-13s %10.3f %10.3f %10.3f\n", nm, HR, lo, hi))
}

## ---- full bound curves across the epsilon grid (for CSV export / plotting) ----
build_bound_curve <- function(nm) {
  HR <- trials[[nm]]$HR_obs
  data.frame(
    trial = nm,
    eps = eps_grid,
    HR_obs = HR,
    lower_kappa2.0 = latency_bound(HR, eps_grid, kappa_hi),
    upper_kappa0.5 = latency_bound(HR, eps_grid, kappa_lo),
    mid_kappa1.5   = latency_bound(HR, eps_grid, kappa_mid)
  )
}
monotonicity_sensitivity_bounds <- do.call(rbind, lapply(names(trials), build_bound_curve))

## write the accompanying CSV referenced in the Supplement
if (!interactive() || Sys.getenv("WRITE_CSV", "0") == "1") {
  out_path <- if (dir.exists("outputs")) "outputs/monotonicity_sensitivity_bounds.csv" else "monotonicity_sensitivity_bounds.csv"
  write.csv(monotonicity_sensitivity_bounds, out_path, row.names = FALSE)
  cat("Wrote", out_path, "\n")
}

## ---- optional band plot ----
if (interactive() || Sys.getenv("MAKE_PLOT", "0") == "1") {
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  for (nm in names(trials)) {
    d <- monotonicity_sensitivity_bounds[monotonicity_sensitivity_bounds$trial == nm, ]
    plot(d$eps, d$mid_kappa1.5, type = "n",
         ylim = range(c(d$lower_kappa2.0, d$upper_kappa0.5)),
         xlab = expression(epsilon~"(violator fraction of treated-arm susceptible pool)"),
         ylab = "Latency HR bound",
         main = nm)
    polygon(c(d$eps, rev(d$eps)), c(d$lower_kappa2.0, rev(d$upper_kappa0.5)),
            col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
    lines(d$eps, d$lower_kappa2.0, col = "steelblue", lwd = 2, lty = 2)
    lines(d$eps, d$upper_kappa0.5, col = "steelblue", lwd = 2, lty = 2)
    lines(d$eps, d$mid_kappa1.5, col = "darkorange", lwd = 2)
    abline(h = d$HR_obs[1], col = "black", lty = 3)
    legend("topleft", legend = c(expression(kappa*"=1.5 (working value)"), "bounds (kappa in [0.5,2.0])", "HR_obs (eps=0)"),
           col = c("darkorange", "steelblue", "black"), lty = c(1, 2, 3), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
  }
  par(op)
}

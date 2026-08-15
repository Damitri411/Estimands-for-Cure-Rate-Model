## run_all.R
##
## Master driver for the PREDOSE-paper CVOT non-susceptibility simulation
## code. Sources all component files in the right order and runs the full
## suite. Adjust N_REP / NBOOT below to trade off fidelity to the manuscript
## (500 reps) against runtime.
##
## Requires the `smcure` package (installed from its CRAN mirror on GitHub;
## see install_smcure.R if not already installed) and `survival` (base R
## recommended package, normally pre-installed).
fp<-"C:/Users/Damitri/Desktop/Cure Rate Model"
suppressMessages({
  library(survival)
  library(smcure)
  library(parallel)
})

here <- function(...) file.path(fp, ...)
source(here("simulate_data.R"))
source(here("fit_wrapper.R"))
source(here("run_simulation_PHS.R"))
source(here("run_simulation_NPHS.R"))
source(here("run_simulation_table4.R"))
source(here("make_figure3.R"))

dir.create("outputs", showWarnings = FALSE)

## ---- configuration ----
## Set to the manuscript's values (500 / 40-100) for a full reproduction;
## the defaults below are lighter so a first run completes in minutes.
N_REP        <- 500   # replications per scenario (manuscript: 500)
N_PER_TRIAL  <- 1000  # subjects per simulated trial (manuscript: 1000, 500/arm)
N_REP_TABLE4 <- 100   # replications for the bootstrap-coverage table (slow path)
NBOOT        <- 40    # bootstrap draws per replicate for Table 4 coverage
SEED         <- 2026

## ============================================================
## 1. Primary PH-S simulation (Section 3.1 / Figure 3 / Table 4 bias columns)
## ============================================================
cat("\n=== 1. PH-S simulation (", N_REP, "reps x 13 scenarios) ===\n")
phs <- run_simulation_PHS(n_rep = N_REP, n_per_trial = N_PER_TRIAL, seed = SEED)
write.csv(phs$summary, "outputs/phs_simulation_summary.csv", row.names = FALSE)
write.csv(phs$raw, "outputs/phs_simulation_raw.csv", row.names = FALSE)

png("outputs/figure3_bias_vs_identifiability.png", width = 1100, height = 750, res = 130)
plot_figure3(phs)
dev.off()

## ============================================================
## 2. NPH-S robustness scenarios
## ============================================================
cat("\n=== 2. NPH-S robustness simulation ===\n")
nphs <- run_simulation_NPHS(n_rep = N_REP, n_per_trial = N_PER_TRIAL, seed = SEED + 1)
write.csv(nphs$summary, "outputs/nphs_simulation_summary.csv", row.names = FALSE)
write.csv(nphs$raw, "outputs/nphs_simulation_raw.csv", row.names = FALSE)

## ============================================================
## 3. Table 4: bias + 95% CI coverage for the 3 representative scenarios
## ============================================================
cat("\n=== 3. Table 4 bias/coverage simulation (bootstrap CIs; slow) ===\n")
tab4 <- run_simulation_table4(n_rep = N_REP_TABLE4, n_per_trial = N_PER_TRIAL,
                               nboot = NBOOT, seed = SEED + 2)
write.csv(tab4$summary, "outputs/table4_bias_coverage_summary.csv", row.names = FALSE)

## ============================================================
## 4. Monotonicity sensitivity analysis (Supplement S1)
## ============================================================
cat("\n=== 4. Monotonicity sensitivity analysis ===\n")
Sys.setenv(WRITE_CSV = "1", MAKE_PLOT = "0")
source(here("monotonicity_sensitivity.R"))  # writes outputs/monotonicity_sensitivity_bounds.csv directly

cat("\nAll done. See the outputs/ directory for CSV summaries and the Figure 3 plot.\n")

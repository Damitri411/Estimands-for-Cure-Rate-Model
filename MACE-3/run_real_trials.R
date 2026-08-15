## run_real_trials.R
##
## Driver script: generates Table 5, Table 6, Figure 1, and Figure 4 for
## SURPASS-CVOT and REWIND, calibrated to the real, published aggregate
## summary statistics of those trials (see R/real_trial_stats.R for sources
## and R/reconstruct_real_trials.R for the important caveats on what
## "calibrated" means here -- no individual patient-level trial data is
## used or available).
##
## Usage:
##   Edit WORKING_DIR below to the folder containing this script and its R/
##   subfolder, then source() this file.

WORKING_DIR <- "C:/Users/Damitri/Desktop/Cure Rate Model"
setwd(WORKING_DIR)

suppressMessages({
  library(survival)
  library(smcure)
})

here <- function(...) file.path("R", ...)
source(here("simulate_data.R"))
source(here("fit_wrapper.R"))
source(here("real_trial_stats.R"))
source(here("reconstruct_real_trials.R"))
source(here("generate_table5_table6.R"))
source(here("generate_figure1_figure4.R"))

dir.create("outputs", showWarnings = FALSE)

cat("Calibrating synthetic datasets and fitting Table 5 / Table 6 ...\n")
res <- run_table5_table6()

cat("\n==== TABLE 5 (Cox vs. PHMC, treatment-policy estimand) ====\n")
print(res$table5, row.names = FALSE)
write.csv(res$table5, "outputs/table5_real_trials.csv", row.names = FALSE)

cat("\n==== TABLE 6 (treatment-policy vs. while-on-treatment) ====\n")
print(res$table6, row.names = FALSE)
write.csv(res$table6, "outputs/table6_real_trials.csv", row.names = FALSE)

cat("\nGenerating Figure 1 (Kaplan-Meier curves) ...\n")
png("outputs/figure1_km_real_trials.png", width = 1300, height = 650, res = 130)
plot_figure1_km(res$calibrations)
dev.off()

cat("Generating Figure 4 (adjusted Cox HR forest plot) ...\n")
png("outputs/figure4_forest_real_trials.png", width = 1300, height = 650, res = 130)
plot_figure4_forest(res$calibrations)
dev.off()

cat("\nDone. See outputs/table5_real_trials.csv, outputs/table6_real_trials.csv,\n",
    "outputs/figure1_km_real_trials.png, outputs/figure4_forest_real_trials.png\n")

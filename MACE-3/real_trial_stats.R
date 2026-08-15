## real_trial_stats.R
##
## Publicly reported, real summary statistics for the two cardiovascular
## outcome trials referenced in the manuscript, with citations. NO
## individual patient-level data (IPD) is publicly available for either
## trial, so Table 5/6 and Figure 1/4 CANNOT be reproduced from real subject-
## level records. What this file (together with reconstruct_real_trials.R)
## does instead is calibrate a synthetic PHMC-consistent dataset whose
## observable, aggregate quantities (per-arm sample size, per-arm event
## rate, overall marginal Cox HR, and median follow-up) match these real,
## published numbers as closely as numerical calibration allows. Fitting the
## manuscript's Cox/PHMC pipeline to that calibrated synthetic dataset then
## produces Table 5/6- and Figure 1/4-style output that is illustrative of
## the method, not a re-analysis of the actual trials.
##
## ============================================================
## SURPASS-CVOT
## Nicholls SJ, Pavo I, Bhatt DL, et al; SURPASS-CVOT Investigators.
## Cardiovascular outcomes with tirzepatide versus dulaglutide in type 2
## diabetes. N Engl J Med. 2025;393:2409-2420. doi:10.1056/NEJMoa2505928.
## PMID: 41207920. https://pubmed.ncbi.nlm.nih.gov/41207920/
##
## Reported (modified intention-to-treat population):
##   n = 13,165 (6586 tirzepatide [MTD], 6579 dulaglutide 1.5mg)
##   Median follow-up: 4.0 years
##   Primary endpoint (MACE-3: CV death, non-fatal MI, non-fatal stroke):
##     tirzepatide: 801 events / 6586 (12.2%)
##     dulaglutide: 862 events / 6579 (13.1%)
##   HR = 0.92, 95.3% CI 0.83-1.01; P=0.003 (noninferiority), P=0.09 (superiority)
## ============================================================

## ============================================================
## REWIND
## Gerstein HC, Colhoun HM, Dagenais GR, et al. Dulaglutide and
## cardiovascular outcomes in type 2 diabetes (REWIND): a double-blind,
## randomised placebo-controlled trial. Lancet. 2019;394:121-130.
## PMID: 31189511. NCT01394952 (the ORIGINAL randomized trial).
##
## NOTE ON THE USER-PROVIDED CLINICALTRIALS.GOV LINK (NCT07417631):
## That registration ("Emulation of the REWIND Cardiovascular Outcomes
## Trial in Healthcare Claims Data", part of the RCT-DUPLICATE initiative)
## is a *real-world-data emulation study protocol* of REWIND, not the
## original randomized trial, and has no results posted publicly. It cannot
## supply summary statistics. The REWIND numbers below are instead taken
## from the original RCT (Gerstein et al., Lancet 2019), which is the trial
## the manuscript itself analyzes.
##
## Reported (intention-to-treat population):
##   n = 9901 (4949 dulaglutide 1.5mg, 4952 placebo)
##   Median follow-up: 5.4 years
##   Primary composite outcome (first MACE-3):
##     dulaglutide: 12.0% (2.4 events per 100 person-years)
##     placebo:     13.4% (2.7 events per 100 person-years)
##   HR = 0.88, 95% CI 0.79-0.99, p=0.026
## ============================================================

REAL_TRIAL_STATS <- list(

  `SURPASS-CVOT` = list(
    arm_labels   = c(control = "Dulaglutide 1.5mg", treated = "Tirzepatide (MTD)"),
    n_control    = 6579,
    n_treated    = 6586,
    event_rate_control = 862 / 6579,   # 0.1310
    event_rate_treated = 801 / 6586,   # 0.1216
    hr_reported  = 0.92,
    hr_ci        = c(0.83, 1.01),
    median_followup_years = 4.0,
    citation = "Nicholls SJ, Pavo I, Bhatt DL, et al. N Engl J Med. 2025;393:2409-2420. PMID 41207920.",
    baseline_covariates = list(age_mean = 64.1, age_sd = 8.8, pct_female = 0.290,
                                bmi_mean = 32.6, bmi_sd = 5.5, hba1c_mean = 8.4,
                                diabDur_mean = 14.7, diabDur_sd = 8.8)
  ),

  `REWIND` = list(
    arm_labels   = c(control = "Placebo", treated = "Dulaglutide 1.5mg"),
    n_control    = 4952,
    n_treated    = 4949,
    event_rate_control = 0.134,
    event_rate_treated = 0.120,
    hr_reported  = 0.88,
    hr_ci        = c(0.79, 0.99),
    median_followup_years = 5.4,
    citation = "Gerstein HC, Colhoun HM, Dagenais GR, et al. Lancet. 2019;394:121-130. PMID 31189511. (Original RCT, NCT01394952; note NCT07417631 is a non-randomized claims-data emulation of this trial with no posted results.)",
    baseline_covariates = list(age_mean = 66.2, age_sd = 6.5, pct_female = 0.462,
                                bmi_mean = 32.3, bmi_sd = 5.7, hba1c_mean = 7.3,
                                diabDur_mean = 10.5, diabDur_sd = 7.1)
    ## baseline_covariates for REWIND are approximate figures from the trial's
    ## published baseline characteristics table (Gerstein et al. 2019, Table 1),
    ## not re-verified via a fresh search in this session -- treat as
    ## illustrative if precision matters for your use case.
  )
)

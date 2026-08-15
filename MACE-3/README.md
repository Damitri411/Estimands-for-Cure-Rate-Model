# Estimand and Causal Framework for Cardiovascular Trials with Non-Susceptible Fractions — Simulation Code

This repository contains a self-contained R implementation of the simulation study, PH mixture-cure (PHMC) model fitting, and monotonicity sensitivity analysis described in the paper below. It reproduces the paper's methodology (Cox vs. PHMC comparison across an identifiability gradient, incidence/latency decomposition, non-proportional-hazards robustness checks, and the ε–κ monotonicity sensitivity bounds), using the same estimation software (`smcure`) referenced in the manuscript.

> **Estimand and Causal Framework for Cardiovascular Trials with Non-Susceptible Fractions**
>
> Cardiovascular outcome trials assume all randomized patients remain susceptible to future events, yet contemporary therapies may render a non-negligible fraction effectively non-susceptible. When this occurs, conventional estimands — hazard ratios and restricted mean survival times — conflate distinct treatment mechanisms, producing bias and limiting causal interpretability. We develop an estimand and causal inference framework for survival data with a latent non-susceptible fraction, aligned with ICH E9(R1) guidelines. Using mixture cure models and potential outcomes theory, treatment effects are decomposed into two components: the effect on long-horizon non-susceptibility (incidence) and the effect on event timing among susceptible individuals (latency), accommodating intercurrent events via strategy-specific counterfactuals. Simulation studies show that Cox proportional hazards estimands exhibit substantial bias in the presence of non-susceptible fractions, even under correct model specification. We additionally evaluate robustness to non-proportional hazards among susceptibles and quantify bias in both incidence and latency estimands across identifiability conditions. Non-susceptibility-aware estimands recover distinct causal mechanisms with superior interpretability. Application to two large cardiovascular outcome trials suggests that, to the extent an effect is present, it is concentrated in the incidence component — an increase in the probability of long-term event-free survival — with little evidence of an effect on event timing among susceptible individuals; however, both components are estimated imprecisely under the observed follow-up, so this decomposition is indicative rather than definitive. This distinction is entirely obscured by conventional hazard ratios. The proposed framework — with explicit formal identification conditions, explicitly layered by causal role and testability, clarified estimand targets, and an open-source interactive dashboard and Statistical Analysis Plan template — equips statisticians and regulators with practical tools to define and interpret treatment effects when sustained event-free survival is biologically plausible.

**Keywords:** Estimands; Non-susceptibility; Causal inference; Directed Acyclic Graphs (DAGs); Mixture Cure Models; Time-to-event analysis

---

## What's in this repository

```
.
├── run_all.R                      # master driver for the simulation study — run this
├── run_real_trials.R              # generates Table 5/6, Figure 1/4 for SURPASS-CVOT & REWIND
├── install_smcure.R               # installs the smcure package if needed
├── R/
│   ├── simulate_data.R            # data-generating process (PH-S and NPH-S)
│   ├── fit_wrapper.R              # wrappers around coxph() and smcure()
│   ├── run_simulation_PHS.R       # primary 13-scenario simulation (Section 3.1 / Fig. 3)
│   ├── run_simulation_NPHS.R      # non-proportional-hazards robustness scenarios
│   ├── run_simulation_table4.R    # bias + bootstrap 95% CI coverage (Table 4)
│   ├── make_figure3.R             # reproduces Figure 3 (bias vs. identifiability)
│   ├── monotonicity_sensitivity.R # ε–κ latency-HR bounds (Supplement S1, Table 1)
│   ├── real_trial_stats.R         # real, published SURPASS-CVOT/REWIND summary statistics
│   ├── reconstruct_real_trials.R  # calibrates synthetic IPD to match those statistics
│   ├── generate_table5_table6.R   # builds Table 5 and Table 6 from the calibrated data
│   ├── generate_figure1_figure4.R # builds Figure 1 (KM curves) and Figure 4 (forest plot)
│   └── optional_custom_em/
│       └── phmc_em.R              # hand-rolled PHMC EM fitter (didactic reference only)
└── outputs/                       # created on first run; all results land here
```

## What each script does

| Script | Reproduces | Output |
|---|---|---|
| `run_simulation_PHS.R` | The primary PH-S simulation: 13 scenarios spanning the Su(τ) identifiability gradient (Section 3.1), bracketing the manuscript's three named cases — well-identified (τ=8yr), REWIND-calibrated (τ=5.4yr), insufficient follow-up (τ=1.5yr) | `phs_simulation_summary.csv`, `phs_simulation_raw.csv` |
| `make_figure3.R` | Figure 3 — estimation bias (log scale) vs. Su(τ) for Cox HR, PHMC incidence OR, and PHMC latency HR | `figure3_bias_vs_identifiability.png` |
| `run_simulation_NPHS.R` | The NPH-S robustness scenarios (crossing hazards among susceptibles) | `nphs_simulation_summary.csv`, `nphs_simulation_raw.csv` |
| `run_simulation_table4.R` | Table 4 — bias and empirical 95% CI coverage (via percentile bootstrap) for Cox HR, PHMC incidence OR, and PHMC latency HR, at the three representative scenarios | `table4_bias_coverage_summary.csv` |
| `monotonicity_sensitivity.R` | Supplement Section 1 — ε–κ perturbation bounds on the latency HR under departures from monotonicity, reproducing Table 1 (SURPASS-CVOT / REWIND bounds at ε=0.10) | `monotonicity_sensitivity_bounds.csv` |
| `run_real_trials.R` | Table 5, Table 6, Figure 1, and Figure 4, calibrated to the **real, published** summary statistics of SURPASS-CVOT and REWIND (see below) | `table5_real_trials.csv`, `table6_real_trials.csv`, `figure1_km_real_trials.png`, `figure4_forest_real_trials.png` |

## Table 5/6 and Figure 1/4 for the real trials

`run_real_trials.R` applies the same Cox/PHMC pipeline to **SURPASS-CVOT** and **REWIND**, calibrated to those trials' real, published aggregate results rather than to a purely simulated scenario:

- **SURPASS-CVOT**: Nicholls SJ, Pavo I, Bhatt DL, et al. *N Engl J Med.* 2025;393:2409-2420. PMID [41207920](https://pubmed.ncbi.nlm.nih.gov/41207920/). n=13,165; HR 0.92 (95.3% CI 0.83–1.01); median follow-up 4.0 years.
- **REWIND**: Gerstein HC, Colhoun HM, Dagenais GR, et al. *Lancet.* 2019;394:121-130. PMID 31189511 (original RCT, NCT01394952). n=9901; HR 0.88 (95% CI 0.79–0.99); median follow-up 5.4 years.

**On the REWIND ClinicalTrials.gov link provided during development** (`NCT07417631`): that registration is *"Emulation of the REWIND Cardiovascular Outcomes Trial in Healthcare Claims Data"* — a non-randomized, claims-data emulation study protocol under the RCT-DUPLICATE initiative, not the original randomized trial, and it has no results posted. The REWIND numbers used here come from the original RCT instead, since that is the trial with real published results and the one the manuscript's own Table 5/6 analyzes.

**No individual patient-level data (IPD) is publicly available for either trial**, and none is used here. `reconstruct_real_trials.R` instead numerically calibrates a synthetic, mixture-cure-consistent dataset so that its *aggregate* properties (per-arm N, per-arm event rate, marginal Cox HR, median follow-up) closely match the real published numbers, then fits the manuscript's Cox/PHMC pipeline to that calibrated synthetic dataset. Table 5/6 and Figure 1/4 built this way are **illustrative reconstructions of the method applied to realistic targets, not a re-analysis of the actual trials** — in particular, the incidence/latency *split* of each trial's reported hazard ratio is not identified by the aggregate HR alone (that split is exactly what the manuscript's PHMC model is designed to estimate from IPD), so the split shown here reflects a documented starting assumption (`assumed_incidence_share` in `reconstruct_real_trials.R`) refined only enough to match the reported aggregate HR — treat the specific cure-rate and latency-HR numbers as illustrative, not as the trials' true decomposition. Figure 4 in particular reproduces the manuscript's forest-plot *display format*, not an extraction of either trial's actual published subgroup analysis (the source papers report subgroup HR forest plots, not a multivariable model's per-covariate coefficients).

## Requirements

- R ≥ 4.0
- Packages: `survival` (ships with base R), `smcure`, `parallel` (ships with base R)

### Installing `smcure`

`smcure` implements the semiparametric PH mixture-cure model used throughout the manuscript (Section 2.5.3). If `install.packages("smcure")` doesn't work in your environment (e.g. a restricted network), run:

```r
source("install_smcure.R")
```

which falls back to installing directly from smcure's official CRAN mirror on GitHub (`https://github.com/cran/smcure`).

## Quick start

1. Clone or download this repository.
2. Open `run_all.R` and set `WORKING_DIR` at the top to the folder containing this repo (with `run_all.R` and `R/` as siblings):

   ```r
   WORKING_DIR <- "path/to/this/repo"
   ```

3. Run it:

   ```r
   source("run_all.R")
   ```

All results are written to an `outputs/` folder created alongside the script.

### Runtime

The default `run_all.R` settings are lighter than the manuscript's design (500 replications) so a first run finishes in a reasonable time. To match the manuscript exactly, edit the configuration block near the top of `run_all.R`:

```r
N_REP        <- 500   # replications per scenario (manuscript default)
N_PER_TRIAL  <- 1000  # subjects per simulated trial (500/arm, manuscript default)
N_REP_TABLE4 <- 100   # replications for the bootstrap-coverage table
NBOOT        <- 40    # bootstrap draws per replicate for Table 4 coverage
```

`run_simulation_table4.R` is the slow path (each replicate itself requires `NBOOT` model refits to build a percentile-bootstrap CI); expect it to take substantially longer than the PH-S/NPH-S scans. Parallel execution via `parallel::mclapply` is used automatically on Linux/Mac (forking across `detectCores() - 1` cores); on Windows it silently falls back to sequential execution.

## What's exact vs. what's approximated

| Component | Status | Detail |
|---|---|---|
| Estimation software | **Exact** | Uses `smcure()` directly, as stated in Section 2.5.3 |
| PHMC model structure | **Exact** | Logistic incidence + Cox PH latency, matching Eq. (1) and Condition (C1) |
| Monotonicity sensitivity bounds (ε–κ) | **Exact** | `monotonicity_sensitivity.R` reproduces the Supplement's formula and Table 1 numbers exactly (SURPASS-CVOT / REWIND bounds at ε=0.10) |
| Covariate distributions | **Exact** | Age, diabetes duration, sex, BMI, SBP drawn from the distributions in Section 3.1 |
| DGP functional form | **Exact** | Logistic incidence model, Weibull-PH latency model, as specified |
| Regression coefficients (α, β) | **Approximated** | Manuscript doesn't report the exact values used to generate its published table; this code uses documented, calibrated coefficients instead |
| Su(τ) identifiability grid | **Calibrated to match** | Root-finding (`calibrate_gamma0`) is used to hit the manuscript's three named Su(τ) targets (≈0.02 / 0.80 / 0.95) exactly, with 10 interpolated scenarios between them |
| Cure-fraction baseline (~20%) | **Approximated** | Chosen to be plausible for a CVOT setting; not stated explicitly in the manuscript |
| smcure sign convention | **Corrected** | `smcure`'s `cureform` internally models P(susceptible); `fit_wrapper.R` negates coefficients to match the manuscript's π(z,x) = P(non-susceptible) convention |
| NPH-S crossing-hazard scenarios | **Illustrative** | Crossing times/magnitudes chosen to demonstrate the qualitative pattern (incidence robust, latency biased); manuscript doesn't report exact NPH-S parameters |
| Table 4 bootstrap coverage | **Reduced scale by default** | Uses fewer replications/bootstrap draws (100 / 40) than the manuscript's 500 for tractable runtime; edit `N_REP_TABLE4`/`NBOOT` in `run_all.R` to scale up |



- **Estimation software**: `smcure()` is used directly (via `R/fit_wrapper.R`), matching Section 2.5.3 ("Maximum likelihood estimation is implemented using standard software (e.g., R package `smcure`)"). The wrapper corrects a sign-convention subtlety: `smcure`'s `cureform` model internally targets P(susceptible), whereas the manuscript's incidence coefficients α are defined for π(z,x) = P(non-susceptible); `fit_wrapper.R` returns coefficients already converted to the manuscript's convention.
- **Data-generating process**: the manuscript specifies the functional form of the DGP and the covariate distributions (Section 3.1) but does not report the exact regression coefficients underlying its published numbers. `R/simulate_data.R` and `R/run_simulation_PHS.R` use documented, calibrated coefficients (via root-finding on the incidence/latency intercepts) chosen to hit the manuscript's headline identifiability targets — Su(τ) ≈ 0.02 / 0.80 / 0.95 at the three representative follow-up horizons. This reproduces the qualitative findings reported (Cox attenuation bias growing with the non-susceptible fraction; PHMC incidence-contrast robustness; PHMC latency-HR degradation as Su(τ) grows) rather than the manuscript's exact, unreported numeric table.
- **`R/optional_custom_em/phmc_em.R`** is a from-scratch EM implementation of the PHMC estimating equations (Peng & Dear, 2000; Sy & Taylor, 2000), kept for reference/teaching purposes. It is not used by the main pipeline — `smcure` is the primary, better-tested estimator.

## R Shiny App

https://surl.li/gmdfeu

## Citation

If you use this code, please cite the manuscript:

> Estimand and Causal Framework for Cardiovascular Trials with Non-Susceptible Fractions. [Author names / journal / year — update with final publication details.]

## License

Add a license of your choice (e.g. MIT) before publishing this repository publicly.

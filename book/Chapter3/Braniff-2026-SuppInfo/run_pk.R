# run_pk.R
# PK-mediated confounding: main results + supplement results/diagnostics.
#
# Reference simulation — DGP summary (main paper Results)
# Main results — RMSE/bias/ER plots/tables (main paper Results)
# Supplement results/diagnostics — positivity, cum. risk/survival, PK diagnostics plots (SI F)

source("simulation/functions.R")

# -- Dry-run mode: DRY_RUN=1 Rscript simulation/run_pk.R --
dry_run <- Sys.getenv("DRY_RUN") == "1"

# -- Parallelism: uses SLURM allocation if available, else sequential --
n_cores     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "0"))
max_workers <- as.integer(Sys.getenv("MAX_FURRR_WORKERS", unset = "5"))  # cap workers (memory-bound, not CPU-bound)
if (n_cores > 4) {
  plan(multicore, workers = min(n_cores - 1, max_workers))  # fork-based: shares reference via copy-on-write
} else {
  plan(sequential)
}

# =============================================================================
# Parameters
# =============================================================================

# DGP specification for PK-mediated confounding.
# Confounding acts through biomarker-driven clearance changes (B -> CL -> E),
# with no dose interruption confounding (gamma.B = 0 disables dose-mediated path).
# See SI A of the paper for full DGP description and parameter tables.

param <- c(
  # -- PK targets (used to derive K, V.targ, v0.lmu below) --
  hl.targ       = 1,          # target half-life (intervals)
  Cl.targ       = 10,         # target median clearance (a.u.)
  ss.hl         = 4,          # half-lives to steady state

  # -- Dose randomization --
  a.frac        = 0.5,        # fraction allocated to low dose
  a.low         = 10,         # low dose (a.u.)
  a.high        = 20,         # high dose (a.u.)

  # -- Biomarker dynamics (indirect response model) --
  b0.lmu        = 0,          # log-mean baseline biomarker (median = exp(0) = 1)
  b0.lsd        = 0.25,       # log-SD baseline biomarker IIV
  IE50.lmu      = log(5),     # log IE50: exposure for half-max biomarker suppression
  IE50.lsd      = 0,          # log-SD IE50 IIV (0 = fixed effect)
  kB.lmu        = log(0.15),  # log convergence rate to setpoint
  kB.lsd        = 0,          # log-SD kB IIV (0 = fixed effect)
  b.sd          = 0.2,        # proportional biomarker SD

  # -- Dose interruption model --
  gamma.0       = -4,         # log-odds intercept (~2% basal rate at median B)
  gamma.B       = 0,          # biomarker effect on dose interruption log-odds (0 = disabled, pk-med. only)

  # -- PK model (1-cpt IV bolus) --
  v0.lsd        = 0.1,        # volume log-SD IIV
  cl0.lsd       = 0.1,        # clearance log-SD IIV
  v0.cl0.lcor   = 0,          # V0-CL0 log-correlation
  phiB.lmu      = log(1),     # CL-biomarker effect, log scale (CL0 * 2^(-phiB * (B0.med.est - B)))
  phiB.lsd      = 0,          # log-SD phiB IIV (0 = fixed effect)
  c.sd          = 0.05,       # proportional concentration measurement error SD

  # -- Response model (per-interval hazard) --
  alpha.0       = -2,         # log-odds intercept
  alpha.B       = -4,         # direct biomarker effect on hazard (note B0.med-B parameterization: higher B → higher risk)
  alpha.E       = -0.5        # direct exposure effect on hazard
)

# Derive PK rate constant and volume from target half-life and clearance.
# V = CL / K ensures the 1-cpt model has the specified half-life.
param['K']      <- log(2) / param['hl.targ']
param['V.targ'] <- param['Cl.targ'] / param['K']
param['v0.lmu'] <- log(param['V.targ'])
param['cl0.lmu'] <- log(param['Cl.targ'])

# =============================================================================
# Shared settings
# =============================================================================

Nref        <- 25000                                          # reference population size
Nints       <- 10                                              # number of intervals per subject
Nexp           <- 5                                            # number of static exposure interventions
exposure.range  <- c(0.5, 2.5)                                  # exposure range for static interventions
exposure.levels <- seq(exposure.range[1], exposure.range[2], length.out = Nexp)  # static intervention levels
interv.list     <- c(NA, exposure.levels)                       # NA = natural course, rest = static interventions
ss.int      <- max(round(param['hl.targ'] * param['ss.hl']), 1)   # interval at which SS is reached
eval.E      <- c(1, 2)                                        # exposure levels for ATE evaluation (= A/CL0 at each dose)
master.seed <- 42                                              # base seed for reproducibility
outdir      <- "paper"                                          # figure output directory (flat, for submission)
datadir     <- file.path(outdir, "data")                        # rds/data output directory

dir.create(outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

# #############################################################################
#
#  Reference simulation
#
# #############################################################################

# Simulate the large reference population under natural course and each static
# intervention. Part 1 and Part 2 both sample from this shared reference to
# avoid re-running the DGP each rep.
set.seed(master.seed)
reference <- run_simulation(param, Nref, Nints, interv.list)

# DGP summary plot: concentration, exposure, fraction-at-risk, biomarker,
# clearance, and dose interruption rate over intervals, stratified by dose group.
# Paper: @fig-sims-pk (Results)
p.summary <- plot_dgp_summary(reference, param)
print(p.summary)
ggsave(file.path(outdir, "pk_dgp_summary.pdf"), p.summary, width = 12, height = 6)


# #############################################################################
#
#  Main Results
#
# #############################################################################

Nsims  <- 250            # number of independent simulated trial
Nids.sim    <- 200       # subjects per simulated trial
use.true.pk <- FALSE     # FALSE = real FOCEI fit per rep; TRUE = use true DGP exposure
if (dry_run) { Nsims <- 2; Nids.sim <- 50 }

plot.methods <- c(
  "Logistic (Early)", "Logistic (SS)",   "Logistic (TAE)",
  "Survival (Early)", "Survival (SS)",   "Survival (TAE)",
  "TV Survival",      "Adj. TV Survival", "G-formula"
)

# Main simulation loop: each sim samples a trial from the reference population,
# giving both the ground truth exposure simulations for the true potential ACERC
# and the observed natural course, which is used to fit a popPK model, impute
# per-interval exposure, then run all estimators.
# Each method's ER curve (risk at each exposure.levels level) is collected into a
# long-format data frame with columns (E, Risk, sim [replicate index], Method).
# try_fit() wraps each estimator so individual failures return NULL without
# crashing the loop. Paper: @fig-risk-pk, @tbl-mse-pk
er.curves <- future_map_dfr(seq_len(Nsims), .options = furrr_options(seed = NULL), function(sim_i) {
  message(sprintf("\n--- Sim rep %d / %d ---", sim_i, Nsims))
  rep.seed <- master.seed + sim_i
  set.seed(rep.seed)

  # Pipeline: sample trial -> fit PK -> impute exposure
  sample.i  <- sample_trial(reference, param, Nids.sim)
  pk.fit    <- fit_pk(sample.i$obs, param, use.true.pk)
  data.prep <- impute_exposure(pk.fit, sample.i$obs, Nints, ss.int)

  # Fit all estimators to imputed data
  results <- list(
    # Ground truth from the full reference (natural course + interventions)
    `Ground Truth`       = compute_ground_truth(sample.i$all, exposure.levels, interv.list),
    # Logistic regression: end-of-study binary outcome, B0-adjusted
    `Logistic (Early)`   = try_fit(fit_logistic(data.prep, "Early", c("B0"), exposure.levels, Nints)),
    `Logistic (SS)`      = try_fit(fit_logistic(data.prep, "SS",    c("B0"), exposure.levels, Nints)),
    `Logistic (TAE)`     = try_fit(fit_logistic(data.prep, "TAE",   c("B0"), exposure.levels, Nints)),
    # Pooled logistic regression (discrete-time survival), B0-adjusted
    `Survival (Early)`   = try_fit(fit_pooled_logistic(data.prep, "Early", c("B0"), exposure.levels, Nints, ss.int)),
    `Survival (SS)`      = try_fit(fit_pooled_logistic(data.prep, "SS",    c("B0"), exposure.levels, Nints, ss.int)),
    `Survival (TAE)`     = try_fit(fit_pooled_logistic(data.prep, "TAE",   c("B0"), exposure.levels, Nints, ss.int)),
    # Time-varying PLR: per-interval imputed exposure, no/with B0 adjustment
    `TV Survival`        = try_fit(fit_pooled_logistic(data.prep, "E.hat", NULL,    exposure.levels, Nints, ss.int)),
    `Adj. TV Survival`   = try_fit(fit_pooled_logistic(data.prep, "E.hat", c("B0"), exposure.levels, Nints, ss.int)),
    # G-formula: causal estimator (no bootstrap CIs in sim study for speed)
    `G-formula`          = try_fit(fit_gformula(data.prep, param, exposure.levels, Nints, Nids.sim,
                                                n.boot = 0, seed = rep.seed, parallel = FALSE))
  )

  # Collect each method's risk vector into long format
  map_dfr(names(results), \(m) {
    if (is.null(results[[m]])) return(tibble(E = NA_real_, Risk = NA_real_, sim = sim_i, Method = m))
    results[[m]]$risk %>% mutate(sim = sim_i, Method = m)
  })
})

## Report fit failures -------------------------------------------------------
# Count methods that returned NULL (via try_fit) across all reps
n.fail <- er.curves %>% filter(is.na(Risk)) %>% count(Method, name = "failures")
if (nrow(n.fail) > 0) {
  message("\nFit failures by method:")
  print(n.fail)
} else {
  message("\nNo fit failures.")
}

## Summarize and plot ----------------------------------------------------------

# Compute median/quantile bands, pointwise RMSE/bias at eval.E, and ATE RMSE/bias
curve.stats <- process_er_curves(er.curves, eval.E)
# Format RMSE table: percentage-scaled, sorted by ATE RMSE. Paper: @tbl-mse-pk
mse.tbl     <- format_mse_table(curve.stats, eval.E, plot.methods)

# Faceted ER curves: mean + 95% band per method vs ground truth. Paper: @fig-risk-pk
p.er <- plot_er_sim_study(curve.stats, plot.methods, exposure.levels)
print(p.er)
ggsave(file.path(outdir, "pk_er_sim_study.pdf"), p.er, width = 12, height = 10)

# Save raw and summarized results for downstream use
saveRDS(er.curves,   file.path(datadir, "pk_er_curves_raw.rds"))
saveRDS(curve.stats, file.path(datadir, "pk_er_curves_stats.rds"))
saveRDS(mse.tbl,     file.path(datadir, "pk_mse_tbl.rds"))


# #############################################################################
#
#  Supplement Results/diagnostics
#
# #############################################################################

## --- PK diagnostics plots --------------
# Sample a trial cohort, fit the popPK model, and impute exposure at N=200.
set.seed(master.seed)
samp.diag      <- sample_trial(reference, param, Nids.sim)
pk.fit.diag    <- fit_pk(samp.diag$obs, param, use.true.pk = FALSE)
data.prep.diag <- impute_exposure(pk.fit.diag, samp.diag$obs, Nints, ss.int)

# Plots for IWRES, parameter accuracy, individual fits
# Paper: @fig-pk-resid-pk, @fig-pk-params-pk, @fig-pk-indfits-pk (SI F)
p.pk <- plot_pk_diagnostics(pk.fit.diag, prep_pk_data(samp.diag$obs), data.prep.diag)
ggsave(file.path(outdir, "pk_diag_resid.pdf"),   p.pk$residuals, width = 12, height = 4)
ggsave(file.path(outdir, "pk_diag_params.pdf"),  p.pk$params,    width = 10, height = 4)
ggsave(file.path(outdir, "pk_diag_indfits.pdf"), p.pk$indfits,   width = 12, height = 8)

## --- Survival and cumulative risk curves (N = 1000, large sample) ---------
Nids.curves <- 1000                    # large sample size to show method bias
n.boot      <- 100                     # bootstrap replicates for CIs
if (dry_run) { Nids.curves <- 50; n.boot <- 5 }  # shrink for fast smoke test
# Sample a trial cohort, fit the popPK model, and impute exposure at N=1000.
set.seed(master.seed)
samp.curves      <- sample_trial(reference, param, Nids.curves)
pk.fit.curves    <- fit_pk(samp.curves$obs, param, use.true.pk = FALSE)
data.prep.curves <- impute_exposure(pk.fit.curves, samp.curves$obs, Nints, ss.int)

# Fit all estimators and the g-formula to the large-N imputed data
results <- list(
  `Ground Truth`     = compute_ground_truth(samp.curves$all, exposure.levels, interv.list),
  `Logistic (Early)` = fit_logistic(data.prep.curves, "Early", c("B0"), exposure.levels, Nints, ci = TRUE, n.boot = n.boot),
  `Logistic (SS)`    = fit_logistic(data.prep.curves, "SS",    c("B0"), exposure.levels, Nints, ci = TRUE, n.boot = n.boot),
  `Logistic (TAE)`   = fit_logistic(data.prep.curves, "TAE",   c("B0"), exposure.levels, Nints, ci = TRUE, n.boot = n.boot),
  `Survival (Early)` = fit_pooled_logistic(data.prep.curves, "Early", c("B0"), exposure.levels, Nints, ss.int, ci = TRUE, n.boot = n.boot),
  `Survival (SS)`    = fit_pooled_logistic(data.prep.curves, "SS",    c("B0"), exposure.levels, Nints, ss.int, ci = TRUE, n.boot = n.boot),
  `Survival (TAE)`   = fit_pooled_logistic(data.prep.curves, "TAE",   c("B0"), exposure.levels, Nints, ss.int, ci = TRUE, n.boot = n.boot),
  `TV Survival`      = fit_pooled_logistic(data.prep.curves, "E.hat", NULL,    exposure.levels, Nints, ss.int, ci = TRUE, n.boot = n.boot),
  `Adj. TV Survival` = fit_pooled_logistic(data.prep.curves, "E.hat", c("B0"), exposure.levels, Nints, ss.int, ci = TRUE, n.boot = n.boot),
  `G-formula`        = fit_gformula(data.prep.curves, param, exposure.levels, Nints, Nids.curves,
                                    ci = TRUE, n.boot = n.boot, seed = master.seed)
)

# Paper: @fig-cumrisk-pk, @fig-surv-pk (SI F)
p.comp <- plot_summary_outcomes(results, exposure.levels, plot.methods)
ggsave(file.path(outdir, "pk_survival.pdf"),   p.comp$surv_step,  width = 5,  height = 12)
ggsave(file.path(outdir, "pk_cumrisk.pdf"),    p.comp$risk_curve, width = 5,  height = 12)

## --- Exposure positivity check (full reference population) ------
# Uses the natural-course from reference population to show the
# population-level positivity structure. Paper: @fig-pos-pk (SI F)
ref.obs <- censor_natural_course(reference, Rep)
p.pos <- plot_positivity(ref.obs, exposure.levels)
ggsave(file.path(outdir, "pk_positivity.pdf"), p.pos, width = 6, height = 10)


# functions.R
# Simulation and analysis functions for:
#   "Causal Exposure Response Analysis: Estimands, Time-varying
#    Confounding, and the Parametric g-Formula" (Braniff)
#
# Code Blocks:
#   1. Utilities & helpers    — try_fit, haz_to_surv, extract_risk, quantile_bands,
#                               bootstrap_surv_ci, glm_fit_predict
#   2. Model factories        — make_pk_model, make_bm_fit
#   3. DGP                    — rep_cols, init_state, sim_step, state_to_long,
#                               run_simulation
#   4. Sampling               — sample_trial
#   5. Ground truth           — compute_ground_truth
#   6. PK imputation          — prep_pk_data, fit_pk, impute_exposure
#   7. Estimators             — fit_logistic, fit_pooled_logistic, fit_gformula
#   8. Results               —  plot_dgp_summary, process_er_curves,
#                               plot_er_sim_study, format_mse_table
#   9. Supplemental          —  plot_positivity,
#                               plot_pk_diagnostics, plot_summary_outcomes
#
# Called by: run_dose.R, run_pk.R

library(purrr)
library(furrr)
library(gridExtra)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(data.table)
library(nlmixr2)
library(gfoRmula)
mvrnorm <- MASS::mvrnorm

# =============================================================================
# Plot defaults
# =============================================================================

theme_set(theme_bw())

METHOD_COLORS <- c(
  'Ground Truth'       = 'black',
  'G-formula'          = '#1a7a1a',
  'Logistic (Early)'   = '#d95f02',
  'Logistic (SS)'      = '#b03a00',
  'Logistic (TAE)'     = '#7f2704',
  'Survival (Early)'   = '#3182bd',
  'Survival (SS)'      = '#1c5a9c',
  'Survival (TAE)'     = '#08306b',
  'TV Survival'        = '#7b6daa',
  'Adj. TV Survival'   = '#54278f'
)

# =============================================================================
# Utilities
# =============================================================================

# Wrap a model fit call; returns NULL on failure
try_fit <- function(expr) tryCatch(expr, error = function(e) { warning(e$message); NULL })

# Restrict to natural-course rows, censored at first event (R == 1), per group
censor_natural_course <- function(df, ...) {
  df %>%
    filter(is.na(Intervention)) %>%
    group_by(...) %>%
    arrange(Ts, .by_group = TRUE) %>%
    slice(1:match(1, R, nomatch = n())) %>%
    ungroup()
}

# Cores allocated to this job. Prefers SLURM's cgroup-aware count over
# parallel::detectCores(), which can over-report the node's full core count
# on some cluster configs rather than just what's allocated to this job.
available_cores <- function() {
  slurm.cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)))
  if (!is.na(slurm.cores)) slurm.cores else parallel::detectCores()
}

# =============================================================================
# Hazard-Survival-Risk & GLM helpers
# =============================================================================

# Compute discrete-time survival and risk from interval hazards
haz_to_surv <- function(haz.df, ...) {
  haz.df %>%
    group_by(...) %>%
    arrange(Ts, .by_group = TRUE) %>%
    mutate(Surv = cumprod(1 - Haz),
           Risk = 1 - Surv) %>%
    ungroup()
}

# Get end-of-follow-up risk per exposure
extract_risk <- function(surv) {
  surv %>%
    group_by(E) %>%
    filter(Te == max(Te)) %>%
    ungroup() %>%
    select(E, Risk, any_of(c("Risk.low", "Risk.up")))
}

# 95% pointwise CIs from bootstrap/simulation replicates
quantile_bands <- function(surv.reps) {
  surv.reps %>%
    group_by(E, Ts, Te) %>%
    summarize(
      Haz.low  = quantile(Haz,  0.025, na.rm = TRUE),
      Haz.up   = quantile(Haz,  0.975, na.rm = TRUE),
      Surv.low = quantile(Surv, 0.025, na.rm = TRUE),
      Surv.up  = quantile(Surv, 0.975, na.rm = TRUE),
      Risk.low = quantile(Risk, 0.025, na.rm = TRUE),
      Risk.up  = quantile(Risk, 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
}

# Subject-level nonparametric bootstrap for survival/risk CI bands.
# Resamples IDs to preserve within-subject correlation in person-period data.
# predict_fn(data) → tibble(E, Ts, Te, Haz); typically wraps glm_fit_predict.
bootstrap_surv_ci <- function(data, n.boot, predict_fn, id_col = "ID") {
  ids <- unique(data[[id_col]])
  boot.haz <- map_dfr(seq_len(n.boot), function(b) {
    # Resample subjects with replacement; re-ID so duplicates are distinct
    boot.ids  <- sample(ids, replace = TRUE)
    boot.data <- tibble(..id.. = boot.ids, new_ID = seq_along(boot.ids)) %>%
      left_join(data, by = setNames(id_col, "..id.."), relationship = "many-to-many") %>%
      mutate(!!id_col := new_ID) %>%
      select(-..id.., -new_ID)
    result <- tryCatch(predict_fn(boot.data), error = function(e) NULL)
    if (!is.null(result)) mutate(result, boot_id = b) else NULL
  })
  boot.surv <- haz_to_surv(boot.haz, E, boot_id)
  quantile_bands(boot.surv)
}

# Fit logistic GLM and predict hazards on a grid; used by fit_logistic and fit_pooled_logistic.
# mask: logical vector — sets hazard to 0 for masked intervals (e.g. pre-steady-state).
glm_fit_predict <- function(fml, d, pg, mask = NULL) {
  fit <- glm(fml, data = d, family = binomial())
  pred.link <- predict(fit, newdata = pg)
  if (!is.null(mask)) pred.link[mask] <- -Inf # zero hazard for masked intervals (e.g. pre-SS)
  list(haz = plogis(pred.link), fit = fit)
}

# =============================================================================
# Model factories 
# =============================================================================

# Build nlmixr2 model function for 1-cmt PK with biomarker-dependent CL.
# sprintf is used because nlmixr2 can't capture variables from the enclosing function scope.
# CL model: CL = CL0 * 2^(-phiB * (B0.med.est - B)): lower B → lower CL.
# When no CL-B effect; lphiB is initialized to -99 because nlmixr2 ini() can't parse -Inf.
make_pk_model <- function(param) {
  eval(parse(text = sprintf(
    'function() {
       ini({
         lCL0        <- %s       ; label("log Cl (vol./cycle)")
         lV0         <- %s       ; label("log Vc (vol.")
         lphiB       <- %s       ; label("CL-B effect")
         eta.V0 + eta.CL0 ~ c(%s, 0, %s)
         prop.pk.err <- %s
       })
       model({
         V    <- exp(lV0 + eta.V0)
         phiB <- exp(lphiB)
         CL0  <- exp(lCL0 + eta.CL0)
         CL   <- CL0 * 2^(-phiB * (B0.med.est - B))
         C   <- linCmt()
         C ~ prop(prop.pk.err)
       })
     }',
    deparse(unname(param["cl0.lmu"])),
    deparse(unname(param["v0.lmu"])),
    deparse(max(unname(param["phiB.lmu"]), -99)),
    deparse(unname(param["v0.lsd"])^2),
    deparse(unname(param["cl0.lsd"])^2),
    deparse(unname(param["c.sd"]))
  )))
}

# Build gnls fit function for biomarker dynamics, used as covfits_custom in gfoRmula.
# fix.IE50 = FALSE: estimate both kB and IE50, default and used for simulation study main results.
# fix.IE50 = TRUE: IE50 hardcoded in formula. Used for SI scenarios where IE50 = 1000; unidentifiable).
make_bm_fit <- function(lkB.start, lIE50.start, fix.IE50 = FALSE) {
  # When fixed, substitute the literal IE50 value into the formula text instead of a
  # parameter name, so gnls never declares it as something to optimize.
  IE50.expr <- if (fix.IE50) sprintf('%.10f', lIE50.start) else 'lIE50.est'
  fml <- sprintf('B ~ (1 - exp(-exp(lkB.est))) * B0 / (1 + lag1_E / exp(%s)) + exp(-exp(lkB.est)) * lag1_B', IE50.expr)

  # params/start declare only the parameters gnls should actually estimate.
  if (fix.IE50) {
    pars  <- 'lkB.est ~ 1'
    start <- sprintf('list(lkB.est = %.10f)', lkB.start)
  } else {
    pars  <- 'lkB.est + lIE50.est ~ 1'
    start <- sprintf('list(lkB.est = %.10f, lIE50.est = %.10f)', lkB.start, lIE50.start)
  }

  # Build the covfits_custom biomarker fitting call for the gfoRmula.
  eval(parse(text = sprintf(
    'function(covparams, covname, obs_data, j) {
       fit <- nlme::gnls(%s,
         data = obs_data, params = %s,
         start = %s,
         weights = nlme::varPower(value = 1, fixed = 1),
         control = nlme::gnlsControl(minScale = 1e-10))
       fit
     }', fml, pars, start)))
}

# =============================================================================
# DGP simulation
# =============================================================================

rep_cols <- function(x, n_cols) matrix(rep(x, n_cols), nrow = length(x), ncol = n_cols)

# Draw per-individual parameters from DGP priors and initialize pre-trial state (interval 0).
# This state seeds sim_step at interval 1 but is not itself observed; run_simulation drops it.
# Returns list of Nref × Ninterv matrices (rows = individuals, cols = treatments)
init_state <- function(Nref, Ninterv, param) {
  # --- Biomarker parameters ---
  B0   <- rlnorm(Nref, param['b0.lmu'],   param['b0.lsd'])
  IE50 <- rlnorm(Nref, param['IE50.lmu'], param['IE50.lsd'])
  kB   <- rlnorm(Nref, param['kB.lmu'],   param['kB.lsd'])

  # --- PK parameters ---
  sd.V0  <- param['v0.lsd']
  sd.CL0 <- param['cl0.lsd']
  Sigma  <- matrix(c(sd.V0^2,
                     sd.V0 * sd.CL0 * param['v0.cl0.lcor'],
                     sd.V0 * sd.CL0 * param['v0.cl0.lcor'],
                     sd.CL0^2), nrow = 2)
  lV0.CL0 <- mvrnorm(Nref, mu = c(param['v0.lmu'], param['cl0.lmu']), Sigma = Sigma)
  V0   <- exp(lV0.CL0[, 1])
  CL0  <- exp(lV0.CL0[, 2])
  phiB <- rlnorm(Nref, param['phiB.lmu'], param['phiB.lsd'])

  # --- Dose assignment ---
  n.low <- round(Nref * param['a.frac'])
  A <- sample(c(rep(param['a.low'],  n.low),
                rep(param['a.high'], Nref - n.low)))

  list(
    Ts   = 0L,  Te   = 0L,             # time interval bounds
    # Per-subject random ICs broadcast across all Ninterv columns via rep_cols
    B0   = rep_cols(B0,   Ninterv),     # baseline biomarker
    IE50 = rep_cols(IE50, Ninterv),     # exposure IE50 for biomarker setpoint
    kB   = rep_cols(kB,   Ninterv),     # biomarker recovery rate
    V0   = rep_cols(V0,   Ninterv),     # volume of distribution
    CL0  = rep_cols(CL0,  Ninterv),     # baseline clearance
    phiB = rep_cols(phiB, Ninterv),     # CL–biomarker effect
    A    = rep_cols(A,    Ninterv),     # assigned dose
    B      = rep_cols(B0,   Ninterv),   # current biomarker
    D      = rep_cols(A,    Ninterv),   # administered dose
    CL     = rep_cols(CL0,  Ninterv),   # current clearance
    # Time-varying trajectories filled in by sim_step; zeroed placeholders here.
    Pi     = matrix(0,  nrow = Nref, ncol = Ninterv),  # P(interruption)
    I      = matrix(0L, nrow = Nref, ncol = Ninterv),  # interruption indicator
    Cmu0   = matrix(0,  nrow = Nref, ncol = Ninterv),  # true conc at interval start
    Cmu1   = matrix(0,  nrow = Nref, ncol = Ninterv),  # true conc at interval mid
    Cmu2   = matrix(0,  nrow = Nref, ncol = Ninterv),  # true conc at interval end
    C0     = matrix(0,  nrow = Nref, ncol = Ninterv),  # observed conc at interval start
    C1     = matrix(0,  nrow = Nref, ncol = Ninterv),  # observed conc at interval mid
    C2     = matrix(0,  nrow = Nref, ncol = Ninterv),  # observed conc at interval end
    E      = matrix(0,  nrow = Nref, ncol = Ninterv),  # exposure (Cavg)
    Pr     = matrix(0,  nrow = Nref, ncol = Ninterv),  # P(event)
    R      = matrix(0L, nrow = Nref, ncol = Ninterv)   # Response event indicator
  )
}

# Advance DGP one interval (Appendix A). Order: biomarker → dose → PK → outcome.
# Interval 1 bypasses biomarker dynamics and interruptions (first dose, baseline B).
sim_step <- function(state, t, interv.list, param) {
  Nref   <- nrow(state$B0)
  Ninterv <- ncol(state$B0)
  N.cells <- Nref * Ninterv
  B0.med <- exp(param['b0.lmu'])  # population median B0, centering reference

  # --- Biomarker dynamics (P_B) ---
  # Discrete autoregressive model: tracks Emax-type setpoint B0/(1+E/IE50),
  # first-order convergence at rate kB, with proportional process noise.
  if (t == 1) {      #t == 1: B in first interval starts at baseline.
    B <- state$B0
  } else {           #t > 1: B evolves dynamically according to model.
    B.setpoint <- state$B0 / (1 + state$E / state$IE50)
    # Δt = 1 per interval, so exp(-kB·Δt) reduces to exp(-kB) here; equivalent form to SI A.
    mu.B       <- B.setpoint - (B.setpoint - state$B) * exp(-state$kB)
    B          <- matrix(pmax(rnorm(N.cells, mu.B, mu.B * param['b.sd']), 0.001),
                         nrow = Nref, ncol = Ninterv)
  }

  # --- Dose interruption ---
  # Suppressed B → higher P(interruption), analogous to haematologic toxicity.
  # I=1 → D=0 this interval; resumes if I=0 next interval.
  # At interval 1: no interruptions possible — assigned dose given regardless of B.
  if (t == 1) {
    Pi     <- matrix(0,  nrow = Nref, ncol = Ninterv)
    I      <- matrix(0L, nrow = Nref, ncol = Ninterv)
    D      <- state$A
  } else {
    Pi     <- 1 / (1 + exp(-(param['gamma.0'] + param['gamma.B'] * (B0.med - B))))
    I      <- matrix(rbinom(N.cells, 1, Pi), nrow = Nref, ncol = Ninterv)
    D      <- (1-I)*state$A
  }

  # --- PK: 1-compartment IV bolus ---
  # CL modulated by biomarker deviation from median (PK-mediated confounding).
  # Three true concentrations per interval (start, mid, end); C observed with prop. error.
  CL   <- state$CL0 * 2^(-(B0.med - B) * state$phiB)
  Cmu0 <- state$Cmu2 + D / state$V0          # bolus on top of residual conc.
  #C(t+δ) = C(t)*exp(-CL/V0 * δ); δ = Δt/2 = 0.5 as is half of the unit interval
  Cmu1 <- Cmu0 * exp(-0.5 * CL / state$V0)   # mid-interval
  Cmu2 <- Cmu1 * exp(-0.5 * CL / state$V0)   # end-of-interval
  C0   <- matrix(rnorm(N.cells, Cmu0, param['c.sd'] * Cmu0), nrow = Nref, ncol = Ninterv)
  C1   <- matrix(rnorm(N.cells, Cmu1, param['c.sd'] * Cmu1), nrow = Nref, ncol = Ninterv)
  C2   <- matrix(rnorm(N.cells, Cmu2, param['c.sd'] * Cmu2), nrow = Nref, ncol = Ninterv)

  # --- Exposure: per-interval average concentration (Cavg) ---
  # interv.list determines the E column structure:
  #    NA = natural course, exposure from simulated PK (Cavg = AUC/1)
  #    numeric = intervention,  exposure fixed at specified counterfactual value
  E       <- matrix(rep(interv.list, each = Nref), nrow = Nref, ncol = Ninterv)
  nc.cols <- which(is.na(interv.list))  # natural-course columns
  # E = Cavg = (Cmu0/Δt) * (V0/CL) * (1 - exp(-CL/V0*Δt)); Δt = 1 so dropped
  E[, nc.cols] <- Cmu0[, nc.cols] * state$V0[, nc.cols] / CL[, nc.cols] *
                   (1 - exp(-CL[, nc.cols] / state$V0[, nc.cols]))  

  # --- Response model (P_R) ---
  # E_ref = population-average idealized steady-state exposure (weighted-mean Ess=dose/CL).
  E_ref <- (param['a.frac'] * param['a.low'] +
            (1 - param['a.frac']) * param['a.high']) /
           exp(param['cl0.lmu'])
  Pr <- 1 / (1 + exp(-(param['alpha.0'] +
                       param['alpha.E'] * (E - E_ref) +
                       param['alpha.B'] * (B0.med - B))))
  R <- matrix(rbinom(N.cells, 1, Pr), nrow = Nref, ncol = Ninterv)

  list(
    Ts = t - 1,  Te = t,
    B0 = state$B0,  IE50 = state$IE50,  kB = state$kB,
    V0 = state$V0,  CL0 = state$CL0,  phiB = state$phiB,  A = state$A,
    B = B,  Pi = Pi,  I = I,  D = D,
    CL = CL,  Cmu0 = Cmu0,  Cmu1 = Cmu1,  Cmu2 = Cmu2,
    C0 = C0,  C1 = C1,  C2 = C2,
    E = E,  Pr = Pr,  R = R
  )
}

# Flatten one time-step's state (list of Nref × Ninterv matrices) to long-format data.frame.
# One row per (individual, treatment); Ts/Te are scalars, everything else is unwound column-major.
state_to_long <- function(state, Nref, interv.list) {
  Ninterv <- length(interv.list)
  df <- data.frame(
    Intervention = rep(interv.list,     each  = Nref),  # exposure value; NA = natural course
    Rep          = rep(seq_len(Nref),   times = Ninterv),
    Ts           = state$Ts,
    Te           = state$Te
  )
  # as.vector unwinds column-major: all reps for trt 1, then trt 2, etc.
  for (v in setdiff(names(state), c('Ts', 'Te'))) df[[v]] <- as.vector(state[[v]])
  df
}

# Run full DGP: init → Nints steps → bind into one long data.frame.
# accumulate returns [init, step1, ..., stepN]; [-1] drops the pre-trial init (interval 0).
# Each step's state is flattened to long format, then row-bound across intervals.
run_simulation <- function(param, Nref, Nints, interv.list) {
  Ninterv <- length(interv.list)
  accumulate(
    seq_len(Nints),
    \(s, t) sim_step(s, t, interv.list, param),   # iterate the DGP sim_step over intervals
    .init = init_state(Nref, Ninterv, param)    # create and use pre-trial (interval 0) IC state
  )[-1] %>%                                    # [-1] drops pre-trial (interval 0) IC state
    map(state_to_long, Nref = Nref, interv.list = interv.list) %>%  # convert to long-format before bind
    bind_rows()
}

# =============================================================================
# Sampling
# =============================================================================

# Sample a trial from the reference simulation. Nids subjects are drawn
# (a.frac to low dose, remainder to high), then split into:
#   all — full counterfactual simulations (all interventions), for ground truth
#   obs — natural-course only, censored at first event; observed, for estimation
sample_trial <- function(reference, param, Nids) {
  n.low  <- round(Nids * param['a.frac'])
  sampled.reps <- c(
    sample(unique(reference$Rep[reference$A == param['a.low']]),  n.low),
    sample(unique(reference$Rep[reference$A == param['a.high']]), Nids - n.low)
  )

  #for computing within-sample ground truth
  all.intervene <- reference %>%
    filter(Rep %in% sampled.reps) %>%
    group_by(Rep) %>% mutate(ID = cur_group_id()) %>% ungroup()

  #for generating the observed data
  obs <- all.intervene %>%
    censor_natural_course(ID) %>%
    select(-Intervention, -Rep) %>%
    relocate(ID) %>%
    ungroup()

  list(all = all.intervene, obs = obs)
}

# =============================================================================
# Ground truth
# =============================================================================

# compute_ground_truth — true risk from full counterfactual simulation.
# Computes per-individual survival curves, then marginalizes over individuals.
compute_ground_truth <- function(sample.all, exposure.levels, interv.list) {

  # Filter to static intervention rows only; natural course (Intervention = NA) is dropped
  sample.ref <- sample.all %>%
    mutate(E = Intervention, Haz = Pr) %>%
    filter(!is.na(E))

  # Per-individual survival curves via haz_to_surv (group by E + Rep)
  sample.ref <- haz_to_surv(sample.ref, E, Rep)

  # Average across individuals per exposure intervention x interval
  surv <- sample.ref %>%
    group_by(E, Ts, Te) %>%
    summarize(
      Haz  = mean(Haz),
      Surv = mean(Surv),
      Risk = mean(Risk),
      .groups = "drop"
    ) %>%
    arrange(E, Ts)

  list(fit = NULL, surv = surv, risk = extract_risk(surv))
}

# =============================================================================
# PK imputation
# =============================================================================

# Reshape observed data into nlmixr2 dosing format.
# Unpivots C0/C1/C2 and Cmu0/Cmu1/Cmu2 into long rows (one per time-after-dose),
# then prepends a dosing event (EVID=1) at the start of each interval.
prep_pk_data <- function(data) {

  # Unpivot concentration columns: C0,C1,C2 → C at observation index K = 0,1,2 (SI's k); same for Cmu
  nl.data <- data %>%
    pivot_longer(cols = c(C0, C1, C2, Cmu0, Cmu1, Cmu2),
                 names_to      = c(".value", "K"),
                 names_pattern = "(C(?:mu)?)(\\d)") %>%
    mutate(K    = as.numeric(K),
           TaD  = K / 2.001,        # time after dose; slight offset avoids tied times
           TIME = Ts + TaD,
           CMT  = 1L) %>%
    select(ID, Te, TIME, TaD, CMT, D, C,
           C.true = Cmu, B,
           CL.true = CL, CL0.true = CL0, V0.true = V0, E.true = E)

  # Add dosing rows (EVID=1, AMT=dose) before each interval's observations
  nl.dosing <- nl.data %>%
    mutate(AMT = 0, EVID = 0, DV = C, DVID = CMT) %>%
    select(ID, Te, TIME, AMT, EVID, D, CMT, DVID, DV,
           C.true, B, CL.true, CL0.true, V0.true, E.true) %>%
    group_by(ID, Te) %>%
    group_modify(~add_row(.x,
      TIME     = first(.x$TIME),
      AMT      = first(.x$D),
      EVID     = 1,  D = 0,  CMT = 1L,  DVID = 1L,
      C.true   = first(.x$C.true),
      B        = first(.x$B),
      CL.true  = first(.x$CL.true),
      CL0.true = first(.x$CL0.true),
      V0.true  = first(.x$V0.true),
      E.true   = first(.x$E.true),
      DV       = NA,
      .before  = TRUE)) %>%
    ungroup() %>%
    select(-D, -Te)

  # Median baseline biomarker — fixed-effect anchor for PK model.
  # B[Ts==0] = B0 by the t=1 initial condition (B1 = B0); just B0's median.
  nl.dosing$B0.med.est <- data %>% filter(Ts == 0) %>% pull(B) %>% median()

  nl.dosing
}

# Fit one-compartment IV bolus PK model via nlmixr2 FOCEI.
# Returns nlmixr2 fit object, or NULL if use.true.pk=TRUE (skip fitting).
# ETAs are warm-started from true individual parameters to aid convergence.
fit_pk <- function(data, param, use.true.pk = FALSE) {
  if (use.true.pk) return(NULL)

  nl.dosing <- prep_pk_data(data)
  pk.model   <- make_pk_model(param)

  # Initial ETAs from true individual values for warm start.
  eta.init <- nl.dosing %>%
    distinct(ID, CL0.true, V0.true) %>%
    arrange(ID) %>%
    transmute(
      eta.V0  = log(V0.true)  - unname(param["v0.lmu"]),
      eta.CL0 = log(CL0.true) - unname(param["cl0.lmu"])
    ) %>%
    as.matrix()

  # Inner threads sized to avoid oversubscribing cores under outer (furrr) parallelism.
  n.rx.threads <- max(1, (available_cores() - 1) %/% future::nbrOfWorkers())
  rxode2::setRxThreads(n.rx.threads)
  fit <- nlmixr2(pk.model, data = nl.dosing, est = "focei",
                 control = foceiControl(
                   covMethod = 0,             # skip covariance — only need individual ETAs
                   etaMat = eta.init,         # warm-start from true individual params
                   maxInnerIterations = 1000  # generous inner limit for small-N convergence
                 ))

  gc()
  fit
}

# Impute per-interval/subject exposure (E.hat) from PK fit (or true E if use.true.pk).
# Returns person-period data.frame with E.hat, summary metrics (TAE, Early, SS), common input to all ER estimators
# Also add time-interval dummy variables for pooled logistic regression (survival).
impute_exposure <- function(pk.fit, data, Nints, ss.int) {
  
  # Impute E.hat: true exposure (oracle) or PK-derived from EBEs
  if (is.null(pk.fit)) {
    # use.true.pk: E.hat = E (true). Used in SI scenarios to isolate
    # estimator behavior from PK estimation error and for debugging.
    data.imputed <- data %>% mutate(E.hat = E)
  } else {
    # Predict individual exposure from PK fit
    pk.impute <- as.data.frame(pk.fit) %>%
      filter((TIME %% 1) == 0) %>%
      transmute(
        ID    = as.integer(ID),
        Ts    = as.integer(TIME),
        Te    = Ts + 1L,
        # E.hat = (IPRED/Δt) * (V/CL) * (1 - exp(-CL/V*Δt)); Δt = 1 so dropped
        E.hat = IPRED * V / CL * (1 - exp(-CL / V))
      )
    data.imputed <- data %>%
      left_join(pk.impute, by = c("ID", "Ts", "Te"))
  }

  # Per-subject summary metrics: time-averaged, early (interval 1), steady-state (interval ss.int)
  subj.exposures <- data.imputed %>%
    group_by(ID) %>%
    summarize(
      TAE   = mean(E.hat),   # obs censored at event, so avg is over observed intervals only
      Early = first(E.hat),
      SS    = E.hat[Te == ss.int][1],  # NA if subject has event before steady state
      .groups = "drop"
    )

  # Append summary metrics and time-interval dummies for pooled logistic regression
  ints      <- seq_len(Nints)
  int.names <- paste0("T", ints)

  # One-hot interval dummies (T1, T2, ...) — allow pooled logistic to estimate
  # a separate baseline hazard per interval (discrete-time analogue of baseline hazard)
  int.dummies <- data.frame(cbind(seq(0, Nints - 1), diag(Nints)))
  colnames(int.dummies) <- c("Ts", int.names)

  data.prep <- data.imputed %>%
    left_join(subj.exposures, by = "ID") %>%
    left_join(int.dummies, by = "Ts") %>%
    # SS=NA for intervals before steady state; don't condition on future exposure
    mutate(SS = if_else(Te < ss.int, NA_real_, SS))

  data.prep
}

# =============================================================================
# Estimators
# =============================================================================

# Standard logistic ER on one-row-per-subject data summarizing all of f/u.
# Predicts response probability at median B0; a conditional ER estimate.
fit_logistic <- function(data, exp.metric, adjust = NULL, exposure.levels, Nints,
                   ci = FALSE, n.boot = 100) {

  # Collapse person-period → one row per subject (discards temporal info)
  data.static <- data %>%
    group_by(ID) %>%
    summarize(
      R     = max(R),  #response event occurs (max(R)=1) or not (max(R)=0) during f/u
      B0    = first(B0),
      Early = first(Early),
      TAE   = first(TAE),
      SS    = SS[!is.na(SS)][1], # first non-NA; NA if event before ss.int
      .groups = "drop"
    )

  # Build formula: R ~ exposure + baseline biomarker
  covs <- c(exp.metric, adjust)
  fml  <- as.formula(paste("R ~", paste(covs, collapse = " + ")))

  # Prediction grid at target exposures with biomarker at sample median
  make_pred_grid <- function(d) {
    pg <- data.frame(x = exposure.levels)
    names(pg) <- exp.metric
    for (v in adjust) pg[[v]] <- median(d[[v]], na.rm = TRUE) # conditional on median B0
    pg
  }

  # Fit model and predict risk (hazard) across exposure grid
  pred.grid <- make_pred_grid(data.static)
  point.est <- glm_fit_predict(fml, data.static, pred.grid)
  glm.fit   <- point.est$fit

  # Logistic reg. 'hazard' is predicted risk P(R=1); reformat for common output format
  haz.df <- tibble(E = exposure.levels, Ts = 0, Te = Nints, Haz = point.est$haz)
  surv   <- haz_to_surv(haz.df, E)

  if (ci) {
    # Bootstrap CI by refitting GLM on resampled subjects, re-derive median B0 per replicate
    predict_fn <- function(d) {
      tibble(E = exposure.levels, Ts = 0, Te = Nints,
             Haz = glm_fit_predict(fml, d, make_pred_grid(d))$haz)
    }
    ci.bands <- bootstrap_surv_ci(data.static, n.boot, predict_fn)
    surv <- surv %>% left_join(ci.bands, by = c("E", "Ts", "Te"))
  }

  list(fit = glm.fit, surv = surv, risk = extract_risk(surv))
}

# Discrete-time survival via pooled logistic regression on person-period data.
# Static metrics use a per-subject summary (Early/SS/TAE); TV variants use per-interval exposure
# Interval-specific intercepts give saturated baseline hazard (≈ Cox PH).
fit_pooled_logistic <- function(data, exp.metric, adjust = NULL, exposure.levels, Nints, ss.int,
                    ci = FALSE, n.boot = 100) {

  Nexp      <- length(exposure.levels)
  ints      <- seq_len(Nints)
  ts        <- seq(0, Nints - 1)
  te        <- seq(1, Nints)
  int.names <- paste0("T", ints)

  # SS model drops pre-steady-state dummies since exposure is undefined there
  if (exp.metric == "SS") {
    model.dummies <- int.names[ints >= ss.int]
  } else {
    model.dummies <- int.names
  }

  # R ~ 0 + exposure + baseline biomarker + interval dummies (one per interval)
  covs <- c(exp.metric, adjust, model.dummies)
  fml  <- as.formula(paste("R ~ 0 +", paste(covs, collapse = " + ")))

  # Dummy matrix for interval intercepts
  int.dummies <- data.frame(cbind(ts, diag(Nints)))
  colnames(int.dummies) <- c("Ts", int.names)

  # Full prediction grid for every exposure level × interval combination
  make_pred_grid <- function(d) {
    pg <- tibble(
      Ts = rep(ts, Nexp),
      Te = rep(te, Nexp),
      E  = rep(exposure.levels, each = Nints)
    )
    pg[[exp.metric]] <- pg$E 
    for (v in adjust) {
      pg[[v]] <- median(d[[v]][d$Ts == 0], na.rm = TRUE) # conditional on median B0
    }
    pg %>% left_join(int.dummies, by = "Ts")
  }

  # Zero out hazard before steady state for SS models
  ss.mask <- if (exp.metric == "SS") rep(te < ss.int, Nexp) else NULL

  # Fit model and predict hazards across exposure grid
  pred.grid <- make_pred_grid(data)
  point.est <- glm_fit_predict(fml, data, pred.grid, mask = ss.mask)
  glm.fit   <- point.est$fit

  haz.df <- tibble(
    E   = pred.grid$E,
    Ts  = pred.grid$Ts,
    Te  = pred.grid$Te,
    Haz = point.est$haz
  )

  # Compute cumulative survival and risk from per-interval hazards 
  surv <- haz_to_surv(haz.df, E)

  if (ci) {
    # Bootstrap CI: refit GLM on resampled subjects, re-derive median B0 per replicate
    predict_fn <- function(d) {
      pg <- make_pred_grid(d)
      tibble(E = pg$E, Ts = pg$Ts, Te = pg$Te,
             Haz = glm_fit_predict(fml, d, pg, mask = ss.mask)$haz)
    }
    ci.bands <- bootstrap_surv_ci(data, n.boot, predict_fn)
    surv <- surv %>% left_join(ci.bands, by = c("E", "Ts", "Te"))
  }

  list(fit = glm.fit, surv = surv, risk = extract_risk(surv))
}

# Parametric g-formula for causal ER estimation under time-varying confounding.
# Uses gfoRmula package to fit biomarker model (gnls) and response model (PLR),
# runs Monte Carlo simulation under static exposure interventions to estimate ACERC.
fit_gformula <- function(data, param, exposure.levels, Nints, Nids,
                         ci = FALSE, n.boot = 100, parallel = TRUE,
                         fix.IE50 = FALSE, seed = 1234) {

  Nexp <- length(exposure.levels)
  gform_boot     <- ci && n.boot > 0
  gform_nsamples <- if (gform_boot) n.boot else 0

  # Subset to columns used by gfoRmula (drop any with structural NAs to avoid na.fail errors).
  # NOTE: overwrite E with the imputed E.hat so the g-formula uses the PK-imputed exposure via E
  gf.cols <- c("ID", "Ts", "Te", "R", "B0", "B", "E", "E.hat", "CL", "D", "Pi")
  data <- data %>% select(any_of(gf.cols)) %>% mutate(E = E.hat) %>% select(-E.hat)

  # Starting values use true params for warm start
  lkB.start   <- unname(param['kB.lmu'])
  lIE50.start <- unname(param['IE50.lmu'])

  # Create custom gnls nonlinear biomarker model fitting func. for use in gfoRmula
  bm_fit <- make_bm_fit(lkB.start, lIE50.start, fix.IE50)

  # Create biomarker prediction function matching the proportional error structure
  # floors at 0.001 to avoid degenerate values
  bm_pred <- function(obs_data, newdf, fitcov, time_name, t, condition, covname) {
    predicted <- as.numeric(predict(object = fitcov, newdata = newdf))
    pmax(rnorm(length(unique(obs_data$newid)),
                      predicted,
                      fitcov$sigma * abs(predicted)),
         0.001)
  }

  # Define static interventions that fix exposure to target values across all intervals.
  interventions <- map(exposure.levels, ~list(c(static, rep(.x, Nints))))
  int_descript  <- paste(exposure.levels, "Exp")

  gform_fit <- gformula(
    obs_data          = as.data.table(data),
    outcome_type      = 'survival',        # PLR for per-interval hazard (P_R)
    id                = 'ID',
    time_points       = Nints,
    time_name         = 'Ts',
    covnames          = c('E', 'B'),
    outcome_name      = 'R',
    covtypes          = c('normal', 'custom'),
    # gfoRmula requires one formula per covariate in covmodels, but neither is used:
    # E is overwritten by the static intervention; B uses covfits_custom (see make_bm_fit)
    covparams         = list(covmodels = c(
                          E ~ 1,
                          B ~ lag1_E + lag1_B + B0)), # dummy; must reference lag vars so gfoRmula creates them
    covfits_custom    = c(NA, bm_fit),
    covpredict_custom = c(NA, bm_pred),
    ymodel            = R ~ E + B,         # P_R: pooled logistic regression (per-interval hazard) model
    intvars           = as.list(rep('E', Nexp)),
    interventions     = interventions,
    int_descript      = int_descript,
    histories         = c(lagged),         # creates lag1_E, lag1_B for biomarker model
    histvars          = list(c('E', 'B')),
    basecovs          = 'B0',              # fixed per-subject throughout MC simulation
    nsimul            = Nids,              # MC trajectories for risk estimation (Step 3, SI C)
    seed              = seed,
    model_fits        = TRUE,              # store fitted P_B, P_R for diagnostics
    boot_diag         = gform_boot,        # bootstrap diagnostics when CI requested
    nsamples          = gform_nsamples,     # bootstrap resamples (0 = no CI)
    parallel          = parallel,
    ncores            = if (parallel) max(1, (available_cores() - 1) %/% future::nbrOfWorkers()) else NA,
    ref_int           = 1                  # reference intervention for risk ratio (unused)
  )

  # Extract results — uniform return structure
  gf.res <- gform_fit$result %>%
    filter(Interv. != 0) %>%                  #drop natural course (not meaningfull with E~1)
    mutate(E = exposure.levels[Interv.]) %>%
    rename(Ts = k, Risk = `g-form risk`) %>%
    rename(any_of(c(Risk.low = "Risk lower 95% CI", Risk.up = "Risk upper 95% CI")))

  # gformula returns cumulative risk directly; derive Surv and per-interval Haz
  surv <- gf.res %>%
    mutate(Te   = Ts + 1,
           Surv = 1 - Risk) %>%
    group_by(E) %>%
    arrange(E, Ts) %>%
    # Haz(t) = (Risk(t)-Risk(t-1))/Surv(t-1), Haz(1)=Risk(1) since Surv(0)=1
    mutate(Haz = c(Risk[1], diff(Risk) / lag(Surv)[-1])) %>%
    ungroup() %>%
    select(E, Ts, Te, Haz, Surv, Risk, any_of(c("Risk.low", "Risk.up")))

  # Derive survival CIs from risk CIs (bounds swap because Surv = 1 - Risk)
  if ("Risk.low" %in% names(surv)) {
    surv <- surv %>%
      mutate(Surv.low = 1 - Risk.up,
             Surv.up  = 1 - Risk.low)
  }

  list(fit = gform_fit, surv = surv, risk = extract_risk(surv))
}

# =============================================================================
# Results
# =============================================================================


# 6-panel overview of DGP trajectories (C, E, N, B, CL, Pi) by dose group.
# Plots the natural course only, censored at first event.
# Mean trajectory with 2.5–97.5% bands; patients-at-risk shown as step function.
plot_dgp_summary <- function(reference, param) {

  plot.vars <- c("C", "E", "N", "B", "CL", "Pi")
  plot.labs <- c(C  = "Concentration (a.u.)",  E  = "Exposure (a.u.)",
                 N  = "Fraction at Risk",       B  = "Biomarker (a.u.)",
                 CL = "Clearance (a.u.)",       Pi = "Dose Interruption Rate")

  # Natural course, censored at first event
  data.nc <- censor_natural_course(reference, Rep)

  traj.df <- bind_rows(
    # Continuous variables: C0/C1/C2 give 3 sub-interval concentration samples,
    # pivoted to fractional time points (/ 2.001 avoids exact overlap at boundaries).
    # Then all plotted variables are pivoted longer for faceting.
    data.nc %>%
      mutate(Group = paste0(A, " a.u.")) %>%
      pivot_longer(cols = c(C0, C1, C2),                       # 3 sub-interval conc samples → 1 "C" column
                   names_to = "K", names_pattern = "C(\\d)", values_to = "C") %>%
      mutate(TaD  = as.integer(K) / 2.001,                     # time after dose (SI's k → TaD)
             Time = Ts + TaD) %>%                               # fractional time within interval
      pivot_longer(cols = c(C, E, B, CL, Pi),                  # wide → long for faceting
                   names_to = "Variable", values_to = "Value") %>%
      group_by(Time, Group, Variable) %>%                      # mean + 95% band per time × group
      summarize(Mean = mean(Value), q025th = quantile(Value, .025),
                q975th = quantile(Value, .975), .groups = "drop"),
    # Patients at risk: normalized count (no uncertainty band)
    data.nc %>%
      mutate(Group = paste0(A, " a.u.")) %>%
      group_by(Group, Time = as.numeric(Ts)) %>%
      summarize(N = n_distinct(Rep), .groups = "drop") %>%     # subjects still alive
      group_by(Group) %>%
      mutate(Mean = N / max(N), q025th = Mean, q975th = Mean, Variable = "N") %>%  # normalize to proportion
      select(-N) %>% ungroup()
  ) %>%
    filter(Variable %in% plot.vars) %>%
    mutate(Variable = factor(recode(Variable, !!!plot.labs), levels = plot.labs[plot.vars]))

  # C and E share a y-axis as they are on the same scale; 
  ce.max <- traj.df %>%
    filter(Variable %in% c(plot.labs[['C']], plot.labs[['E']])) %>%
    pull(q975th) %>% max()
  # B and CL axis are set so init medians are aligned, allows easier tracking of trend
  b.max    <- traj.df %>% filter(Variable == plot.labs[['B']])  %>% pull(q975th) %>% max()
  cl.max   <- traj.df %>% filter(Variable == plot.labs[['CL']]) %>% pull(q975th) %>% max()
  bc.ratio <- exp(param['cl0.lmu']) / exp(param['b0.lmu'])
  b.upper  <- max(b.max,  cl.max  / bc.ratio)
  cl.upper <- max(cl.max, b.max   * bc.ratio)

  # Per-facet y-axis overrides (positional, must match plot.vars order):
  # C, E, N, B, CL, Pi
  traj.df %>%
    ggplot(aes(Time, Mean, color = Group, fill = Group)) +
    geom_ribbon(data = ~filter(., Variable != "Fraction at Risk"),
                aes(ymin = q025th, ymax = q975th), alpha = .2, color = NA) +
    geom_line(data = ~filter(., Variable != "Fraction at Risk")) +
    geom_step(data = ~filter(., Variable == "Fraction at Risk"), linewidth = 0.7) +
    scale_y_continuous(limits = \(x) c(0, x[2])) +
    facet_wrap2(~Variable, scales = "free_y", ncol = 3) +
    facetted_pos_scales(y = list(
      scale_y_continuous(limits = c(0, ce.max)),   # C
      scale_y_continuous(limits = c(0, ce.max)),   # E
      NULL,                                         # N (default)
      scale_y_continuous(limits = c(0, b.upper)),  # B
      scale_y_continuous(limits = c(0, cl.upper)), # CL
      scale_y_continuous(limits = c(0, 1))          # Pi
    )) +
    labs(x = "Time (interval)", y = NULL, color = "Dose (a.u.)", fill = "Dose (a.u.)")
}

# Summarize ER curves across simulation reps: mean/CI bands for plotting,
# point-wise bias/RMSE at two evaluation exposures, and ATE between them.
process_er_curves <- function(er.curves, eval.E) {

   # Mean risk and 95% interval (bands) across reps, per method
  bands <- er.curves %>%
    group_by(Method, E) %>%
    summarize(
      Risk.mean = mean(Risk),
      Risk.lo  = quantile(Risk, 0.025),
      Risk.hi  = quantile(Risk, 0.975),
      .groups  = "drop"
    )

  # Per-rep ground truth for RMSE (varies across reps due to different simulated datasets)
  truth <- er.curves %>%
    filter(Method == "Ground Truth") %>%
    select(sim, E, Risk.true = Risk)

  # Estimator risks at the two eval exposures, joined with per-rep ground truth for comparison
  eval.base <- er.curves %>%
    filter(abs(E - eval.E[1]) < 1e-9 | abs(E - eval.E[2]) < 1e-9) %>%
    left_join(truth, by = c("sim", "E")) %>%
    mutate(col = if_else(abs(E - eval.E[1]) < 1e-9, "R1", "R2"))

  # Point-wise mean risk, RMSE, and signed bias vs ground truth at each eval exposure
  # (ground truth included with RMSE = Bias = 0)
  eval.stats <- eval.base %>%
    group_by(Method, col) %>%
    summarize(
      Mean = mean(Risk, na.rm = TRUE),
      RMSE = sqrt(mean((Risk - Risk.true)^2, na.rm = TRUE)),
      Bias = mean(Risk - Risk.true, na.rm = TRUE),
      .groups = "drop"
    )

  # ATE: risk difference between high and low eval exposures (R2 - R1),
  # with RMSE/Bias comparing estimated vs true difference across reps
  ate.stats <- eval.base %>%
    select(sim, Method, col, Risk, Risk.true) %>%
    pivot_wider(names_from = col, values_from = c(Risk, Risk.true)) %>%
    mutate(Diff = Risk_R2 - Risk_R1,
           true.Diff = Risk.true_R2 - Risk.true_R1) %>%
    group_by(Method) %>%
    summarize(
      ATE      = mean(Diff, na.rm = TRUE),
      ATE.RMSE = sqrt(mean((Diff - true.Diff)^2, na.rm = TRUE)),
      ATE.Bias = mean(Diff - true.Diff, na.rm = TRUE),
      .groups  = "drop"
    )

  list(bands = bands, eval = eval.stats, ate = ate.stats)
}

# Formats sim study results into a comparison table (Table 1 in paper).
# One row per method with ATE and risk at two evaluation points, each with RMSE and signed Bias.
format_mse_table <- function(curve.stats, eval.E, plot.methods = NULL) {

  r1.lab  <- paste0("R(E=", eval.E[1], ")")
  r2.lab  <- paste0("R(E=", eval.E[2], ")")
  r1r.lab <- paste0("R(E=", eval.E[1], ") RMSE")
  r2r.lab <- paste0("R(E=", eval.E[2], ") RMSE")
  r1b.lab <- paste0("R(E=", eval.E[1], ") Bias")
  r2b.lab <- paste0("R(E=", eval.E[2], ") Bias")

  # Pivot eval points wide so each method is one row; join ATE.
  # Ground truth RMSE/Bias are 0 by construction — display as NA.
  eval.wide <- curve.stats$eval %>%
    pivot_wider(names_from = col, values_from = c(Mean, RMSE, Bias))

  all.rows <- eval.wide %>%
    left_join(curve.stats$ate, by = "Method") %>%
    transmute(
      Method,
      ATE        = ATE * 100,
      `ATE RMSE` = if_else(Method == "Ground Truth", NA_real_, ATE.RMSE * 100),
      `ATE Bias` = if_else(Method == "Ground Truth", NA_real_, ATE.Bias * 100),
      !!r1.lab  := Mean_R1 * 100,
      !!r1r.lab := if_else(Method == "Ground Truth", NA_real_, RMSE_R1 * 100),
      !!r1b.lab := if_else(Method == "Ground Truth", NA_real_, Bias_R1 * 100),
      !!r2.lab  := Mean_R2 * 100,
      !!r2r.lab := if_else(Method == "Ground Truth", NA_real_, RMSE_R2 * 100),
      !!r2b.lab := if_else(Method == "Ground Truth", NA_real_, Bias_R2 * 100)
    )

  # Ground truth first, then methods sorted by ATE RMSE (best first)
  gt.row      <- all.rows %>% filter(Method == "Ground Truth")
  method.rows <- all.rows %>% filter(Method != "Ground Truth") %>% arrange(`ATE RMSE`)

  # Keep table methods in sync with whatever the paper's figure plots; Ground
  # Truth is exempt (filters method.rows only) so it's always retained above.
  if (!is.null(plot.methods)) {
    method.rows <- method.rows %>% filter(Method %in% plot.methods)
  }

  bind_rows(gt.row, method.rows)  # Ground Truth row + (filtered) method rows
}

# Faceted ER plot — each panel shows one estimator's mean and 95% simulation
# interval overlaid on the ground truth (black dashed + grey ribbon).
plot_er_sim_study <- function(curve.stats, plot.methods, exposure.levels) {

  # Replicate ground truth bands into every estimator facet; convert risk to percent
  truth.overlay <- curve.stats$bands %>%
    filter(Method == "Ground Truth") %>%
    select(-Method) %>%
    crossing(tibble(Method = factor(plot.methods, levels = plot.methods))) %>%
    mutate(across(starts_with("Risk"), ~ .x * 100))

  # Subset and reorder estimator bands; convert risk to percent
  est.summary <- curve.stats$bands %>%
    filter(Method %in% plot.methods) %>%
    mutate(Method = factor(Method, levels = plot.methods)) %>%
    mutate(across(starts_with("Risk"), ~ .x * 100))

  # Grey ribbon + dashed line for ground truth; colored ribbon + solid line for estimator
  ggplot(est.summary, aes(x = E)) +
    geom_ribbon(data = truth.overlay,
                aes(ymin = Risk.lo, ymax = Risk.hi), fill = "black", alpha = 0.20) +
    geom_ribbon(aes(ymin = Risk.lo, ymax = Risk.hi, fill = Method), alpha = 0.25) +
    geom_line(aes(y = Risk.mean, color = Method)) +
    geom_line(data = truth.overlay,
              aes(y = Risk.mean), color = "black", linetype = "dashed") +
    facet_wrap(~Method, ncol = 3) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_fill_manual(values = METHOD_COLORS) +
    scale_x_continuous(breaks = exposure.levels) +
    coord_cartesian(ylim = c(0, 100)) +
    guides(color = "none", fill = "none") +
    labs(x = "Exposure (a.u.)",
         y = "Event Probability by End of Follow-Up (%)") +
    theme(panel.spacing = unit(0.8, "lines"))
}

# =============================================================================
# Supplemental
# =============================================================================

# Positivity diagnostic boxplots of exposure by biomarker tertile at each interval.
# Red dashed lines mark the exposure intervention range.
plot_positivity <- function(data, exposure.levels) {
  data %>%
    group_by(Ts) %>%
    mutate(B.tert = factor(cut(B,                              # within-interval biomarker tertiles
                               breaks = quantile(B, c(0, .33, .66, 1)),
                               labels = c('1st', '2nd', '3rd'),
                               include.lowest = TRUE))) %>%
    ggplot(aes(E, B.tert, group = B.tert, color = B.tert)) +
    geom_boxplot(alpha = .7, width = .3) +
    geom_vline(xintercept = c(min(exposure.levels), max(exposure.levels)),   # intervention range bounds
               color = 'red', linetype = 'dashed') +
    facet_grid(rows = vars(Ts), switch = "y") +                # one row per interval
    labs(y = "Interval", x = "E (a.u.)", color = "Biomarker\nTertile") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}


# PK model diagnostic plots. Residuals, parameter recovery, and individual fits.
# Outliers (IWRES > 6) are flagged and excluded from all panels.
plot_pk_diagnostics <- function(pk.fit, nl.dosing, data) {

  # Extract observed and true concentrations from dosing dataset (non-dosing rows only)
  pk.conc <- nl.dosing %>%
    filter(EVID == 0) %>%
    select(ID, TIME, C.obs = DV, C.true)

  # Flag extreme outliers that would distort diagnostic plots
  outlier.ids <- as.data.frame(pk.fit) %>%
    group_by(ID) %>%
    summarize(IWRES.max = max(IWRES), .groups = "drop") %>%
    filter(IWRES.max > 6) %>%
    pull(ID)
  message(paste("PK outliers (IWRES > 6):", length(outlier.ids)))

  # Merge fitted predictions with observed/true concentrations, excluding outliers
  pk.fit.df <- as.data.frame(pk.fit) %>%
    select(ID, TIME, IPRED, CL, V, IWRES) %>%
    mutate(ID = as.integer(ID)) %>%
    left_join(pk.conc, by = c("ID", "TIME")) %>%
    rename(C.pred = IPRED) %>%
    filter(!(ID %in% outlier.ids))

  # Residuals vs time — checks for time-dependent bias in the PK fit
  p.res.time <- ggplot(pk.fit.df, aes(TIME, IWRES)) +
    geom_point(fill = "steelblue", shape = 21, alpha = .4) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(x = "Time (interval)", y = "IWRES")

  # Residuals vs predicted concentration — checks for variance heterogeneity
  p.res.pred <- ggplot(pk.fit.df, aes(C.pred, IWRES)) +
    geom_point(fill = "steelblue", shape = 21, alpha = .4) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(x = "Fit Concentration (a.u.)", y = "IWRES")

  # Predicted vs true concentration — overall accuracy of individual predictions
  p.conc <- ggplot(pk.fit.df, aes(C.true, C.pred)) +
    geom_point(fill = "steelblue", shape = 21, alpha = .4) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(x = "True Concentration (a.u.)", y = "Fit Concentration (a.u.)")

  # Recover per-interval PK parameters (CL, V) and derived exposure from the fit,
  # then compare against true simulation values to assess parameter recovery.
  # E is back-calculated from fitted IPRED, V, CL using the analytical Cavg formula.
  pk.parm.long <- as.data.frame(pk.fit) %>%
    select(ID, TIME, IPRED, CL, V) %>%
    mutate(ID = as.integer(ID),
           Ts = as.integer(TIME),
           # E = Cavg = (IPRED/Δt) * (V/CL) * (1 - exp(-CL/V*Δt)); V=V0, Δt = 1 so dropped
           E  = IPRED * V / CL * (1 - exp(-CL / V))) %>%
    filter((TIME %% 1) == 0) %>%
    group_by(ID) %>% mutate(int.count = n()) %>% ungroup() %>%
    left_join(data %>% select(ID, Ts, E.true = E, CL.true = CL, V0.true = V0),
              by = c("ID", "Ts")) %>%
    filter(!(ID %in% outlier.ids)) %>%
    pivot_longer(cols = c(E, CL, V), names_to = "param", values_to = "val.fit") %>%
    mutate(val.true = case_when(param == "E"  ~ E.true,
                                param == "CL" ~ CL.true,
                                param == "V"  ~ V0.true),
           param    = factor(recode(param,
                                    E = "Exposure (a.u.)", CL = "CL (a.u.)", V = "V0 (a.u.)"),
                             levels = c("V0 (a.u.)", "CL (a.u.)", "Exposure (a.u.)"))) %>%
    # V0 is constant per subject (not time-varying like CL/E); collapse its
    # per-interval duplicate rows to one point per subject.
    distinct(ID, param, val.fit, val.true, int.count)

  # Fitted vs true parameter scatter, faceted by parameter, colored by on-study duration
  p.param <- ggplot(pk.parm.long, aes(val.true, val.fit, fill = int.count)) +
    geom_point(shape = 21, alpha = .5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    facet_wrap(~param, scales = "free") +
    labs(x = "True value", y = "Fit value", fill = "Intervals\nOn-study")

  # Individual concentration-time profiles for a random subset of 16 subjects.
  # Observed (points), true (grey), and fitted (red dashed) overlaid.
  set.seed(42)
  ind.ids <- sample(setdiff(unique(pk.fit.df$ID), outlier.ids), 16)

  p.ind <- pk.fit.df %>%
    filter(ID %in% ind.ids) %>%
    mutate(ID = factor(ID)) %>%
    ggplot(aes(x = TIME)) +
    geom_point(aes(y = C.obs),  fill = "steelblue", shape = 21, size = 1.2) +
    geom_line(aes(y = C.true),  color = "gray30") +
    geom_line(aes(y = C.pred),  color = "firebrick", linetype = "dashed") +
    facet_wrap(~ID, ncol = 4, scales = "free_y") +
    labs(x = "Time (interval)", y = "Concentration (a.u.)")

  # Return as named list for flexible arrangement in the report
  list(residuals = arrangeGrob(p.res.time, p.res.pred, p.conc, nrow = 1),
       params    = p.param,
       indfits   = p.ind)
}


# Per-scenario diagnostic plots comparing all estimators against ground truth.
# Returns survival (KM) step plot (p1), cumulative risk over time (p2), and
# end-of-follow-up ER curves faceted by estimator (p3; not used in SI).
plot_summary_outcomes <- function(results, exposure.levels, plot.methods = NULL) {

  method.order <- if (!is.null(plot.methods)) {
    intersect(plot.methods, names(results))
  } else {
    setdiff(names(results), "Ground Truth")
  }

  # Facet order for survival/cumrisk plots: GT first, then g-formula, then others
  surv.facet.order <- c("Ground Truth", "G-formula", "Adj. TV Survival", "TV Survival",
                         "Survival (Early)", "Survival (SS)", "Survival (TAE)")

  # LR only has final risk (one row per E), so exclude from time-series plots
  surv.methods <- names(results)[
    map_lgl(results, ~nrow(.x$surv) > length(exposure.levels))
  ]

  # Stack every method's surv table into one long data frame, tagged by Method;
  # factor levels preserve results' order for legends/facets.
  surv.all <- map_dfr(names(results), function(m) {
    results[[m]]$surv %>% mutate(Method = m)
  }) %>%
    mutate(Method = factor(Method, levels = names(results)))

  # Same stacking for the end-of-study risk tables.
  risk.all <- map_dfr(names(results), function(m) {
    results[[m]]$risk %>% mutate(Method = m)
  })

  # Extract ground truth to replicate into every estimator facet for overlay
  gt.risk <- risk.all %>% filter(Method == "Ground Truth")

  risk.faceted <- risk.all %>%
    filter(Method != "Ground Truth") %>%
    mutate(Facet = Method) %>%
    bind_rows(
      expand_grid(E = exposure.levels, Facet = unique(risk.all$Method[risk.all$Method != "Ground Truth"])) %>%
        left_join(gt.risk %>% select(-Method), by = "E") %>%
        mutate(Method = "Ground Truth")
    ) %>%
    mutate(Facet = factor(Facet, levels = method.order))

  trt.lab <- "E (a.u.)"

  # Shared prep for p1/p2: restrict to survival methods.
  # Adds baseline row w/ Te=0: Surv=1, Risk=0.
  surv.plot.df <- surv.all %>%
    filter(Method %in% surv.methods) %>%
    bind_rows(
      expand_grid(E = exposure.levels, Method = unique(.$Method)) %>%
        mutate(Ts = -1, Te = 0, Haz = 0, Surv = 1, Risk = 0)
    ) %>%
    mutate(E = factor(E),
           Method = factor(Method, levels = surv.facet.order),
           across(c(starts_with("Surv"), starts_with("Risk")), ~ .x * 100))

  # Survival step plot, survival methods only
  p1 <- surv.plot.df %>%
    ggplot() +
    geom_step(aes(Te, Surv, group = interaction(E, Method), color = E),
              linewidth = .8, alpha = 0.8) +
    # rect (not ribbon) to match the stepped Surv curve: one constant-height block
    # per interval, capped at xmax so the last interval doesn't overshoot the plot.
    {if ("Surv.low" %in% names(surv.plot.df))
      geom_rect(aes(xmin = Te, xmax = pmin(Te + 1, max(Te)), ymin = Surv.low, ymax = Surv.up,
                     group = interaction(E, Method), fill = E),
                alpha = .1, na.rm = TRUE)
    } +
    facet_grid(rows = vars(Method)) +
    labs(x = "Time (interval)", y = "Survival Probability (%)",
         color = trt.lab, fill = trt.lab)

  # Cumulative risk curves, survival methods only
  p2 <- surv.plot.df %>%
    ggplot() +
    geom_line(aes(Te, Risk, group = interaction(E, Method), color = E),
              alpha = 0.8) +
    # Only draw a band if CIs were actually computed (Risk.low present)
    {if ("Risk.low" %in% names(surv.plot.df))
      geom_ribbon(aes(Te, ymin = Risk.low, ymax = Risk.up,
                       group = interaction(E, Method), fill = E),
                  alpha = .1, na.rm = TRUE)
    } +
    facet_grid(rows = vars(Method)) +
    labs(x = "Time (interval)", y = "Cumulative Risk (%)",
         color = trt.lab, fill = trt.lab)

  # End-of-FU risk vs exposure, faceted by estimator with ground truth overlay
  p3 <- risk.faceted %>%
    mutate(Label = ifelse(Method == 'Ground Truth', 'Ground Truth', 'Causal Estimate'),
           across(starts_with("Risk"), ~ .x * 100)) %>%
    ggplot() +
    geom_line(aes(E, Risk, group = Method, color = Method, linetype = Label)) +
    # Only draw a band if CIs were actually computed (Risk.low present)
    {if ("Risk.low" %in% names(risk.faceted))
      geom_ribbon(aes(E, ymin = Risk.low, ymax = Risk.up, fill = Method),
                  alpha = .3, na.rm = TRUE)
    } +
    facet_wrap(~Facet) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_fill_manual(values = METHOD_COLORS) +
    coord_cartesian(ylim = c(0, 100)) +
    guides(color = "none", fill = "none", linetype = "none") +
    labs(x = "Exposure (a.u.)", y = "Event Probability by End of Follow-Up (%)")

  invisible(list(surv_step = p1, risk_curve = p2, risk_er = p3))
}

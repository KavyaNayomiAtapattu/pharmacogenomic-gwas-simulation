# =============================================================================
# Addressing Statistical Complexities in Modeling Treatment Effects from
# Repeated Measurements: A Simulation Study in Pharmacogenomics
#
# Author      : Kavya Atapattu Mudiyanselage
# Supervisor  : Prof. Matti Pirinen
# Institution : University of Helsinki, Faculty of Science
#               Master's Programme in Mathematics and Statistics
# Date        : May 2026
# =============================================================================
#
# OVERVIEW
#   This script implements the full Monte Carlo simulation study supporting
#   the three research objectives of the thesis. The motivating application
#   is statin-induced LDL-C change in a pharmacogenomic GWAS setting, where
#   treatment is assigned deterministically to individuals whose baseline
#   LDL-C exceeds a clinical threshold (top 10%). The simulation evaluates
#   how this threshold-based selection distorts genetic association estimates
#   through regression to the mean (RTM) and collider bias.
#
# RESEARCH OBJECTIVES
#   Objective 1 : Evaluate the unadjusted change score model M1 (DeltaY ~ G)
#                 under threshold and randomised treatment assignment.
#   Objective 2 : Compare M1 against the baseline-adjusted model M2
#                 (DeltaY ~ G + Y0), assessing collider bias from baseline
#                 conditioning.
#   Objective 3 : Evaluate the Dudbridge adjustment (applied to M2) and the
#                 Slope-Hunter correction via hunt() (applied to M1) as
#                 candidate bias-correction methods.
#
# SIMULATION DESIGN
#   N = 20,000 individuals per iteration | 1,000 Monte Carlo replications
#   MAF = 0.30 | Baseline heritability h2_Y0 = 0.50
#   220 SNPs across four classes: G0 (n=100), GD (n=10), G0D (n=10),
#   Null (n=100)
#   Variance decomposition: 50% genetic, 30% confounder, 20% residual noise
#
# SNP CLASSES
#   G0   : constant effect variants — affect Y0 and Y1 equally; no drug
#           interaction; true effect on DeltaY is zero
#   GD   : pharmacogenetic variants — drug interaction only; no main effect
#           on Y0; primary target of pharmacogenomic discovery
#   G0D  : pleiotropic variants — both a main effect on Y0 and a drug
#           interaction effect on Y1
#   Null : no effect on Y0 or Y1; type I error reference class
#
# TREATMENT MECHANISMS (evaluated in parallel within each iteration)
#   Threshold  : D = 1 if Y0 > Q_0.90(Y0), reflecting clinical statin
#                prescribing practice; induces RTM in the treated stratum
#   Randomised : D ~ Bernoulli(0.10), independent of Y0; serves as the
#                unbiased reference condition
#
# NAMING CONVENTIONS
#   Y0      pre-treatment baseline LDL-C
#   Y1      post-treatment follow-up LDL-C
#   DeltaY  change score: Y1 - Y0
#   M1      unadjusted change score model: DeltaY ~ G
#   M2      baseline-adjusted model: DeltaY ~ G + Y0
#   Dub     Dudbridge adjustment applied to M2 (treated stratum only)
#   SH      Slope-Hunter correction via hunt() applied to M1 (treated only)
#
# OUTPUTS  (written to OUTPUT_DIR)
#   Fig01_Histogram_Baseline_Y0.png
#   Fig02_Y0_by_TreatmentStatus.png
#   Fig03_VarDecomp_Y0.png
#   Fig04_VarDecomp_Y1.png
#   Fig05_VarDecomp_DeltaY.png
#   Fig06_Scatter_Y0vsY1.png
#   Fig07_Scatter_Y0vsDY_Combined.png
#   Fig08_Scatter_Y0vsDY_Randomised.png
#   Fig09_GenotypeMean_G0.png
#   Fig10_GenotypeMean_GD.png
#   Fig11_GenotypeMean_G0D.png
#   Fig12_GenotypeMean_GNull.png
#   Fig13_M1_Boxplot.png
#   Fig14_M1_QQ_Large.png
#   Fig15_G0_QQ_M1vsM2.png
#   Fig16_G0_QQ_M1vsM2vsDubvsSH.png
#   Fig17_M2_QQ.png
#   Fig18_Dub_QQ.png
#   Fig19_SH_QQ.png
#   Fig20_M2_Boxplot.png
#   Fig21_Dub_Boxplot.png
#   Fig22_SH_Boxplot.png
#   Table_Bias_M1_M2.csv
#   Table_Bias_Obj3.csv
#   simulation_workspace.RData   (if SAVE_WORKSPACE = TRUE)
# =============================================================================


# =============================================================================
# SECTION 0: PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
  library(grid)
  library(SlopeHunter)   #hunt() — official Slope-Hunter correction
})


# =============================================================================
# SECTION 1: GLOBAL SETTINGS
# =============================================================================

#1a. Workflow flags

SAVE_WORKSPACE <- TRUE

OUTPUT_DIR <- "simulation_outputs"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#1b.Core simulation parameters

N_ITERATIONS   <- 1000
N_patients     <- 20000
MAF            <- 0.30
target_h2_Y0   <- 0.50
fixed_beta_D   <- -3.00       #average drug effect (negative = LDL reduction)
TREAT_QUANTILE <- 0.90        #top 10% treated under Threshold scenario

n_G0 <- 100; n_GD <- 10; n_G0D <- 10; n_Null <- 100
N_snps <- n_G0 + n_GD + n_G0D + n_Null   # 220

#Variance calibration: 50% genetic, 30% confounder, 20% residual noise
var_g_single    <- 2 * MAF * (1 - MAF)
beta_0          <- sqrt(target_h2_Y0 / ((n_G0 + n_G0D) * var_g_single))
gamma_fixed     <- beta_0
beta_U_fixed    <- sqrt(0.30)
sigma_eps_fixed <- sqrt(0.20)

#SNP index sets
idx_G0   <- seq_len(n_G0)
idx_GD   <- (n_G0 + 1):(n_G0 + n_GD)
idx_G0D  <- (n_G0 + n_GD + 1):(n_G0 + n_GD + n_G0D)
idx_Null <- (n_G0 + n_GD + n_G0D + 1):N_snps

BETA_VEC   <- rep(0, N_snps); BETA_VEC[c(idx_G0, idx_G0D)]   <- beta_0
GAMMA_MASK <- rep(0, N_snps); GAMMA_MASK[c(idx_GD, idx_G0D)] <- 1

snp_meta <- data.frame(
  SNP      = seq_len(N_snps),
  SNP_Type = factor(
    c(rep("G0", n_G0), rep("GD", n_GD), rep("G0D", n_G0D), rep("Null", n_Null)),
    levels = c("G0", "GD", "G0D", "Null")
  ),
  gamma_true = GAMMA_MASK * gamma_fixed
)

#Dudbridge bootstrap replicates
DUB_N_BOOT <- 100

#Slope-Hunter: p-value threshold for incidence filter (passed to hunt())
#Bootstrapping enabled (M = 100) for proper SE(b) in Slope-Hunter correction.
SH_XPTHRESH <- 0.05

#Lambda-GC: exact expected median of chi-squared(1)
LAMBDA_GC_DIVISOR <- qchisq(0.5, df = 1)

#Stratum labels (used consistently throughout)
STRAT_TRT <- "Treated (D = 1)"
STRAT_UNT <- "Untreated (D = 0)"
STRAT_TRT_N <- "Treated (D = 1) [n\u22482,000]"
STRAT_UNT_N <- "Untreated (D = 0) [n\u224818,000]"

#1c. Colour palettes

MODEL_COLORS <- c(
  "M1"  = "#4477AA",   # steel blue
  "M2"  = "#EE6677",   # salmon red
  "Dub" = "#CCBB44",   # gold
  "SH"  = "#AA3377"    # magenta
)

MODEL_LABELS <- c(
  "M1"  = expression(paste("M1: ", Delta, "Y ~ G")),
  "M2"  = expression(paste("M2: ", Delta, "Y ~ G + Y"[0])),
  "Dub" = "Dudbridge (corrects M2)",
  "SH"  = "Slope-Hunter (corrects M1)"
)

SNP_COLORS <- c(
  "G0"   = "#E69F00",
  "GD"   = "#0072B2",
  "G0D"  = "#009E73",
  "Null" = "#999999"
)

STRATUM_COLORS <- c(
  "Treated (D = 1)"   = "#D6604D",
  "Untreated (D = 0)" = "#4A7FA5"
)

SCENARIO_COLS_2 <- c("Randomised" = "#4A7FA5", "Threshold" = "#D6604D")

VARDECOMP_Y0 <- c("Genetic" = "#4477AA", "Confounder" = "#EE6677",
                  "Noise"   = "#CCBB44")
VARDECOMP_Y1 <- c("Genetic" = "#4477AA", "Confounder" = "#EE6677",
                  "Drug (main effect)" = "#AA3377", "Noise" = "#CCBB44")
VARDECOMP_DY <- c("Genetic (interaction)" = "#4477AA",
                  "Confounder"            = "#EE6677",
                  "Drug (main effect)"    = "#AA3377",
                  "Noise"                 = "#CCBB44")

#1d. Font sizes 

FS_BASE       <- 13
FS_AXIS_TITLE <- 12
FS_AXIS_TEXT  <- 11
FS_STRIP      <- 12
FS_LEG        <- 10
FS_ANNOT_SM   <- 3.0

#1e. ggplot theme 

thesis_theme <- theme_classic(base_size = FS_BASE) +
  theme(
    strip.background   = element_rect(fill = "#F0F0F0", colour = "grey70"),
    strip.text         = element_text(face = "bold", size = FS_STRIP),
    axis.title         = element_text(size = FS_AXIS_TITLE),
    axis.text          = element_text(size = FS_AXIS_TEXT),
    legend.title       = element_text(face = "bold", size = FS_LEG + 1),
    legend.text        = element_text(size = FS_LEG),
    panel.grid.major.y = element_line(colour = "grey92"),
    plot.margin        = margin(6, 12, 6, 8)
  )

cat(sprintf(
  "Parameters: N=%d | MAF=%.2f | h2_Y0=%.2f | beta_0=%.5f | gamma=%.5f\n",
  N_patients, MAF, target_h2_Y0, beta_0, gamma_fixed
))


# =============================================================================
# SECTION 2: DATA-GENERATING FUNCTIONS
# =============================================================================

simulate_data <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G  <- matrix(rbinom(N_patients * N_snps, 2, MAF), nrow = N_patients)
  U  <- rnorm(N_patients)
  Y0 <- as.numeric(G %*% BETA_VEC) +
    beta_U_fixed * U +
    rnorm(N_patients, 0, sigma_eps_fixed)
  list(G = G, U = U, Y0 = Y0)
}

assign_treatment <- function(Y0, scenario) {
  if      (scenario == "Rand")      rbinom(N_patients, 1, 0.10)
  else if (scenario == "Threshold") as.integer(Y0 > quantile(Y0, TREAT_QUANTILE))
  else stop(paste("Unknown scenario:", scenario))
}

#generate_Y1: fresh independent noise eps1 is drawn each call.
#This independence between eps0 and eps1 is the mechanism for RTM under threshold.
generate_Y1 <- function(G, U, Y0, D) {
  as.numeric(G %*% BETA_VEC) +
    as.numeric(G %*% (GAMMA_MASK * gamma_fixed)) * D +
    beta_U_fixed * U +
    fixed_beta_D * D +
    rnorm(N_patients, 0, sigma_eps_fixed)
}

#compute_var_decomp: partition Var(Y0), Var(Y1), and Var(DeltaY) into
#interpretable components. Called per iteration to accumulate empirical means.
compute_var_decomp <- function(G, U, Y0, D, Y1) {
  
  decomp_y1 <- function(idx) {
    if (length(idx) < 2L) return(rep(NA_real_, 5L))
    g_st <- as.numeric(G[idx, ] %*% BETA_VEC)
    g_in <- as.numeric(G[idx, ] %*% (GAMMA_MASK * gamma_fixed)) * D[idx]
    cf   <- beta_U_fixed * U[idx]
    dr   <- fixed_beta_D * D[idx]
    eps  <- Y1[idx] - g_st - g_in - cf - dr
    tot  <- var(Y1[idx])
    c(var(g_st) / tot, var(g_in) / tot, var(cf) / tot, var(dr) / tot,
      var(eps) / tot)
  }
  
  decomp_dy <- function(idx) {
    if (length(idx) < 2L) return(rep(NA_real_, 3L))
    dY   <- Y1[idx] - Y0[idx]
    g_in <- as.numeric(G[idx, ] %*% (GAMMA_MASK * gamma_fixed)) * D[idx]
    dr   <- fixed_beta_D * D[idx]
    ns   <- dY - g_in - dr
    tot  <- var(dY)
    if (tot < 1e-10) return(rep(0, 3L))
    c(var(g_in) / tot, var(dr) / tot, var(ns) / tot)
  }
  
  g0   <- as.numeric(G %*% BETA_VEC)
  cf   <- beta_U_fixed * U
  eps0 <- Y0 - g0 - cf
  tot0 <- var(Y0)
  
  idx_t <- which(D == 1L); idx_u <- which(D == 0L)
  
  #Full-sample DeltaY decomposition (all N individuals, treated + untreated).
  #The confounder (beta_U * U) cancels exactly in Y1 - Y0, so conf = 0.
  #Components are normalised to sum to 1 (absorbing covariance terms).
  dY_all  <- Y1 - Y0
  g_in_fs <- as.numeric(G %*% (GAMMA_MASK * gamma_fixed)) * D
  dr_fs   <- fixed_beta_D * D
  ns_fs   <- dY_all - g_in_fs - dr_fs   # residual = eps1 - eps0
  raw_fs  <- c(var(g_in_fs), 0, var(dr_fs), var(ns_fs))
  DYfs    <- if (sum(raw_fs) > 1e-10) raw_fs / sum(raw_fs) else rep(0, 4L)
  
  list(
    Y0   = c(var(g0) / tot0, var(cf) / tot0, var(eps0) / tot0),
    Y1t  = decomp_y1(idx_t), Y1u = decomp_y1(idx_u),
    DYt  = decomp_dy(idx_t), DYu = decomp_dy(idx_u),
    DYfs = DYfs
  )
}


# =============================================================================
# SECTION 3: ANALYSIS FUNCTIONS
# =============================================================================

#3a. Dudbridge correction (applied to M2, treated stratum)
#
#Implements the Hedges-Olkin regression dilution correction for the
#regression of prognosis betas on incidence betas.  Bootstrap resampling
#estimates Var(b).  Applied to M2 estimates.

dudbridge_correct <- function(beta_prog, se_prog, beta_inc, se_inc,
                              n_boot = DUB_N_BOOT) {
  idx  <- seq_along(beta_prog)
  raw  <- sum(beta_prog * beta_inc) / sum(beta_inc^2)
  vb   <- var(beta_inc); ms2 <- mean(se_inc^2); denom <- vb - ms2
  b    <- if (abs(denom) > 1e-10) raw * vb / denom else raw
  
  b_bs <- replicate(n_boot, {
    r  <- sample(idx, replace = TRUE)
    rw <- sum(beta_prog[r] * beta_inc[r]) / sum(beta_inc[r]^2)
    vr <- var(beta_inc[r]); m2r <- mean(se_inc[r]^2); dr <- vr - m2r
    if (abs(dr) > 1e-10) rw * vr / dr else rw
  })
  var_b <- var(b_bs, na.rm = TRUE)
  
  ba   <- beta_prog - b * beta_inc
  se_a <- sqrt(se_prog^2 + b^2 * se_inc^2 +
                 beta_inc^2 * var_b + se_inc^2 * var_b)
  p_a  <- 2 * pnorm(-abs(ba / se_a))
  list(beta_adj = ba, se_adj = se_a, p_adj = p_a, b = b)
}

#3b. Slope-Hunter correction via SlopeHunter::hunt()
#
#Uses the official SlopeHunter package with Bootstrapping = TRUE so that
#proper SE(b) is available.  hunt() uses print()/cat() internally, so all
#console output is suppressed with capture.output() in addition to
#suppressMessages()/suppressWarnings().  The returned slope b and its
#bootstrap SE are applied manually:
# beta_adj = beta_M1 - b * beta_GY0
# se_adj   = sqrt(se_M1^2 + b^2*se_GY0^2 + beta_GY0^2*var_b + se_GY0^2*var_b)
#Applied to M1 estimates (treated stratum).

slopehunter_correct <- function(beta_inc, se_inc, beta_prog, se_prog, snp_ids) {
  p_inc  <- 2 * pnorm(-abs(beta_inc  / se_inc))
  p_prog <- 2 * pnorm(-abs(beta_prog / se_prog))
  
  sh_dat <- data.frame(
    SNP            = snp_ids,
    BETA.incidence = beta_inc,  SE.incidence = se_inc,   Pval.incidence = p_inc,
    BETA.prognosis = beta_prog, SE.prognosis = se_prog,  Pval.prognosis = p_prog
  )
  
  sh_out <- NULL
  invisible(capture.output(
    sh_out <- tryCatch(
      suppressMessages(suppressWarnings(
        hunt(sh_dat,
             snp_col          = "SNP",
             xbeta_col        = "BETA.incidence", xse_col = "SE.incidence",
             xp_col           = "Pval.incidence",
             ybeta_col        = "BETA.prognosis", yse_col = "SE.prognosis",
             yp_col           = "Pval.prognosis",
             xp_thresh        = SH_XPTHRESH,
             Bootstrapping    = TRUE,
             M                = 100,
             Plot             = FALSE,
             show_adjustments = FALSE)
      )),
      error = function(e) NULL
    )
  ))
  
  #Fall back to uncorrected M1 if hunt() fails or returns NA slope
  if (is.null(sh_out) || is.na(sh_out$b)) {
    return(list(beta_adj = beta_prog, se_adj = se_prog,
                p_adj    = p_prog,    b      = NA_real_,  bse = NA_real_))
  }
  
  b     <- sh_out$b
  bse   <- if (!is.null(sh_out$bse) && !is.na(sh_out$bse)) sh_out$bse else 0
  var_b <- bse^2
  
  ba   <- beta_prog - b * beta_inc
  se_a <- pmax(
    sqrt(se_prog^2 + b^2 * se_inc^2 + beta_inc^2 * var_b + se_inc^2 * var_b),
    1e-10
  )
  list(beta_adj = ba, se_adj = se_a,
       p_adj    = 2 * pnorm(-abs(ba / se_a)),
       b        = b,
       bse      = bse)
}


# =============================================================================
# SECTION 4: SEED-FIXED EXPLORATORY DATASET
# =============================================================================

cat("Generating seed-fixed dataset (seed = 2024) for exploratory figures ...\n")
dat_seed     <- simulate_data(seed = 2024)
D_thr_seed   <- assign_treatment(dat_seed$Y0, "Threshold")
D_rand_seed  <- assign_treatment(dat_seed$Y0, "Rand")
Y1_thr_seed  <- generate_Y1(dat_seed$G, dat_seed$U, dat_seed$Y0, D_thr_seed)
Y1_rand_seed <- generate_Y1(dat_seed$G, dat_seed$U, dat_seed$Y0, D_rand_seed)
tau_seed     <- quantile(dat_seed$Y0, TREAT_QUANTILE)
cat(sprintf("  n_trt(Threshold) = %d | n_trt(Randomised) = %d\n\n",
            sum(D_thr_seed), sum(D_rand_seed)))


# =============================================================================
# SECTION 5: MONTE CARLO LOOP
# =============================================================================

cat(sprintf("Starting Monte Carlo loop (%d iterations) ...\n", N_ITERATIONS))

res_list <- vector("list", N_ITERATIONS * 2L)
vd_list  <- vector("list", N_ITERATIONS * 2L)
k        <- 1L

for (iter in seq_len(N_ITERATIONS)) {
  if (iter %% 100L == 0L)
    cat(sprintf("  Iteration %d / %d\n", iter, N_ITERATIONS))
  
  dat_i <- simulate_data(seed = iter)
  
  for (sc in c("Rand", "Threshold")) {
    sc_lab <- if (sc == "Rand") "Randomised" else "Threshold"
    
    D_i  <- assign_treatment(dat_i$Y0, sc)
    Y1_i <- generate_Y1(dat_i$G, dat_i$U, dat_i$Y0, D_i)
    dY_i <- Y1_i - dat_i$Y0
    
    idx_t <- which(D_i == 1L); n_t <- length(idx_t)
    idx_u <- which(D_i == 0L); n_u <- length(idx_u)
    
    if (n_t < 20L || n_u < 20L) { k <- k + 1L; next }
    
    bM1  <- seM1  <- pM1  <-
      bM1u <- seM1u <- pM1u <-
      bM2  <- seM2  <- pM2  <-
      bM2u <- seM2u <- pM2u <-
      bGY0 <- seGY0 <- numeric(N_snps)
    
    for (s in seq_len(N_snps)) {
      gs <- dat_i$G[, s]
      
      m <- coef(summary(lm(dY_i[idx_t] ~ gs[idx_t])))
      bM1[s]  <- m[2,1]; seM1[s]  <- m[2,2]; pM1[s]  <- m[2,4]
      
      m <- coef(summary(lm(dY_i[idx_t] ~ gs[idx_t] + dat_i$Y0[idx_t])))
      bM2[s]  <- m[2,1]; seM2[s]  <- m[2,2]; pM2[s]  <- m[2,4]
      
      m <- coef(summary(lm(dY_i[idx_u] ~ gs[idx_u])))
      bM1u[s] <- m[2,1]; seM1u[s] <- m[2,2]; pM1u[s] <- m[2,4]
      
      m <- coef(summary(lm(dY_i[idx_u] ~ gs[idx_u] + dat_i$Y0[idx_u])))
      bM2u[s] <- m[2,1]; seM2u[s] <- m[2,2]; pM2u[s] <- m[2,4]
      
      m <- coef(summary(lm(dat_i$Y0 ~ gs)))
      bGY0[s] <- m[2,1]; seGY0[s] <- m[2,2]
    }
    
    #Dudbridge: corrects M2 in treated stratum
    dub <- tryCatch(
      dudbridge_correct(bM2, seM2, bGY0, seGY0),
      error = function(e) list(beta_adj = rep(NA, N_snps),
                               se_adj   = rep(NA, N_snps),
                               p_adj    = rep(NA, N_snps), b = NA_real_)
    )
    
    #Slope-Hunter: corrects M1 in treated stratum
    sh  <- slopehunter_correct(bGY0, seGY0, bM1, seM1, snp_meta$SNP)
    
    res_list[[k]] <- data.frame(
      iter = iter, scenario = sc_lab, snp = seq_len(N_snps),
      n_trt = n_t, n_unt = n_u,
      bM1  = bM1,       seM1  = seM1,     pM1  = pM1,
      bM1u = bM1u,      seM1u = seM1u,    pM1u = pM1u,
      bM2  = bM2,       seM2  = seM2,     pM2  = pM2,
      bM2u = bM2u,      seM2u = seM2u,    pM2u = pM2u,
      bDub = dub$beta_adj, seDub = dub$se_adj, pDub = dub$p_adj,
      bSH  = sh$beta_adj,  seSH  = sh$se_adj,  pSH  = sh$p_adj,
      b_dub = dub$b, b_sh = sh$b
    )
    
    vd <- compute_var_decomp(dat_i$G, dat_i$U, dat_i$Y0, D_i, Y1_i)
    vd_list[[k]] <- data.frame(
      iter = iter, scenario = sc_lab,
      Y0_gen   = vd$Y0[1],  Y0_conf  = vd$Y0[2],  Y0_noise = vd$Y0[3],
      Y1t_gen  = vd$Y1t[1], Y1t_gint = vd$Y1t[2], Y1t_conf = vd$Y1t[3],
      Y1t_drug = vd$Y1t[4], Y1t_noise= vd$Y1t[5],
      Y1u_gen  = vd$Y1u[1], Y1u_gint = vd$Y1u[2], Y1u_conf = vd$Y1u[3],
      Y1u_drug = vd$Y1u[4], Y1u_noise= vd$Y1u[5],
      DYt_gint = vd$DYt[1], DYt_drug = vd$DYt[2], DYt_noise= vd$DYt[3],
      DYu_gint = vd$DYu[1], DYu_drug = vd$DYu[2], DYu_noise= vd$DYu[3],
      DYfs_gint= vd$DYfs[1],DYfs_conf= vd$DYfs[2],
      DYfs_drug= vd$DYfs[3],DYfs_noise=vd$DYfs[4]
    )
    k <- k + 1L
  }
}

cat("Monte Carlo loop complete.\n\n")


# =============================================================================
# SECTION 6: POST-LOOP DATA PROCESSING
# =============================================================================

cat("Processing results ...\n")

full_res <- bind_rows(res_list) %>%
  left_join(snp_meta, by = c("snp" = "SNP")) %>%
  mutate(scenario = factor(scenario, levels = c("Randomised", "Threshold")))

vd_df <- bind_rows(vd_list) %>%
  mutate(scenario = factor(scenario, levels = c("Randomised", "Threshold")))

#Create tidy long-format subsets for each model × stratum combination.
#These are used directly by all QQ plot and boxplot figure functions.

mfl <- function(df, b_col, se_col, p_col, stratum_val, model_val) {
  df %>%
    select(iter, scenario, SNP_Type,
           beta = all_of(b_col), se = all_of(se_col), p = all_of(p_col)) %>%
    mutate(model = model_val, stratum = stratum_val)
}

fl_M1t  <- mfl(full_res, "bM1",  "seM1",  "pM1",  STRAT_TRT, "M1")
fl_M1u  <- mfl(full_res, "bM1u", "seM1u", "pM1u", STRAT_UNT, "M1")
fl_M2t  <- mfl(full_res, "bM2",  "seM2",  "pM2",  STRAT_TRT, "M2")
fl_M2u  <- mfl(full_res, "bM2u", "seM2u", "pM2u", STRAT_UNT, "M2")
fl_Dubt <- mfl(full_res, "bDub", "seDub", "pDub", STRAT_TRT, "Dub")
fl_SHt  <- mfl(full_res, "bSH",  "seSH",  "pSH",  STRAT_TRT, "SH")

fl_M1_all <- bind_rows(fl_M1t, fl_M1u) %>%
  mutate(stratum = factor(stratum, levels = c(STRAT_TRT, STRAT_UNT)))
fl_M2_all <- bind_rows(fl_M2t, fl_M2u) %>%
  mutate(stratum = factor(stratum, levels = c(STRAT_TRT, STRAT_UNT)))

cat("Done.\n\n")


# =============================================================================
# SECTION 7: EMPIRICAL VARIANCE DECOMPOSITION (averaged over iterations)
# =============================================================================

cat("Computing empirical variance decomposition ...\n")

#Ordered stratum levels used by Fig04 and Fig05
STRAT_LEVELS_3 <- c("Full Sample", STRAT_TRT, STRAT_UNT)

#Y0: averaged per scenario (panels will be identical; shown in two facets
#to make the scenario-independence of Y0 explicit)
vd_y0_bar <- vd_df %>%
  group_by(scenario) %>%
  summarise(Genetic    = mean(Y0_gen,   na.rm = TRUE) * 100,
            Confounder = mean(Y0_conf,  na.rm = TRUE) * 100,
            Noise      = mean(Y0_noise, na.rm = TRUE) * 100,
            .groups = "drop") %>%
  pivot_longer(-scenario, names_to = "Component", values_to = "Pct") %>%
  mutate(
    Component = factor(Component, levels = c("Genetic", "Confounder", "Noise")),
    Stratum   = "Y0"           # single x position per panel
  ) %>%
  group_by(scenario) %>%
  arrange(Component) %>%
  mutate(ypos  = cumsum(Pct) - Pct / 2,
         label = ifelse(Pct > 2, sprintf("%.1f%%", Pct), "")) %>%
  ungroup()

#Y1: three strata per scenario 
#Full Sample uses Y0 decomp as baseline reference (pre-treatment).
#Treated and Untreated combine the stable (G0) and interaction genetic terms
#into a single "Genetic" component, as both drive Y1 in their stratum.
vd_y1_bar <- vd_df %>%
  group_by(scenario) %>%
  summarise(
    fs_gen   = mean(Y0_gen,             na.rm = TRUE),
    fs_conf  = mean(Y0_conf,            na.rm = TRUE),
    fs_drug  = 0,                        # no drug at Y0 baseline
    fs_noise = mean(Y0_noise,           na.rm = TRUE),
    trt_gen  = mean(Y1t_gen + Y1t_gint, na.rm = TRUE),
    trt_conf = mean(Y1t_conf,           na.rm = TRUE),
    trt_drug = mean(Y1t_drug,           na.rm = TRUE),
    trt_noise= mean(Y1t_noise,          na.rm = TRUE),
    unt_gen  = mean(Y1u_gen + Y1u_gint, na.rm = TRUE),
    unt_conf = mean(Y1u_conf,           na.rm = TRUE),
    unt_drug = mean(Y1u_drug,           na.rm = TRUE),
    unt_noise= mean(Y1u_noise,          na.rm = TRUE),
    .groups = "drop") %>%
  pivot_longer(-scenario, names_to = "key", values_to = "Proportion") %>%
  mutate(
    Proportion = Proportion * 100,
    Stratum    = case_when(
      startsWith(key, "fs_")  ~ "Full Sample",
      startsWith(key, "trt_") ~ STRAT_TRT,
      startsWith(key, "unt_") ~ STRAT_UNT
    ),
    Component  = case_when(
      endsWith(key, "_gen")   ~ "Genetic",
      endsWith(key, "_conf")  ~ "Confounder",
      endsWith(key, "_drug")  ~ "Drug (main effect)",
      endsWith(key, "_noise") ~ "Noise"
    ),
    Stratum   = factor(Stratum,   levels = STRAT_LEVELS_3),
    Component = factor(Component, levels = c("Genetic", "Confounder",
                                             "Drug (main effect)", "Noise"))
  ) %>%
  group_by(scenario, Stratum) %>%
  arrange(Component) %>%
  mutate(ypos  = cumsum(Proportion) - Proportion / 2,
         label = ifelse(Proportion > 2, sprintf("%.1f%%", Proportion), "")) %>%
  ungroup()

#DeltaY: three strata per scenario
#Full Sample uses the full-sample DeltaY decomposition (DYfs).
#Confounder = 0 for Treated and Untreated strata because beta_U*U cancels
#exactly in Y1 - Y0.  It is retained in the legend for completeness and to
#make this cancellation visible.
vd_dy_bar <- vd_df %>%
  group_by(scenario) %>%
  summarise(
    fs_gint  = mean(DYfs_gint,  na.rm = TRUE),
    fs_conf  = mean(DYfs_conf,  na.rm = TRUE),   # = 0 by construction
    fs_drug  = mean(DYfs_drug,  na.rm = TRUE),
    fs_noise = mean(DYfs_noise, na.rm = TRUE),
    trt_gint = mean(DYt_gint,   na.rm = TRUE),
    trt_drug = mean(DYt_drug,   na.rm = TRUE),
    trt_noise= mean(DYt_noise,  na.rm = TRUE),
    unt_gint = mean(DYu_gint,   na.rm = TRUE),
    unt_drug = mean(DYu_drug,   na.rm = TRUE),
    unt_noise= mean(DYu_noise,  na.rm = TRUE),
    .groups = "drop") %>%
  mutate(trt_conf = 0, unt_conf = 0) %>%
  pivot_longer(-scenario, names_to = "key", values_to = "Proportion") %>%
  mutate(
    Proportion = Proportion * 100,
    Stratum    = case_when(
      startsWith(key, "fs_")  ~ "Full Sample",
      startsWith(key, "trt_") ~ STRAT_TRT,
      startsWith(key, "unt_") ~ STRAT_UNT
    ),
    Component  = case_when(
      endsWith(key, "_gint")  ~ "Genetic (interaction)",
      endsWith(key, "_conf")  ~ "Confounder",
      endsWith(key, "_drug")  ~ "Drug (main effect)",
      endsWith(key, "_noise") ~ "Noise"
    ),
    Stratum   = factor(Stratum,   levels = STRAT_LEVELS_3),
    Component = factor(Component, levels = c("Genetic (interaction)", "Confounder",
                                             "Drug (main effect)", "Noise"))
  ) %>%
  group_by(scenario, Stratum) %>%
  arrange(Component) %>%
  mutate(ypos  = cumsum(Proportion) - Proportion / 2,
         label = ifelse(Proportion > 2, sprintf("%.1f%%", Proportion), "")) %>%
  ungroup()

cat("Done.\n\n")


# =============================================================================
# SECTION 8: QQ AND BOXPLOT HELPER FUNCTIONS
# =============================================================================

#8a. 95% pointwise confidence envelope for a QQ plot

qq_envelope <- function(n_snps, alpha = 0.95) {
  k    <- seq_len(n_snps)
  p_lo <- qbeta((1 - alpha) / 2, k, n_snps - k + 1)
  p_hi <- qbeta((1 + alpha) / 2, k, n_snps - k + 1)
  data.frame(
    expected = rev(-log10((k - 0.5) / n_snps)),
    lo       = rev(-log10(p_hi)),
    hi       = rev(-log10(p_lo))
  )
}

#8b. Lambda-GC

compute_lambda <- function(p_vec) {
  chi2 <- qchisq(1 - p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1], df = 1)
  if (length(chi2) == 0L) return(NA_real_)
  median(chi2) / LAMBDA_GC_DIVISOR
}

#8c. Add stratum label with approximate n

add_stratum_lab <- function(df) {
  df %>% mutate(
    stratum_lab = dplyr::case_when(
      stratum == STRAT_TRT ~ STRAT_TRT_N,
      stratum == STRAT_UNT ~ STRAT_UNT_N,
      TRUE ~ as.character(stratum)
    ),
    stratum_lab = factor(stratum_lab, levels = c(STRAT_TRT_N, STRAT_UNT_N))
  )
}

#8d. Build QQ data frame with expected and lambda

build_qq_df <- function(fl) {
  fl %>%
    filter(!is.na(p), p > 0, p <= 1) %>%
    group_by(SNP_Type, stratum, scenario) %>%
    arrange(p) %>%
    mutate(n        = n(),
           expected = -log10((seq_along(p) - 0.5) / n),
           observed = -log10(p)) %>%
    ungroup()
}

build_qq_df_model <- function(fl) {
  fl %>%
    filter(!is.na(p), p > 0, p <= 1) %>%
    group_by(model, stratum, scenario) %>%
    arrange(p) %>%
    mutate(n        = n(),
           expected = -log10((seq_along(p) - 0.5) / n),
           observed = -log10(p)) %>%
    ungroup()
}

#8e. Build confidence envelope for a QQ (based on one class, one iter)

build_envelope <- function(qq_df, group_vars = c("stratum_lab", "scenario")) {
  nn <- qq_df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(n_e = round(first(n) / length(unique(qq_df$SNP_Type))),
              .groups = "drop")
  nn %>% rowwise() %>%
    mutate(env = list(qq_envelope(n_e))) %>%
    unnest(env) %>% ungroup()
}

build_envelope_model <- function(qq_df, group_vars = c("stratum_lab", "scenario")) {
  nn <- qq_df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(n_e = first(n), .groups = "drop")
  nn %>% rowwise() %>%
    mutate(env = list(qq_envelope(n_e))) %>%
    unnest(env) %>% ungroup()
}

#8f. QQ plot: all SNP classes, one model, all strata × scenarios
#Used for: Fig14 (M1 large), Fig17 (M2), Fig18 (Dub treated), Fig19 (SH treated)

make_qq_allclass <- function(fl_all, treated_only = FALSE,
                             annot_sz = FS_ANNOT_SM) {
  df <- fl_all %>%
    {if (treated_only) filter(., stratum == STRAT_TRT) else .} %>%
    add_stratum_lab() %>%
    mutate(SNP_Type = factor(SNP_Type, levels = c("G0","GD","G0D","Null")))
  
  qq_df  <- build_qq_df(df)
  env_df <- build_envelope(qq_df)
  
  max_obs <- max(qq_df$observed, na.rm = TRUE)
  max_exp <- max(qq_df$expected, na.rm = TRUE)
  
  line_tips <- qq_df %>%
    group_by(SNP_Type, stratum_lab, scenario) %>%
    slice_max(observed, n = 1, with_ties = FALSE) %>%
    summarise(x_tip = first(expected), y_tip = first(observed), .groups = "drop")
  
  lam_ann <- qq_df %>%
    group_by(SNP_Type, stratum_lab, scenario) %>%
    summarise(lam = compute_lambda(p), .groups = "drop") %>%
    group_by(stratum_lab, scenario) %>%
    arrange(desc(lam)) %>%
    mutate(rank_i = row_number(),
           y_pos  = max_obs * 0.97 - (rank_i - 1) * max_obs * 0.12,
           x_pos  = max_exp * 0.97,
           lbl    = sprintf("lambda[GC] == %.2f", lam)) %>%
    ungroup() %>%
    left_join(line_tips, by = c("SNP_Type", "stratum_lab", "scenario"))
  
  ggplot(qq_df, aes(expected, observed, colour = SNP_Type)) +
    geom_ribbon(data = env_df,
                aes(x = expected, ymin = lo, ymax = hi),
                inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_segment(data = lam_ann,
                 aes(x = x_pos, y = y_pos, xend = x_tip, yend = y_tip,
                     colour = SNP_Type),
                 linewidth = 0.25, alpha = 0.6,
                 inherit.aes = FALSE, show.legend = FALSE) +
    geom_label(data = lam_ann,
               aes(x = x_pos, y = y_pos, label = lbl,
                   colour = SNP_Type, fill = SNP_Type),
               parse = TRUE, size = annot_sz, hjust = 1,
               label.size = 0.25, label.padding = unit(0.15, "lines"),
               show.legend = FALSE) +
    scale_colour_manual(values = SNP_COLORS, name = "SNP class") +
    scale_fill_manual(values  = scales::alpha(SNP_COLORS, 0.15), guide = "none") +
    facet_grid(stratum_lab ~ scenario) +
    labs(x = expression("Expected" ~ -log[10](p)),
         y = expression("Observed" ~ -log[10](p))) +
    thesis_theme +
    theme(strip.text.y    = element_text(angle = -90, face = "bold"),
          legend.position = "bottom")
}

#8g. QQ comparison: G0 only, multiple models, all strata × scenarios
#Used for Fig15 (M1 vs M2) and Fig16 (all four)

make_qq_model_comparison <- function(fl_named_list, models,
                                     annot_sz = FS_ANNOT_SM) {
  combined <- bind_rows(lapply(models, function(m) {
    fl_named_list[[m]] %>%
      filter(SNP_Type == "G0") %>%
      mutate(model_f = factor(m, levels = models))
  })) %>%
    filter(!is.na(p), p > 0, p <= 1) %>%
    add_stratum_lab() %>%
    mutate(stratum_lab = droplevels(stratum_lab))
  
  qq_df <- combined %>%
    group_by(model_f, stratum_lab, scenario) %>%
    arrange(p) %>%
    mutate(n = n(), expected = -log10((seq_along(p) - 0.5) / n),
           observed = -log10(p)) %>% ungroup()
  
  env_df <- build_envelope_model(qq_df, c("stratum_lab", "scenario"))
  
  max_obs <- max(qq_df$observed, na.rm = TRUE)
  max_exp <- max(qq_df$expected, na.rm = TRUE)
  
  line_tips <- qq_df %>%
    group_by(model_f, stratum_lab, scenario) %>%
    slice_max(observed, n = 1, with_ties = FALSE) %>%
    summarise(x_tip = first(expected), y_tip = first(observed), .groups = "drop")
  
  lam_ann <- qq_df %>%
    group_by(model_f, stratum_lab, scenario) %>%
    summarise(lam = compute_lambda(p), .groups = "drop") %>%
    group_by(stratum_lab, scenario) %>%
    arrange(desc(lam)) %>%
    mutate(rank_i = row_number(),
           y_pos  = max_obs * 0.97 - (rank_i - 1) * max_obs * 0.12,
           x_pos  = max_exp * 0.97,
           lbl    = sprintf("lambda[GC] == %.2f", lam)) %>%
    ungroup() %>%
    left_join(line_tips, by = c("model_f", "stratum_lab", "scenario"))
  
  ggplot(qq_df, aes(expected, observed, colour = model_f)) +
    geom_ribbon(data = env_df,
                aes(x = expected, ymin = lo, ymax = hi),
                inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_segment(data = lam_ann,
                 aes(x = x_pos, y = y_pos, xend = x_tip, yend = y_tip,
                     colour = model_f),
                 linewidth = 0.25, alpha = 0.6,
                 inherit.aes = FALSE, show.legend = FALSE) +
    geom_label(data = lam_ann,
               aes(x = x_pos, y = y_pos, label = lbl,
                   colour = model_f, fill = model_f),
               parse = TRUE, size = annot_sz, hjust = 1,
               label.size = 0.25, label.padding = unit(0.15, "lines"),
               show.legend = FALSE) +
    scale_colour_manual(values = MODEL_COLORS[models],
                        labels = MODEL_LABELS[models], name = "Model") +
    scale_fill_manual(values  = scales::alpha(MODEL_COLORS[models], 0.15),
                      guide = "none") +
    facet_grid(stratum_lab ~ scenario) +
    labs(x = expression("Expected" ~ -log[10](p)),
         y = expression("Observed" ~ -log[10](p))) +
    thesis_theme +
    theme(strip.text.y    = element_text(angle = -90, face = "bold"),
          legend.position = "bottom")
}

#8h. Boxplot: distribution of beta-hat across iterations

make_boxplot <- function(fl_all, treated_only = FALSE) {
  df <- fl_all %>%
    {if (treated_only) filter(., stratum == STRAT_TRT) else .} %>%
    add_stratum_lab() %>%
    filter(!is.na(beta)) %>%
    mutate(SNP_Type = factor(SNP_Type, levels = c("G0","GD","G0D","Null")))
  
  p <- ggplot(df, aes(x = SNP_Type, y = beta,
                      fill = SNP_Type, colour = SNP_Type)) +
    geom_boxplot(alpha = 0.70, outlier.size = 0.40, linewidth = 0.50) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.55) +
    scale_fill_manual(values   = SNP_COLORS, name = "SNP class") +
    scale_colour_manual(values = SNP_COLORS, guide = "none") +
    labs(x = "SNP class", y = expression(hat(beta))) +
    thesis_theme +
    theme(legend.position = "bottom")
  
  if (treated_only) {
    p + facet_wrap(~ scenario, nrow = 1) +
      labs(subtitle = STRAT_TRT_N)
  } else {
    p + facet_grid(stratum_lab ~ scenario) +
      theme(strip.text.y = element_text(angle = -90, face = "bold"))
  }
}


# =============================================================================
# SECTION 9: EXPLORATORY FIGURES (Fig01 – Fig08)
# =============================================================================

cat("Generating exploratory figures (Fig01 - Fig08) ...\n")

#Fig01: Histogram of baseline Y0

fig01 <- ggplot(data.frame(Y0 = dat_seed$Y0), aes(x = Y0)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "#4477AA", colour = "white", alpha = 0.75) +
  geom_density(colour = "#D6604D", linewidth = 0.8) +
  geom_vline(xintercept = tau_seed, linetype = "dashed",
             colour = "black", linewidth = 0.7) +
  annotate("text", x = tau_seed + 0.06, y = Inf, vjust = 1.6,
           label = sprintf("Threshold\n(%.2f)", tau_seed),
           size = 3.4, hjust = 0) +
  labs(x = expression("Baseline LDL (" * Y[0] * ")"), y = "Density") +
  thesis_theme
ggsave(file.path(OUTPUT_DIR, "Fig01_Histogram_Baseline_Y0.png"),
       fig01, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig01_Histogram_Baseline_Y0.png\n")

#Fig02: Y0 distribution by treatment status — both scenarios
#Left panel: Randomised assignment.  Right panel: Threshold assignment.
#The threshold vertical line (tau_90) is drawn only in the Threshold panel.

fig02_df <- rbind(
  data.frame(Y0 = dat_seed$Y0,
             Status   = factor(ifelse(D_rand_seed == 1,
                                      "Treated (D = 1)", "Untreated (D = 0)"),
                               levels = c("Treated (D = 1)", "Untreated (D = 0)")),
             Scenario = factor("Randomised",
                               levels = c("Randomised", "Threshold"))),
  data.frame(Y0 = dat_seed$Y0,
             Status   = factor(ifelse(D_thr_seed == 1,
                                      "Treated (D = 1)", "Untreated (D = 0)"),
                               levels = c("Treated (D = 1)", "Untreated (D = 0)")),
             Scenario = factor("Threshold",
                               levels = c("Randomised", "Threshold")))
)

#Threshold line data: appears only in the Threshold facet
fig02_vline <- data.frame(
  Scenario = factor("Threshold", levels = c("Randomised", "Threshold")),
  tau      = tau_seed
)

fig02 <- ggplot(fig02_df, aes(x = Y0, fill = Status, colour = Status)) +
  geom_density(alpha = 0.40, linewidth = 0.7) +
  geom_vline(data = fig02_vline,
             aes(xintercept = tau),
             linetype = "dashed", colour = "black", linewidth = 0.7,
             inherit.aes = FALSE) +
  scale_fill_manual(values   = STRATUM_COLORS, name = "Treatment status") +
  scale_colour_manual(values = STRATUM_COLORS, guide = "none") +
  facet_wrap(~ Scenario) +
  labs(x = expression("Baseline LDL (" * Y[0] * ")"), y = "Density") +
  thesis_theme + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "Fig02_Y0_by_TreatmentStatus.png"),
       fig02, width = 10, height = 5, dpi = 300, bg = "white")

cat("  Saved: Fig02_Y0_by_TreatmentStatus.png\n")

#Fig03: Empirical variance decomposition of Y0 — two scenario panels
#Both panels are identical because Y0 is generated before treatment is
#assigned; the two-panel layout makes this scenario-independence explicit.

fig03 <- ggplot(vd_y0_bar,
                aes(x = Stratum, y = Pct, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55) +
  geom_text(aes(y = ypos, label = label),
            size = 3.4, colour = "white", fontface = "bold") +
  scale_fill_manual(values = VARDECOMP_Y0, name = "Component") +
  facet_wrap(~ scenario) +
  labs(x = NULL,
       y = expression("Percentage of Var(" * Y[0] * ") (%)")) +
  thesis_theme +
  theme(legend.position  = "bottom",
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())
ggsave(file.path(OUTPUT_DIR, "Fig03_VarDecomp_Y0.png"),
       fig03, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig03_VarDecomp_Y0.png\n")

#Fig04: Empirical variance decomposition of Y1
#x: Full Sample (Y0 baseline reference), Treated (D=1), Untreated (D=0).
#Facets: Randomised | Threshold.
#"Genetic" combines the stable G0 effect and the interaction term.

fig04 <- ggplot(vd_y1_bar,
                aes(x = Stratum, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = ypos, label = label),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = VARDECOMP_Y1, name = "Component") +
  facet_wrap(~ scenario) +
  labs(x = NULL,
       y = expression("Percentage of Var(" * Y[1] * ") (%)")) +
  thesis_theme +
  theme(legend.position = "bottom",
        axis.text.x     = element_text(size = 9))
ggsave(file.path(OUTPUT_DIR, "Fig04_VarDecomp_Y1.png"),
       fig04, width = 10, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig04_VarDecomp_Y1.png\n")

#Fig05: Empirical variance decomposition of DeltaY
#x: Full Sample (full-sample DeltaY), Treated (D=1), Untreated (D=0).
#Facets: Randomised | Threshold.
#Confounder = 0 for Treated and Untreated because beta_U*U cancels in
#Y1 - Y0; it appears in Full Sample at near-zero and is retained in the
#legend to make this cancellation visible.

fig05 <- ggplot(vd_dy_bar,
                aes(x = Stratum, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = ypos, label = label),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = VARDECOMP_DY, name = "Component") +
  facet_wrap(~ scenario) +
  labs(x = NULL,
       y = expression("Percentage of Var(" * Delta * Y * ") (%)")) +
  thesis_theme +
  theme(legend.position = "bottom",
        axis.text.x     = element_text(size = 9))
ggsave(file.path(OUTPUT_DIR, "Fig05_VarDecomp_DeltaY.png"),
       fig05, width = 10, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig05_VarDecomp_DeltaY.png\n")

#Fig06: Scatter Y0 vs Y1, both scenarios

fig06_df <- rbind(
  data.frame(Y0 = dat_seed$Y0, Y1 = Y1_thr_seed,
             Scenario = "Threshold",
             Status   = ifelse(D_thr_seed == 1, "Treated","Untreated")),
  data.frame(Y0 = dat_seed$Y0, Y1 = Y1_rand_seed,
             Scenario = "Randomised",
             Status   = ifelse(D_rand_seed == 1, "Treated","Untreated"))
) %>%
  mutate(Scenario = factor(Scenario, levels = c("Randomised","Threshold")))

fig06 <- ggplot(fig06_df %>% slice_sample(n = 4000),
                aes(x = Y0, y = Y1, colour = Status)) +
  geom_point(alpha = 0.35, size = 0.6) +
  geom_smooth(aes(fill = Status), method = "lm", se = TRUE,
              linewidth = 0.8, alpha = 0.20) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("Treated" = "#D6604D", "Untreated" = "#4A7FA5"),
                      name = "Treatment status") +
  scale_fill_manual(values = c("Treated" = "#D6604D", "Untreated" = "#4A7FA5"),
                    guide = "none") +
  facet_wrap(~ Scenario) +
  labs(x = expression("Baseline LDL (" * Y[0] * ")"),
       y = expression("Follow-up LDL (" * Y[1] * ")")) +
  thesis_theme + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "Fig06_Scatter_Y0vsY1.png"),
       fig06, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig06_Scatter_Y0vsY1.png\n")

#Fig07: Scatter Y0 vs DeltaY, both scenarios (combined)
#Randomised is the left panel, Threshold the right, matching Fig06 layout.
#The two previously separate figures (Threshold-only and Randomised-only) are
#merged here; Fig08 is superseded.

fig07_df <- rbind(
  data.frame(
    Y0       = dat_seed$Y0,
    DY       = Y1_rand_seed - dat_seed$Y0,
    Status   = factor(ifelse(D_rand_seed == 1, "Treated", "Untreated"),
                      levels = c("Treated", "Untreated")),
    Scenario = factor("Randomised", levels = c("Randomised", "Threshold"))
  ),
  data.frame(
    Y0       = dat_seed$Y0,
    DY       = Y1_thr_seed - dat_seed$Y0,
    Status   = factor(ifelse(D_thr_seed == 1, "Treated", "Untreated"),
                      levels = c("Treated", "Untreated")),
    Scenario = factor("Threshold", levels = c("Randomised", "Threshold"))
  )
)

fig07 <- ggplot(fig07_df %>% slice_sample(n = 8000),
                aes(x = Y0, y = DY, colour = Status)) +
  geom_point(alpha = 0.35, size = 0.6) +
  geom_smooth(aes(fill = Status), method = "lm", se = TRUE,
              linewidth = 0.8, alpha = 0.20) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c("Treated" = "#D6604D", "Untreated" = "#4A7FA5"),
                      name = "Treatment status") +
  scale_fill_manual(values = c("Treated" = "#D6604D", "Untreated" = "#4A7FA5"),
                    guide = "none") +
  facet_wrap(~ Scenario) +
  labs(x = expression("Baseline LDL (" * Y[0] * ")"),
       y = expression(Delta * Y ~ "=" ~ Y[1] - Y[0])) +
  thesis_theme + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "Fig07_Scatter_Y0vsDY_Combined.png"),
       fig07, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig07_Scatter_Y0vsDY_Combined.png\n\n")

#Fig08: superseded by Fig07 combined panel
#fig08_df and fig08 are no longer generated; the Randomised scenario is
#now the left panel of Fig07_Scatter_Y0vsDY_Combined.png.


# =============================================================================
# SECTION 10: GENOTYPE-MEAN FIGURES (Fig09 – Fig12)
# =============================================================================
# For each SNP class, compute grand mean(DeltaY | genotype = 0/1/2) within
# each stratum across 1000 iterations, then plot as violin + boxplot.
# Uses the same seed offset (i + 10000L) as the main loop to ensure
# independent datasets for this set of figures.

cat("Generating genotype-mean DeltaY figures (Fig09 - Fig12) ...\n")

GENO_SPECS <- list(
  list(label = "G0",   idxs = idx_G0,   fignum = "Fig09"),
  list(label = "GD",   idxs = idx_GD,   fignum = "Fig10"),
  list(label = "G0D",  idxs = idx_G0D,  fignum = "Fig11"),
  list(label = "GNull",idxs = idx_Null, fignum = "Fig12")
)

for (sp in GENO_SPECS) {
  cat(sprintf("  Building %s_GenotypeMean_DeltaY_%s ...\n", sp$fignum, sp$label))
  
  geno_rows <- vector("list", N_ITERATIONS * 2L * 2L * 3L)
  row_k     <- 1L
  
  for (i in seq_len(N_ITERATIONS)) {
    dat_i <- simulate_data(seed = i + 10000L)
    
    for (sc in c("Rand", "Threshold")) {
      D_i  <- assign_treatment(dat_i$Y0, sc)
      Y1_i <- generate_Y1(dat_i$G, dat_i$U, dat_i$Y0, D_i)
      dY_i <- Y1_i - dat_i$Y0
      
      for (stratum_lbl in c(STRAT_TRT, STRAT_UNT)) {
        idx_s <- if (stratum_lbl == STRAT_TRT) which(D_i == 1L) else which(D_i == 0L)
        if (length(idx_s) < 10L) next
        dY_s <- dY_i[idx_s]
        
        snp_means <- vapply(sp$idxs, function(k) {
          g_s <- dat_i$G[idx_s, k]
          c(mean(dY_s[g_s == 0L], na.rm = TRUE),
            mean(dY_s[g_s == 1L], na.rm = TRUE),
            mean(dY_s[g_s == 2L], na.rm = TRUE))
        }, numeric(3L))
        
        gm       <- rowMeans(snp_means, na.rm = TRUE)
        sc_label <- if (sc == "Rand") "Randomised" else "Threshold"
        
        for (g_idx in 1:3) {
          geno_rows[[row_k]] <- data.frame(
            Iter     = i,
            Scenario = sc_label,
            Stratum  = stratum_lbl,
            Genotype = c("g = 0", "g = 1", "g = 2")[g_idx],
            MeanDY   = gm[g_idx]
          )
          row_k <- row_k + 1L
        }
      }
    }
  }
  
  geno_df <- bind_rows(geno_rows) %>%
    filter(!is.na(MeanDY)) %>%
    mutate(
      Scenario = factor(Scenario, levels = c("Randomised","Threshold")),
      Stratum  = factor(Stratum,  levels = c(STRAT_TRT, STRAT_UNT)),
      Genotype = factor(Genotype, levels = c("g = 0","g = 1","g = 2"))
    )
  
  pad <- 0.05
  trt_d  <- filter(geno_df, Stratum == STRAT_TRT)
  unt_d  <- filter(geno_df, Stratum == STRAT_UNT)
  trt_q  <- quantile(trt_d$MeanDY, c(0.002, 0.998), na.rm = TRUE)
  unt_q  <- quantile(unt_d$MeanDY, c(0.002, 0.998), na.rm = TRUE)
  trt_lim <- c(trt_q[1] - pad * diff(trt_q), trt_q[2] + pad * diff(trt_q))
  unt_lim <- c(unt_q[1] - pad * diff(unt_q), unt_q[2] + pad * diff(unt_q))
  
  base_aes <- aes(x = Genotype, y = MeanDY, fill = Scenario, colour = Scenario)
  common_layers <- list(
    scale_fill_manual(values   = SCENARIO_COLS_2, name = "Scenario"),
    scale_colour_manual(values = SCENARIO_COLS_2, guide = "none"),
    thesis_theme,
    theme(legend.position = "none")
  )
  
  p_trt <- ggplot(trt_d, base_aes) +
    geom_violin(position = position_dodge(0.8), alpha = 0.35,
                linewidth = 0.55, scale = "width", trim = TRUE) +
    geom_boxplot(position = position_dodge(0.8), width = 0.18,
                 alpha = 0.80, outlier.shape = NA, linewidth = 0.6,
                 colour = "grey25") +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.55) +
    coord_cartesian(ylim = trt_lim) +
    labs(x = sprintf("%s genotype dosage", sp$label),
         y = expression("Grand mean " * Delta * Y ~ "(treated, D = 1)")) +
    common_layers
  
  p_unt <- ggplot(unt_d, base_aes) +
    geom_violin(position = position_dodge(0.8), alpha = 0.35,
                linewidth = 0.55, scale = "width", trim = TRUE) +
    geom_boxplot(position = position_dodge(0.8), width = 0.18,
                 alpha = 0.80, outlier.shape = NA, linewidth = 0.6,
                 colour = "grey25") +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.55) +
    coord_cartesian(ylim = unt_lim) +
    labs(x = sprintf("%s genotype dosage", sp$label),
         y = expression("Grand mean " * Delta * Y ~ "(untreated, D = 0)")) +
    common_layers
  
  get_legend <- function(p) {
    g   <- ggplot_gtable(ggplot_build(
      p + theme(legend.position = "bottom")))
    pos <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(pos) == 0L) return(NULL)
    g$grobs[[pos]]
  }
  leg <- get_legend(p_trt + theme(legend.position = "bottom"))
  
  fig_out <- gridExtra::arrangeGrob(
    gridExtra::arrangeGrob(p_trt, p_unt, ncol = 2),
    leg, nrow = 2, heights = c(10, 1)
  )
  fname <- sprintf("%s_GenotypeMean_DeltaY_%s.png", sp$fignum, sp$label)
  ggsave(file.path(OUTPUT_DIR, fname),
         fig_out, width = 10, height = 6, dpi = 300, bg = "white")
  cat(sprintf("    Saved: %s\n", fname))
}
cat("\n")


# =============================================================================
# SECTION 11: OBJECTIVE 1 FIGURES (Fig13 – Fig14)
# =============================================================================

cat("Generating Objective 1 figures (Fig13 - Fig14) ...\n")

#Fig13: M1 coefficient distributions, all strata

fig13 <- make_boxplot(fl_M1_all)
ggsave(file.path(OUTPUT_DIR, "Fig13_M1_Boxplot.png"),
       fig13, width = 10, height = 7, dpi = 300, bg = "white")
cat("  Saved: Fig13_M1_Boxplot.png\n")

#Fig14: M1 QQ plot, all SNP classes, all strata × scenarios

fig14 <- make_qq_allclass(fl_M1_all, annot_sz = FS_ANNOT_LG)
ggsave(file.path(OUTPUT_DIR, "Fig14_M1_QQ_Large.png"),
       fig14, width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig14_M1_QQ_Large.png\n\n")


# =============================================================================
# SECTION 12: OBJECTIVE 2 FIGURES (Fig15, Fig17, Fig20)
# =============================================================================

cat("Generating Objective 2 figures (Fig15, Fig17, Fig20) ...\n")

# Named list for the model-comparison QQ function
fl_by_model <- list(
  "M1"  = fl_M1_all,
  "M2"  = fl_M2_all,
  "Dub" = fl_Dubt %>% mutate(stratum = factor(stratum, levels = c(STRAT_TRT))),
  "SH"  = fl_SHt  %>% mutate(stratum = factor(stratum, levels = c(STRAT_TRT)))
)

#Fig15: G0 QQ comparing M1 and M2, both strata 

fig15 <- make_qq_model_comparison(fl_by_model, models = c("M1","M2"))
ggsave(file.path(OUTPUT_DIR, "Fig15_G0_QQ_M1vsM2.png"),
       fig15, width = 9, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig15_G0_QQ_M1vsM2.png\n")

#Fig17: M2 QQ, all SNP classes, both strata 

fig17 <- make_qq_allclass(fl_M2_all)
ggsave(file.path(OUTPUT_DIR, "Fig17_M2_QQ.png"),
       fig17, width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig17_M2_QQ.png\n")

#Fig20: M2 coefficient distributions, both strata 

fig20 <- make_boxplot(fl_M2_all)
ggsave(file.path(OUTPUT_DIR, "Fig20_M2_Boxplot.png"),
       fig20, width = 10, height = 7, dpi = 300, bg = "white")
cat("  Saved: Fig20_M2_Boxplot.png\n\n")


# =============================================================================
# SECTION 13: OBJECTIVE 3 FIGURES (Fig16, Fig18, Fig19, Fig21, Fig22)
# =============================================================================

cat("Generating Objective 3 figures (Fig16, Fig18, Fig19, Fig21, Fig22) ...\n")

#Fig16: G0 QQ all four approaches, treated + untreated rows
#Dub and SH are only applicable in the treated stratum.
#The untreated row shows M1 and M2 for reference.

fl_obj3_all <- bind_rows(
  # Treated: all four models
  fl_by_model[["M1"]]  %>% filter(stratum == STRAT_TRT),
  fl_by_model[["M2"]]  %>% filter(stratum == STRAT_TRT),
  fl_Dubt,
  fl_SHt,
  # Untreated: M1 and M2 only
  fl_by_model[["M1"]]  %>% filter(stratum == STRAT_UNT),
  fl_by_model[["M2"]]  %>% filter(stratum == STRAT_UNT)
) %>% mutate(
  model_f = factor(model, levels = c("M1","M2","Dub","SH")),
  stratum  = factor(stratum, levels = c(STRAT_TRT, STRAT_UNT))
)

fig16 <- make_qq_model_comparison(
  fl_named_list = list(
    "M1"  = fl_obj3_all %>% filter(model == "M1"),
    "M2"  = fl_obj3_all %>% filter(model == "M2"),
    "Dub" = fl_obj3_all %>% filter(model == "Dub"),
    "SH"  = fl_obj3_all %>% filter(model == "SH")
  ),
  models = c("M1","M2","Dub","SH")
)
ggsave(file.path(OUTPUT_DIR, "Fig16_G0_QQ_M1vsM2vsDubvsSH.png"),
       fig16, width = 9, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig16_G0_QQ_M1vsM2vsDubvsSH.png\n")

#Fig18: Dudbridge QQ, all SNP classes, treated stratum only

fig18 <- make_qq_allclass(fl_Dubt %>%
                            mutate(stratum = STRAT_TRT) %>%
                            add_stratum_lab(),
                          treated_only = TRUE)
ggsave(file.path(OUTPUT_DIR, "Fig18_Dub_QQ.png"),
       fig18, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig18_Dub_QQ.png\n")

#Fig19: Slope-Hunter QQ, all SNP classes, treated stratum only

fig19 <- make_qq_allclass(fl_SHt %>%
                            mutate(stratum = STRAT_TRT) %>%
                            add_stratum_lab(),
                          treated_only = TRUE)
ggsave(file.path(OUTPUT_DIR, "Fig19_SH_QQ.png"),
       fig19, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig19_SH_QQ.png\n")

#Fig21: Dudbridge boxplot, treated stratum

fig21 <- make_boxplot(fl_Dubt, treated_only = TRUE)
ggsave(file.path(OUTPUT_DIR, "Fig21_Dub_Boxplot.png"),
       fig21, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig21_Dub_Boxplot.png\n")

#Fig22: Slope-Hunter boxplot, treated stratum

fig22 <- make_boxplot(fl_SHt, treated_only = TRUE)
ggsave(file.path(OUTPUT_DIR, "Fig22_SH_Boxplot.png"),
       fig22, width = 9, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig22_SH_Boxplot.png\n\n")


# =============================================================================
# SECTION 14: SIGNED BIAS TABLES
# =============================================================================
#Convention
# Treated stratum   : bias = beta_hat - gamma_true
#                     (gamma_true from snp_meta; 0 for G0/Null, gamma_fixed for GD/G0D)
# Untreated stratum : bias = beta_hat - 0
#                     (true pharmacogenetic effect = 0 for D = 0)
# =============================================================================

cat("Computing signed bias tables ...\n")

fmt_val <- function(x) {
  ifelse(abs(x) < 0.001,
         formatC(x, format = "e", digits = 2),
         sprintf("%.4f", x))
}


compute_bias <- function(fl, model_label, strat_val) {
  true_gamma <- if (strat_val == STRAT_TRT) {
    snp_meta %>% select(SNP_Type, gamma_true)
  } else {
    snp_meta %>% select(SNP_Type) %>% mutate(gamma_true = 0)
  }
  fl %>%
    filter(stratum == strat_val, !is.na(beta)) %>%
    left_join(true_gamma, by = "SNP_Type") %>%
    mutate(signed_bias = beta - gamma_true) %>%
    group_by(SNP_Type, scenario) %>%
    summarise(
      mean_b = mean(signed_bias, na.rm = TRUE),
      se_b   = sd(signed_bias,  na.rm = TRUE) / sqrt(sum(!is.na(signed_bias))),
      .groups = "drop"
    ) %>%
    mutate(
      cell  = paste0(fmt_val(mean_b), " (", fmt_val(se_b), ")"),
      Model = model_label,
      Stratum = strat_val
    )
}

pivot_wide <- function(bias_df) {
  bias_df %>%
    mutate(
      sc_short = case_when(
        grepl("Rand", scenario) ~ "Rand",
        TRUE                    ~ "Thr"),
      col = paste(sc_short, "Treated", sep = "_")
    ) %>%
    select(SNP_Type, Model, col, cell) %>%
    pivot_wider(names_from = col, values_from = cell) %>%
    mutate(SNP_Type = factor(SNP_Type, levels = c("G0","GD","G0D","Null"))) %>%
    arrange(SNP_Type, Model)
}

#Table 1: M1 vs M2, both strata

bias_m1m2 <- bind_rows(
  compute_bias(fl_M1_all %>% mutate(stratum = stratum),  "M1", STRAT_TRT),
  compute_bias(fl_M1_all %>% mutate(stratum = stratum),  "M1", STRAT_UNT),
  compute_bias(fl_M2_all %>% mutate(stratum = stratum),  "M2", STRAT_TRT),
  compute_bias(fl_M2_all %>% mutate(stratum = stratum),  "M2", STRAT_UNT)
) %>%
  mutate(
    sc_short = case_when(grepl("Rand", scenario) ~ "Rand", TRUE ~ "Thr"),
    st_short = case_when(Stratum == STRAT_TRT ~ "Treated", TRUE ~ "Untreated"),
    col = paste0(sc_short, "_", st_short)
  ) %>%
  select(SNP_Type, Model, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  mutate(SNP_Type = factor(SNP_Type, levels = c("G0","GD","G0D","Null"))) %>%
  arrange(SNP_Type, Model)

write.csv(bias_m1m2,
          file.path(OUTPUT_DIR, "Table_Bias_M1_M2.csv"),
          row.names = FALSE)
cat("  Saved: Table_Bias_M1_M2.csv\n")
print(bias_m1m2, n = Inf)

#Table 2: M1, SH, M2, Dub — treated stratum only

bias_obj3 <- bind_rows(
  compute_bias(fl_M1t, "M1",  STRAT_TRT),
  compute_bias(fl_SHt, "SH",  STRAT_TRT),
  compute_bias(fl_M2t, "M2",  STRAT_TRT),
  compute_bias(fl_Dubt,"Dub", STRAT_TRT)
) %>%
  mutate(
    sc_short = case_when(grepl("Rand", scenario) ~ "Rand", TRUE ~ "Thr"),
    col = paste0(sc_short, "_Treated")
  ) %>%
  select(SNP_Type, Model, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  mutate(SNP_Type = factor(SNP_Type, levels = c("G0","GD","G0D","Null")),
         Model    = factor(Model,    levels = c("M1","SH","M2","Dub"))) %>%
  arrange(SNP_Type, Model)

write.csv(bias_obj3,
          file.path(OUTPUT_DIR, "Table_Bias_Obj3.csv"),
          row.names = FALSE)
cat("  Saved: Table_Bias_Obj3.csv\n\n")
print(bias_obj3, n = Inf)


# =============================================================================
# SECTION 15: SAVE WORKSPACE
# =============================================================================

if (SAVE_WORKSPACE) {
  ws_path <- file.path(OUTPUT_DIR, "simulation_workspace.RData")
  save.image(file = ws_path)
  cat(sprintf("\nWorkspace saved to: %s\n", ws_path))
}

cat(strrep("=", 65), "\n")
cat(sprintf("All outputs written to: %s/\n", OUTPUT_DIR))
cat(strrep("=", 65), "\n")



summary(full_res$b_sh[full_res$scenario == "Threshold"])

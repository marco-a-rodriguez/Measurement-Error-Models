################################################################################
# R code to generate results of the section "Results: Estimates of movement
# scale and fractionation in three fish species" in "Measurement error models
# reveal the scale of consumer movements along an isoscape gradient"
#
# Reads data from MEMdataSalmonRasmussen2009.csv, MEMdataBluesharkBird2018.csv,
# or MEMdataCatsharkBird2018.csv
#
# Calls Stan programs MEMgaussian.stan, MEMlaplace.stan, or MEMstudent.stan
#
# Loads required functions from MEMfunctions.r
#
# M.A. Rodriguez
################################################################################

# Load required libraries and functions
source("MEMfunctions.r")

# Specify MEM model type used for parameter estimation:
# one of "gaussian", "laplace", or "student"
model_type = "gaussian"
if (model_type == "gaussian") {
  stan_model = "MEMgaussian.stan"
  f_sd = function(x) ftgaussian_sig(x, nu_mcmc[i], xl, xu) * dunif(x, xl, xu)
} else if (model_type == "laplace") {
  stan_model = "MEMlaplace.stan"
  f_sd = function(x) fdtlaplace_sig(x, nu_mcmc[i], xl, xu) * dunif(x, xl, xu)
} else if (model_type == "student") {
  stan_model = "MEMstudent.stan"
  f_sd = function(x) fdtstudent_sig(x, nu_mcmc[i], xl, xu) * dunif(x, xl, xu)
}

# Specify data type: one of "salmon", "blueshark", or "catshark"
data_type = "salmon"
if (data_type == "salmon") {
  data_file = "MEMdataSalmonRasmussen2009.csv"
} else if (data_type == "blueshark") {
  data_file = "MEMdataBluesharkBird2018.csv"
} else if (data_type == "catshark") {
  data_file = "MEMdataCatsharkBird2018.csv"
}

# Read data and specify bounds of study region
dat = read.csv(data_file, head=TRUE)
yb = na.omit(dat$yb)
yc = na.omit(dat$yc)
xb = dat$x[!is.na(dat$yb)]
xc = dat$x[!is.na(dat$yc)]
xl = min(xb, xc)
xu = max(xb, xc)
Nb = length(na.omit(dat$yb))
Nc = length(na.omit(dat$yc))

# Scale and center x and y variables using baseline mean and sd
mean_yb = mean(yb)
sd_yb = sd(yb)
mean_xb = mean(xb)
sd_xb = sd(xb)
yb_sc = (yb - mean_yb) / sd_yb
yc_sc = (yc - mean_yb) / sd_yb
xb_sc = (xb - mean_xb) / sd_xb
xc_sc = (xc - mean_xb) / sd_xb
xl_sc = (xl - mean_xb) / sd_xb
xu_sc = (xu - mean_xb) / sd_xb

# Obtain OLS regression estimates for QAD and SAURM models
regybxb = lm(yb_sc ~ xb_sc)
regycxc = lm(yc_sc ~ xc_sc)
alphab = summary(regybxb)$coef[1, 1]
betab = summary(regybxb)$coef[2, 1]
alphac = summary(regycxc)$coef[1, 1]
betac = summary(regycxc)$coef[2, 1]
serror = summary(regybxb)$sigma
var_x = var(xb_sc)
var_qad_sc = (abs(betab / betac) - 1) * var_x
range_prop_saurm = (1 - betac / betab)
gamma_qad_sc = alphac - (alphab + betab * mean(xb_sc) * var_qad_sc / (var_x + var_qad_sc))
gamma_saurm_sc = alphac - mean(yb_sc) * (1 - betac / betab) - alphab * (betac / betab)

# Obtain Bayesian estimates for MEM models
stan_dat = list(xb_sc=xb_sc, yb_sc=yb_sc, Nb=Nb, xc_sc=xc_sc, yc_sc=yc_sc, Nc=Nc, xl_sc=xl_sc, xu_sc=xu_sc)
pars = c('alpha','beta','gamma','sigmab','nu','xstar_sc','sigmac')
iter = 100000
chains = 5
thin = 75
warmup = 25000
fit_mem = stan(file=stan_model, data=stan_dat,
               control=list(adapt_delta=0.99, adapt_gamma=0.05, max_treedepth=20),
               iter=iter, chains=chains, thin=thin, warmup=warmup, seed=73)

# Effective sample sizes and Gelman-Rubin Rhat diagnostics (parameter estimates are for scaled variables)
print(fit_mem, pars = pars, digits_summary = 3)

# Plot chain autocorrelation
quietgg(stan_ac(fit_mem, pars = c('alpha','beta','gamma','sigmab','nu','sigmac')))

# Extract model diagnostics
check_hmc_diagnostics(fit_mem)

# Obtain score values (CRPS and PPL) for model comparisons
# Observation vector
y_obs = c(yb_sc, yc_sc)
# Posterior predictive vector
fit = extract(fit_mem)
y_rep = t(fit$y_rep)

# CRPS
crps = mean(crps_sample(y_obs, y_rep))

# Posterior predictive loss with k, the relative weight of losses associated
# with fit and variance terms, set to 1 (see Gelfand & Ghosh, 1998)
k = 1
pplf = function(yo, yr) {
  mean_yr = apply(yr, 1, mean)
  var_yr = apply(yr, 1, var)
  L2_fit = sum((mean_yr - yo) ^ 2)
  L2_var = sum(var_yr)
  (k / (k + 1)) * L2_fit + L2_var
}
PPL = pplf(y_obs, y_rep)

# Point estimates of scale (nu) and fractionation (gamma) for QAD and SAURM models
# after back-transformation to original scales
lambda_hat = betac / betab
if(abs(lambda_hat) > 1) {var_qad_sc = 0}
nu_qad = sqrt(var_qad_sc) * sd_xb
nu_saurm = range_prop_saurm * (xu - xl) / sqrt(12)
if(lambda_hat < 0) {
  nu_saurm = 1 * (xu - xl) / sqrt(12)
} else if (lambda_hat > 1) {
  nu_saurm = 0 * (xu - xl) / sqrt(12)
}

gamma_qad = gamma_qad_sc * sd_yb
gamma_saurm = gamma_saurm_sc * sd_yb

# Medians and highest posterior density (HPD) 95% credible intervals for MEM model parameters
# after back-transformation to original scales

# Expected value for standard deviation of truncated density averaged over
# entire study range, omega
nu_mcmc = extract(fit_mem, pars = c('nu'))$nu * sd_xb
omega = numeric(length(nu_mcmc))
for(i in 1:length(nu_mcmc)) {omega[i] = hcubature(f_sd, xl, xu)$integral}
omega = as.mcmc(omega)

sprintf(
  "Expected value for scale parameter omega --- Median: %.3f; 95%% CI: %.3f - %.3f",
  median(omega),
  HPDinterval(omega)[1],
  HPDinterval(omega)[2]
)

# Fractionation parameter gamma
gamma_mcmc = as.mcmc(as.numeric(extract(fit_mem, pars = c('gamma'))$gamma)) * sd_yb
sprintf(
  "Fractionation parameter gamma --- Median: %.2f; 95%% CI: %.2f - %.2f",
  median(gamma_mcmc),
  as.vector(HPDinterval(gamma_mcmc))[1],
  as.vector(HPDinterval(gamma_mcmc))[2]
)

# Baseline regression slope parameter beta
beta_mcmc = as.mcmc(as.numeric(extract(fit_mem, pars = c('beta'))$beta)) * sd_yb / sd_xb
sprintf(
  "Baseline regression slope parameter beta --- Median: %.5f; 95%% CI: %.5f - %.5f",
  median(beta_mcmc),
  as.vector(HPDinterval(beta_mcmc))[1],
  as.vector(HPDinterval(beta_mcmc))[2]
)

# Baseline regression intercept parameter alpha
alpha_mcmc = mean_yb + as.mcmc(as.numeric(extract(fit_mem, pars = c('alpha'))$alpha)) * sd_yb -
  mean_xb * as.mcmc(as.numeric(extract(fit_mem, pars = c('beta'))$beta)) * sd_yb / sd_xb
sprintf(
  "Baseline regression intercept parameter alpha --- Median: %.2f; 95%% CI: %.2f - %.2f",
  median(alpha_mcmc),
  as.vector(HPDinterval(alpha_mcmc))[1],
  as.vector(HPDinterval(alpha_mcmc))[2]
)

# Baseline regression SE parameter sigmab
sigmab_mcmc = as.mcmc(as.numeric(extract(fit_mem, pars = c('sigmab'))$sigmab)) * sd_yb
sprintf(
  "Baseline regression SE parameter sigmab --- Median: %.2f; 95%% CI: %.2f - %.2f",
  median(sigmab_mcmc),
  as.vector(HPDinterval(sigmab_mcmc))[1],
  as.vector(HPDinterval(sigmab_mcmc))[2]
)

# Baseline regression SE parameter sigmac
sigmac_mcmc = as.mcmc(as.numeric(extract(fit_mem, pars = c('sigmac'))$sigmac)) * sd_yb
sprintf(
  "Baseline regression SE parameter sigmac --- Median: %.2f; 95%% CI: %.2f - %.2f",
  median(sigmac_mcmc),
  as.vector(HPDinterval(sigmac_mcmc))[1],
  as.vector(HPDinterval(sigmac_mcmc))[2]
)


################################################################################
# R sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
#
# Matrix products: default
# BLAS/LAPACK: /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/imkl/2018.3.222/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so
#
# locale:
# LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C
# LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8
# LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8
# LC_PAPER=en_CA.UTF-8       LC_NAME=C
# LC_ADDRESS=C               LC_TELEPHONE=C
# LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
# stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
# cubature_2.0.4     crch_1.0-4         extraDistr_1.8.11  scoringRules_1.0.0
# coda_0.19-3        rstan_2.19.3       ggplot2_3.3.0      StanHeaders_2.19.2
#
# loaded via a namespace (and not attached):
# Rcpp_1.0.4.6       pillar_1.4.3       compiler_3.6.1     prettyunits_1.1.1
# tools_3.6.1        pkgbuild_1.0.7     lifecycle_0.2.0    tibble_3.0.1
# gtable_0.3.0       lattice_0.20-41    pkgconfig_2.0.3    rlang_0.4.5
# cli_2.0.2          parallel_3.6.1     xfun_0.13          loo_2.2.0
# gridExtra_2.3      withr_2.2.0        knitr_1.28         vctrs_0.2.4
# stats4_3.6.1       grid_3.6.1         glue_1.4.0         inline_0.3.15
# R6_2.4.1           processx_3.4.2     fansi_0.4.1        Formula_1.2-3
# callr_3.4.3        magrittr_1.5       scales_1.1.0       ps_1.3.2
# ellipsis_0.3.0     matrixStats_0.56.0 MASS_7.3-51.6      assertthat_0.2.1
# colorspace_1.4-1   sandwich_2.5-1     munsell_0.5.0      crayon_1.3.4
# zoo_1.8-7

################################################################################
# R code to generate results of the section "Results: Performance of estimators
# in simulated scenarios" in "Measurement error models reveal the scale of
# consumer movements along an isoscape gradient"
#
# Calls Stan program MEMsims.stan
#
# Loads required functions from MEMfunctions.r
#
# M.A. Rodriguez
################################################################################

# Load required libraries and functions and set seed for R generators
source("MEMfunctions.r")
set.seed(173)

# Select simulation scenario: sample size (N) and movement range as a proportion of
# the study range (S)
N = 300
S = 0.05

# Set isoscape parameters
Nb = Nc = N
r2 = 0.7
range_yb = 12
sdyb = range_yb / sqrt(12)
sigmab = sqrt(1 - r2) * sdyb * sqrt(Nb / (Nb - 1))
sigmac = 0.5 * sigmab

xl = 0
xu = 1
sdxb = (xu - xl) / sqrt(12)
alpha = -30
beta = sqrt(r2) * sdyb / sdxb
gamma = 1

# Determine number of simulations (iters) and initialize variables
iters = 250
var_qad_sc = numeric(iters)
range_prop_saurm = numeric(iters)
nu_median_sc = numeric(iters)
mean_xb = numeric(iters)
sd_xb = numeric(iters)
mean_yb = numeric(iters)
sd_yb = numeric(iters)
gamma_qad_sc = numeric(iters)
gamma_saurm_sc = numeric(iters)
gamma_median_sc = numeric(iters)
alphab = numeric(iters)
alphac = numeric(iters)
betab = numeric(iters)
betac = numeric(iters)
alpha_median_sc = numeric(iters)
beta_median_sc = numeric(iters)
sigmab_median_sc = numeric(iters)
sigmac_median_sc = numeric(iters)

pars = c('alpha','beta','gamma','sigmab','nu','sigmac')
npars = length(pars)
neff = matrix(0, iters, npars)
rhat = matrix(0, iters, npars)

# Enter simulation loop
for(i in 1:iters) {

# Generate baseline and consumer observations
xb = runif(Nb, xl, xu)
yb = rnorm(Nb, alpha + beta * xb, sigmab)
xstar = runif(Nc, xl, xu)
yc = rnorm(Nc, alpha + gamma + beta * xstar, sigmac)

# Generate random movements from specified movement model
# ("gaussian", "laplace", or "uniform")
movement_model = "gaussian"

xc = numeric(Nc)
if (movement_model == "gaussian") {
  nu = (xu - xl) * S
  for(j in 1:Nc) {
    xc[j] = rtnorm(1, xstar[j], nu, xl, xu)
  }
} else if (movement_model == "laplace") {
  nu = (xu - xl) * S / sqrt(2)
  for(j in 1:Nc) {
    xc[j] = rtlaplace(1, xstar[j], nu, xl, xu)
  }
} else if (movement_model == "uniform") {
  nu = (xu - xl) * S / sqrt(12)
  delta = (xu - xl) * S / 2
  for(j in 1:Nc){
    xc[j] = runif(1, max(xstar[j] - delta, xl), min(xstar[j] + delta, xu))
  }
}

# Scale and center x and y variables using baseline mean and sd
mean_yb[i] = mean(yb)
sd_yb[i] = sd(yb)
mean_xb[i] = mean(xb)
sd_xb[i] = sd(xb)
yb_sc = (yb - mean_yb[i]) / sd_yb[i]
yc_sc = (yc - mean_yb[i]) / sd_yb[i]
xb_sc = (xb - mean_xb[i]) / sd_xb[i]
xc_sc = (xc - mean_xb[i]) / sd_xb[i]
xl_sc = (xl - mean_xb[i]) / sd_xb[i]
xu_sc = (xu - mean_xb[i]) / sd_xb[i]

# OLS Regression estimates for QAD and SAURM models
regybxb = lm(yb_sc ~ xb_sc)
regycxc = lm(yc_sc ~ xc_sc)
alphab[i] = summary(regybxb)$coef[1, 1]
betab[i] = summary(regybxb)$coef[2, 1]
alphac[i] = summary(regycxc)$coef[1, 1]
betac[i] = summary(regycxc)$coef[2, 1]
var_x = var(xb_sc)
var_qad_sc[i] = (abs(betab[i] / betac[i]) - 1) * var_x
range_prop_saurm[i] = (1 - betac[i] / betab[i])
gamma_qad_sc[i] = alphac[i] - (alphab[i] + betab[i] * mean(xb_sc) *
               var_qad_sc[i] / (var_x + var_qad_sc[i]))
gamma_saurm_sc[i] = alphac[i] - mean(yb_sc) * (1 - betac[i] / betab[i]) -
                 alphab[i] * (betac[i] / betab[i])

# mcmcian estimation for MEM models
stan_dat = list(xb_sc=xb_sc, yb_sc=yb_sc, Nb=Nb, xc_sc=xc_sc, yc_sc=yc_sc,
                Nc=Nc, xl_sc=xl_sc, xu_sc=xu_sc)
samples = 20000
chains = 1
thin = 15
warmup = 5000
fit_mem = stan(file="MEMsims.stan", data=stan_dat,
               control=list(adapt_delta=0.99, max_treedepth=20),
               init="random",init_r=0.5,iter=samples, chains=chains, thin=thin, warmup=warmup)

# Extract MCMC samples from the posterior distribution and obtain medians
nu_mcmc_sc = extract(fit_mem, pars=c('nu'))$nu
alpha_mcmc_sc = extract(fit_mem, pars=c('alpha'))$alpha
beta_mcmc_sc = extract(fit_mem, pars=c('beta'))$beta
gamma_mcmc_sc = extract(fit_mem, pars=c('gamma'))$gamma
sigmab_mcmc_sc = extract(fit_mem, pars=c('sigmab'))$sigmab
sigmac_mcmc_sc = extract(fit_mem, pars=c('sigmac'))$sigmac
nu_median_sc[i] = median(nu_mcmc_sc)
alpha_median_sc[i] = median(alpha_mcmc_sc)
beta_median_sc[i] = median(beta_mcmc_sc)
gamma_median_sc[i] = median(gamma_mcmc_sc)
sigmab_median_sc[i] = median(sigmab_mcmc_sc)
sigmac_median_sc[i] = median(sigmac_mcmc_sc)

# Obtain effective sample size and Gelman-Rubin diagnostics
neff[i, ] = summary(fit_mem)$summary[1:npars, "n_eff"]
rhat[i, ] = summary(fit_mem)$summary[1:npars, "Rhat"]

}
# End of simulation loop

# Model diagnostics
pairs(fit_mem, pars = c('alpha','beta','gamma','sigmab','nu','sigmac'))

# Examination of pairs() plots for individual chains shows that
# divergent transition warnings can occur when the chain is exploring
# very small values of nu, suggesting that the NUTS sampler is not fully
# exploring the parameter space near those values. This can lead to bias
# in estimates of nu, especially at smaller sample sizes
# (see "Discussion: Performance of estimators in simulated scenarios", and
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup)

# Effective sample size (neff) and Gelman-Rubin diagnostics (rhat)
# for MEM parameters
colnames(neff) = dimnames(summary(fit_mem)$summary)[[1]][1:6]
colnames(rhat) = dimnames(summary(fit_mem)$summary)[[1]][1:6]
summary(neff)
summary(rhat)


# Scale (nu) and fractionation (gamma) for QAD and SAURM models
# Set unfeasible QAD and SAURM nu estimates to zero
lambda_hat = betac / betab

for(i in 1:length(lambda_hat)) {
  if(abs(lambda_hat[i]) > 1) {
    var_qad_sc[i] = 0
  }
}
nu_qad = sqrt(var_qad_sc) * sd_xb

nu_saurm = range_prop_saurm * (xu - xl) / sqrt(12)
for(i in 1:length(lambda_hat)) {
  if(lambda_hat[i] < 0) {
    nu_saurm[i] = 1 * (xu - xl) / sqrt(12)
  } else if (lambda_hat[i] > 1) {
    nu_saurm[i] = 0 * (xu - xl) / sqrt(12)
  }
}

# Fractionation parameter gamma
gamma_qad = gamma_qad_sc * sd_yb
gamma_saurm = gamma_saurm_sc * sd_yb
gamma_median = gamma_median_sc * sd_yb
boxplot((cbind(gamma_median, gamma_qad, gamma_saurm)),
        las=1,main=paste("N =",Nb,"  S =",S),
        ylab="Fractionation parameter gamma")
abline(h=gamma,lty=1,col=4,lwd=2)

# Medians and highest posterior density (HPD) 95% credible intervals for MEM model parameters
# after back-transformation to original scales
# Scale parameter nu
nu_median = nu_median_sc * sd_xb
summary(nu_median)

# Expected value for standard deviation of truncated density averaged over
# entire study range, omega
omega = numeric(length(nu_median))
for(i in 1:length(nu_median)) {
  f_sd = function(x) ftgaussian_sig(x, nu_median[i], xl, xu) * dunif(x, xl, xu)
  omega[i] = hcubature(f_sd, xl, xu)$integral
}
summary(omega)
boxplot((cbind(omega, nu_qad, nu_saurm)),
        las=1,main=paste("N =",Nb,"  S =",S),ylab="Scale parameter omega")
abline(h=nu,lty=1,col=4,lwd=2)

# Regression slope parameter beta
beta_median = beta_median_sc * sd_yb / sd_xb
boxplot(beta_median,las=1,main=paste("N =",Nb,"  S =",S),
        ylab="Regression slope parameter beta")
abline(h=beta,lty=1,col=4,lwd=2)

# Regression intercept parameter alpha
alpha_median = mean_yb + alpha_median_sc * sd_yb -
                mean_xb * beta_median_sc * sd_yb / sd_xb
boxplot(alpha_median,las=1,main=paste("N =",Nb,"  S =",S),
        ylab="Regression intercept parameter alpha")
abline(h=alpha,lty=1,col=4,lwd=2)

# Regression SE parameter sigmab
sigmab_median = sd_yb * sigmab_median_sc
boxplot(sigmab_median,las=1,main=paste("N =",Nb,"  S =",S),
        ylab="Regression SE parameter sigmab")
abline(h=sigmab,lty=1,col=4,lwd=2)

# Regression SE parameter sigmac
sigmac_median = sd_yb * sigmac_median_sc
boxplot(sigmac_median,las=1,main=paste("N =",Nb,"  S =",S),
        ylab="Regression SE parameter sigmac")
abline(h=sigmac,lty=1,col=4,lwd=2)


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


################################################################################
# R libraries and functions sourced by R programs in "Measurement error
# models reveal the scale of consumer movements along an isoscape gradient"
#
# M.A. Rodriguez
################################################################################

################################################################################
# Function arguments:
# mu: location parameter
# s: scale parameter
# lb: lower bound of distribution
# ub: upper bound of distribution
# range: range of uniform distribution
################################################################################

# Load liibraries
library(rstan)
library(coda)
library(scoringRules)
library(extraDistr)
library(crch)
library(cubature)

# Set rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Define functions
# Function to obtain s.d. of truncated normal
ftgaussian_sig = function(mu, s, lb, ub){
  alp = (lb - mu) / s
  bet = (ub - mu) / s
  phi_alp = dnorm(alp, 0, 1)
  phi_bet = dnorm(bet, 0, 1)
  Z = pnorm(bet, 0, 1) - pnorm(alp, 0, 1)
  sqrt(s ^ 2 * (1 + (alp * phi_alp - bet * phi_bet) / Z - ((phi_alp - phi_bet) / Z) ^ 2))
}

# Obtain s.d. of truncated Laplace
fdtlaplace_sig = function(mu, s, lb, ub) {
  m1 = function(x) dtlaplace(x, mu, s, lb, ub) * x
  tlm = hcubature(m1, lb, mu)$integral + hcubature(m1, mu, ub)$integral
  m2 = function(x) dtlaplace(x, mu, s, lb, ub) * (x - tlm) ^ 2
  sqrt(hcubature(m2, lb, mu)$integral + hcubature(m2, mu, ub)$integral)
}

# Obtain s.d. of truncated Student t_3
fdtstudent_sig = function(mu, s, lb, ub) {
  m1 = function(x) dtt(x, mu, s, 3, lb, ub) * x
  tlm = hcubature(m1, lb, ub)$integral
  m2 = function(x) dtt(x, mu, s, 3, lb, ub) * (x - tlm) ^ 2
  sqrt(hcubature(m2, lb, ub)$integral)
}

# Obtain s.d. of uniform
funif_sig = function(mu, range, lb, ub){
  (min(mu + range / 2, ub) - max(mu - range / 2, lb)) / sqrt(12)
}

# Truncated Laplace density function
dtlaplace = function(x, mu, s, lb, ub) {
  alpha = plaplace(lb, mu, s)
  beta = plaplace(ub, mu, s)
  dlaplace(x, mu, s) / (beta - alpha)
}

# Generate random samples from truncated Laplace
rtlaplace = function(N, mu, b, lb, ub) {
  alpha = plaplace(lb, mu, b)
  beta = plaplace(ub, mu, b)
  p = runif(N, alpha, beta)
  mu - b * sign(p - 1/2) * log(1 - 2 * abs(p - 1/2))
}


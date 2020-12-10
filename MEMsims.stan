/* 
Stan code to generate results of the section "Results: Performance of
estimators in simulated scenarios" in "Measurement error models reveal the
scale of consumer movements along an isoscape gradient"

Is called by R program MEM_simulated_data.r

M.A. Rodriguez
*/

data {
  int<lower=1> Nb;
  int<lower=1> Nc;
  real xl_sc;
  real xu_sc;
  vector[Nb] yb_sc;
  vector[Nc] yc_sc;
  vector<lower=xl_sc, upper=xu_sc>[Nb] xb_sc;
  vector<lower=xl_sc, upper=xu_sc>[Nc] xc_sc;
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real<lower=0> sigmab;
  real<lower=0> sigmac;
  real<lower=0> nu;
  vector<lower=xl_sc, upper=xu_sc>[Nc] xstar_sc;
}
model {
  alpha ~ normal(0, 1e0);
  beta ~ normal(0, 1e1);
  gamma ~ normal(0, 1e1); 
  sigmab ~ student_t(3, 0, 1);
  sigmac ~ student_t(3, 0, 1);
  nu ~ student_t(3, 0, 1);
  xstar_sc ~ uniform(xl_sc, xu_sc);
  for(i in 1:Nc) {
    xc_sc[i] ~ normal(xstar_sc[i], nu) T[xl_sc, xu_sc];
  }
  yb_sc ~ normal(alpha +         beta * xb_sc,    sigmab);
  yc_sc ~ normal(alpha + gamma + beta * xstar_sc, sigmac);
}

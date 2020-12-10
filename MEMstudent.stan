/* 
Stan code to generate results for the section "Applications of the model: quantifying the movements of three fish species" in "Measurement error models reveal the scale of consumer movements along an isoscape gradient"

Is called by R program MEM_empirical_data.r

M.A. Rodriguez
*/

functions {
  // define randon number generator function for truncated Student-t
  real student_tt_rng (real mu, real sig, real lb, real ub) {
  real rt;
  int k = 0;
  while (k == 0) {
    rt = student_t_rng(3., mu, sig);
    if (rt < lb || rt > ub) continue;
      return rt;
    }
  }
}
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
  for(i in 1:Nc){
     xc_sc[i] ~ student_t(3, xstar_sc[i], nu) T[xl_sc, xu_sc];
  }
  yb_sc ~ normal(alpha +         beta * xb_sc,    sigmab);
  yc_sc ~ normal(alpha + gamma + beta * xstar_sc, sigmac);
}
generated quantities {
  vector[Nb] yb_sc_rep;
  vector[Nc] yc_sc_rep;
  vector[Nc] xc_sc_rep;
  vector[Nb + Nc] y_rep;

  // posterior predictive
  for(i in 1:Nb) yb_sc_rep[i] = normal_rng(alpha + beta * xb_sc[i], sigmab);
  for(j in 1:Nc) {
    yc_sc_rep[j] = normal_rng(alpha + gamma + beta * xstar_sc[j], sigmac);
    xc_sc_rep[j] = student_tt_rng(xstar_sc[j], nu, xl_sc, xu_sc);
    y_rep = append_row(yb_sc_rep, yc_sc_rep);
   }

}





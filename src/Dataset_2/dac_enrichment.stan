data {
  int<lower=0> N;
  int y [N];
  int n [N];
  int g [N];
}

parameters {
  real alpha;
  vector [max(g)] mu;
  vector <lower=0> [max(g)] sigma;
  vector [N] z;
}

transformed parameters {
    vector [N] beta;
    for(i in 1:N) {
        beta[i] = mu[g[i]] + sigma[g[i]] * z[i];
    }
}

model {
    alpha ~ normal(0, 5);
    mu ~ normal(0, 1);
    sigma ~ normal(0, 1);
    z ~ normal(0, 1);
    
    for(i in 1:N) {
        y[i] ~ binomial_logit(n[i], alpha + beta[i]);
    }
}

generated quantities {
    int yhat [N];
    real phat [N];
    real log_lik [N];
    real yhat_group [max(g)];
    real phat_group [max(g)];
    real d [6];
    real r [6];
    
    
    for(i in 1:N) {
        yhat[i] = binomial_rng(n[i], inv_logit(alpha + beta[i]));
        phat[i] = inv_logit(alpha + beta[i]);
        log_lik[i] = binomial_lpmf(y[i] | n[i], inv_logit(alpha + beta[i]));
    }
    
    for(i in 1:max(g)) {
        yhat_group[i] = inv_logit(normal_rng(alpha + mu[i], sigma[i]));
        phat_group[i] = inv_logit(alpha + mu[i]);
    }
    
    d[1] = phat_group[1]-phat_group[2];
    d[2] = phat_group[1]-phat_group[3];
    d[3] = phat_group[1]-phat_group[4];
    d[4] = phat_group[2]-phat_group[3];
    d[5] = phat_group[2]-phat_group[4];
    d[6] = phat_group[3]-phat_group[4];
    
    r[1] = phat_group[1]/phat_group[2];
    r[2] = phat_group[1]/phat_group[3];
    r[3] = phat_group[1]/phat_group[4];
    r[4] = phat_group[2]/phat_group[3];
    r[5] = phat_group[2]/phat_group[4];
    r[6] = phat_group[3]/phat_group[4];
}


data {
  // dimensions and slicing
  int<lower = 1> N;                 // observations
  int<lower = 2> T;                 // time points
  int<lower = N * T> NT;            // product
  int startstop[N, 2];              // slicer
  int<lower = 2> K;                 // N(states)

  // predictor dimensions
  int<lower = 0> Md_sha;            // N(shared parameters of discrete process)
  int<lower = 0> Md_var;            // N(varying parameters of discrete process)
  int<lower = 0> Mc_sha;            // N(varying parameters of continuous process)
  int<lower = 0> Mc_var;            // N(varying parameters of continuous process)

  // predictors
  matrix[NT, Md_sha] z_sha;         // matrix of shared parameters of discrete process
  matrix[NT, Md_var] z_var;         // matrix of shared parameters of discrete process
  matrix[NT, Mc_sha] x_sha;         // matrix of shared continuous of discrete process
  matrix[NT, Mc_var] x_var;         // matrix of shared continuous of discrete process

  // output
  vector[NT] y;                     // vector of output

  // priors
  // real prior_mean[Md];
}

parameters {
  // Discrete state model
  simplex[K] pi1[N];                 // initial state probabilities
  vector[Md_sha] gamma;              // tvtp coefficient
  vector[Md_var] lambda[K, K];          // tvtp coefficient

  // Continuous observation model
  ordered[K] mu;                     // observation means
  ordered[K] phi;                    // observation AR
  vector[Mc_sha] zeta;               // continuous observation model coefficients
  matrix[Mc_var, K] beta;            // continuous observation model coefficients
  real<lower=0> sigma;               // observation standard deviations
}

transformed parameters {
  vector[K] logalpha[T, N];
  simplex[K] A[N, T, K];             // A[t][i][j] = p(z_t = j | z_{t-1} = i)
  for (nn in 1:N){
    for (t in 1:T){
      for (i in 1:K){
        real den[K - 1];
        real den_sum;
        for (j in 1:K){
          if (i == j){
            den[j] = exp(z_var[startstop[nn, 1] + t - 1] * lambda[i, j]);
          } else {
            den[j] = 0;
          }
        }
        den_sum = log_sum_exp(den);
        for (j in 1:K){
          if (i == j){
            A[nn, t, i, j] = exp(log(1) - log(1 + den_sum));
          }
          A[nn, t, i, j] = exp(log(exp(z_var[startstop[nn, 1] + t - 1] * lambda[i, j])) - log(den_sum));
        }
      }
    }
  }
  { // Forward algorithm log p(z_t = j | x_{1:t})
  // logalpha[t, j] = Pr(z_t = j| x_{1:t})
  // logalpha[t-1, i] = Pr(z_t = i| x_{1:t-1})
  // normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma):  p(x_t|z_t = j) local evidence at t for j

  for (nn in 1:N){
    for (j in 1:K){
      logalpha[1, nn, j] = log(pi1[nn][j]) + normal_lpdf(y[startstop[nn, 1]] | mu[j], sigma);
    }
    for (t in 2:T){
      for (j in 1:K) {
        real accumulator[K];
        for (i in 1:K) {
          accumulator[i] = logalpha[t - 1, nn, i] +
                            log(A[nn, t, i, j]) +
                            normal_lpdf(y[startstop[nn, 1] + t - 1] | mu[j] +
                                        phi[j] * y[startstop[nn, 1] + t - 2] +
                                        x_sha[startstop[nn, 1] + t - 1] * zeta +
                                        x_var[startstop[nn, 1] + t - 1] * beta[:,j], sigma);
        }
        logalpha[t, nn, j] = log_sum_exp(accumulator);
      }
    }
  }
  }
}

model {
  // structural
  target += normal_lpdf(mu | 0, 4);
  target += normal_lpdf(phi | 0.5, 0.25);

  // shared
  target += normal_lpdf(gamma | 0, 0.5);
  target += normal_lpdf(zeta | 0, 4);

  // varying
  for (i in 1:Md_var){
    for (j in 1:K){
      target += normal_lpdf(lambda[:,j, i] | 0, 2);
    }
  }
  for (i in 1:K)
    target += normal_lpdf(beta[:,i] | 0, 5);
  }

  // likelihood
  for (n in 1:N){
    target += log_sum_exp(logalpha[T, n]);
  }
}

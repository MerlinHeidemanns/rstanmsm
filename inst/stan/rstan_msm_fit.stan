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
}

parameters {
  // Discrete state model
  simplex[K] pi1[N];                    // initial state probabilities
  vector[Md_sha] gamma;              // tvtp coefficient
  matrix[Md_var, K] lambda;              // tvtp coefficient

  // Continuous observation model
  ordered[K] mu;                     // observation means
  ordered[K] phi;                    // observation AR
  vector[Mc_sha] zeta;                // continuous observation model coefficients
  matrix[Mc_var, K] beta;                // continuous observation model coefficients
  real<lower=0> sigma;               // observation standard deviations
}

transformed parameters {
  vector[2] logalpha[T, N];
  simplex[2] A[N, T, 2];                // A[t][i][j] = p(z_t = j | z_{t-1} = i)
  for (nn in 1:N){
    for (t in 1:T){
      A[nn, t, 1, 1] = normal_cdf(z_sha[startstop[nn, 1] + t - 1] *   gamma +
                                  z_var[startstop[nn, 1] + t - 1] *  lambda[:,1] , 0, 1);
      A[nn, t, 2, 2] = normal_cdf(z_sha[startstop[nn, 1] + t - 1] *   gamma +
                                  z_var[startstop[nn, 1] + t - 1] *  lambda[:,2], 0, 1);
      A[nn, t, 1, 2] = 1 - A[nn, t, 1, 1];
      A[nn, t, 2, 1] = 1 - A[nn, t, 2, 2];
    }
  }
  { // Forward algorithm log p(z_t = j | x_{1:t})
  // logalpha[t, j] = Pr(z_t = j| x_{1:t})
  // logalpha[t-1, i] = Pr(z_t = i| x_{1:t-1})
  // normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma):  p(x_t|z_t = j) local evidence at t for j

  for (nn in 1:N){
    for (j in 1:2){
      logalpha[1, nn, j] = log(pi1[nn][j]) + normal_lpdf(y[startstop[nn, 1]] | mu[j], sigma);
    }
    for (t in 2:T){
      for (j in 1:2) {
        real accumulator[2];
        for (i in 1:2) {
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
  for (i in 1:2){
    target += normal_lpdf(lambda[:,i] | 0, 2);
    target += normal_lpdf(beta[:,i] | 0, 5);
  }

  // likelihood
  for (n in 1:N){
    target += log_sum_exp(logalpha[T, n]);
  }
}

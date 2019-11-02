data {
  // dimensions and slicing
  int<lower = 1> N;                 // observations
  int<lower = 2> T;                 // time points
  int<lower = N * T> NT;            // product
  int startstop[N, 2];              // slicer
  int<lower = 2> K;                 // N(states)

  // predictor dimensions
  int<lower = 0> Mx_sha;            // N(varying parameters of continuous process)
  int<lower = 0> Mx_var;            // N(varying parameters of continuous process)

  // predictors
  // matrix[NT, Md_sha] z_sha;         // matrix of i->-i transience parameters parameters
  matrix[NT, Mx_sha] x_sha;         // matrix of shared continuous of discrete process
  matrix[NT, Mx_var] x_var;         // matrix of varying continuous of discrete process

  // output
  vector[NT] y;                     // vector of output

  // intercepts
  int has_intercept;             // 1: common, 2: state-specific, 3: state-state specific
}

parameters {
  // Discrete state model
  simplex[K] pi1[N];                 // initial state probabilities
  simplex[K] lambda[K];              // tvtp coefficient (only initialization)

  // Continuous observation model
  ordered[K] mu;                     // observation means
  ordered[K] phi;                    // observation AR
  vector[Mx_sha] alpha;                 // continuous observation model coefficients
  matrix[Mx_var, K] beta;            // continuous observation model coefficients
  real<lower=0> sigma;               // observation standard deviations
}

transformed parameters {
  vector[K] logalpha[T, N];

  { // Forward algorithm log p(z_t = j | x_{1:t})
  // logalpha[t, j] = Pr(z_t = j| x_{1:t})
  // logalpha[t-1, i] = Pr(z_t = i| x_{1:t-1})
  // normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma):  p(x_t|z_t = j) local evidence at t for j

  for (nn in 1:N){
    for (j in 1:K){
      logalpha[1, nn, j] = log(pi1[nn, j]) + normal_lpdf(y[startstop[nn, 1]] | mu[j], sigma);
    }
    for (t in 2:T){
      for (j in 1:K) {
        real accumulator[K];
        for (i in 1:K) {
          accumulator[i] = logalpha[t - 1, nn, i] +
                            log(lambda[i, j]) +
                            normal_lpdf(y[startstop[nn, 1] + t - 1] | mu[j] +
                                        phi[j] * y[startstop[nn, 1] + t - 2] +
                                        x_sha[startstop[nn, 1] + t - 1] * alpha +
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
  target += normal_lpdf(mu | 0, 5);
  target += beta_lpdf(phi | 5, 5);

  // varying
  for (k in 1:K){
    target += dirichlet_lpdf(lambda[k] | rep_vector(1, K));
  }

  target += normal_lpdf(alpha | 0, 2);
  for (i in 1:K){
    target += normal_lpdf(beta[:,i] | 0, 2);
  }

  // likelihood
  for (n in 1:N){
    target += log_sum_exp(logalpha[T, n]);
  }
}

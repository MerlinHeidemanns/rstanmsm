data {
  // dimensions and slicing
  int<lower = 1> N;                 // observations
  int<lower = 2> T;                 // time points
  int<lower = N * T> NT;            // product
  int startstop[N, 2];              // slicer
  int<lower = 2> K;                 // N(states)

  // predictor dimensions
  int<lower = 0> Mx_d;            // N(varying parameters of continuous process)
  int<lower = 0> Mx_e;            // N(varying parameters of continuous process)

  // predictors
  // matrix[NT, Md_sha] z_sha;         // matrix of i->-i transience parameters parameters
  matrix[NT, Mx_d] x_d;         // matrix of shared continuous of discrete process
  matrix[NT, Mx_e] x_e;         // matrix of varying continuous of discrete process

  // output
  vector[NT] y;                     // vector of output

  // intercepts
  int has_intercept[5];             // 1: alpha, 2: beta
}

transformed data {
  // add intercepts
  matrix[NT, has_intercept[1] + Mx_d] x_d_ = has_intercept[1] ? append_col(rep_vector(1.0, NT), x_d): x_d;
  matrix[NT, has_intercept[2] + Mx_e] x_e_ = has_intercept[2] ? append_col(rep_vector(1.0, NT), x_e): x_e;
  int Mx_d_ = Mx_d + has_intercept[1];
  int Mx_e_ = Mx_e + has_intercept[2];
}

parameters {
  // Discrete state model
  simplex[K] pi1[N];                 // initial state probabilities
  simplex[K] A[K];              // tvtp coefficient (only initialization)

  // Continuous observation model
  ordered[K] mu;                     // observation means
  ordered[K] phi;                    // observation AR
  vector [Mx_d_] alpha;                 // continuous observation model coefficients
  matrix[Mx_e_, K] beta;            // continuous observation model coefficients
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
                            log(A[i, j]) +
                            normal_lpdf(y[startstop[nn, 1] + t - 1] | mu[j] +
                                        phi[j] * y[startstop[nn, 1] + t - 2] +
                                        x_d_[startstop[nn, 1] + t - 1] * alpha +
                                        x_e_[startstop[nn, 1] + t - 1] * beta[:,j], sigma);
        }
        logalpha[t, nn, j] = log_sum_exp(accumulator);
      }
    }
  }
  }
}

model {

  target += normal_lpdf(sigma | 0, 5);

  // structural
  target += normal_lpdf(mu | 0, 5);
  target += normal_lpdf(phi | 0.5, 0.5);

  // varying
  for (k in 1:K){
    target += dirichlet_lpdf(A[k] | rep_vector(1, K));
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

generated quantities {
  // quantities
  vector[K] prS[T, N];
  real log_lik[NT];
  real y_hat[NT];

  // log likelihood
  for (nn in 1:N){
    for (t in 1:T){
      log_lik[startstop[nn, 1] + t - 1] = log_sum_exp(logalpha[t, nn]);
    }
  }

  // state_pr
  for (n in 1:N){
    for (t in 1:T){
      prS[t, n] = softmax(logalpha[t, n]);
    }
  }

  // y_hat
  for (nn in 1:N){
    y_hat[startstop[nn, 1]] = log_sum_exp(prS[1, nn] .* to_vector(normal_rng(mu, sigma)));
    for (t in 2:T){
      y_hat[startstop[nn, 1] + t - 1] = log_sum_exp(prS[t, nn] .* to_vector(normal_rng(mu +
                                        phi * y[startstop[nn, 1] + t - 2] +
                                        x_d_[startstop[nn, 1] + t - 1] * alpha +
                                        to_vector(x_e_[startstop[nn, 1] + t - 1] * beta), sigma)));
    }
  }
}

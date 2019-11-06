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
  matrix[NT, Mx_d] x_d;         // matrix of shared continuous of discrete process
  matrix[NT, Mx_e] x_e;         // matrix of varying continuous of discrete process

  // output
  vector[NT] y;                     // vector of output

  // intercepts
  int has_intercept[5];             // 1: alpha, 2: beta

  // shared
  int shared_TP;                       // 0: individual, 1: shared

  // ordering
  int order_x_e[Mx_e + has_intercept[2]]; // 0: unordered, 1: ordered

}

transformed data {
  // add intercepts
  int Mx_d_ = has_intercept[1] ? Mx_d + has_intercept[1] : (Mx_d == 0 ? 1 : Mx_d);
  int Mx_e_ = has_intercept[2] ? Mx_e + has_intercept[2] : (Mx_e == 0 ? 1 : Mx_e);
  matrix[NT, Mx_d_] x_d_ = has_intercept[1] ? append_col(rep_vector(1.0, NT), x_d):
                                             (Mx_d == 0       ? rep_matrix(0.0, NT, 1) : x_d);
  matrix[NT, Mx_e_] x_e_ = has_intercept[2] ? append_col(rep_vector(1.0, NT), x_e):
                                             (Mx_e == 0       ? rep_matrix(0.0, NT, 1) : x_e);

  // shared_TP
  int N_ = shared_TP ? 1 : N;

  // ordereing on x_e
  int Mx_e_un = Mx_e_ - sum(order_x_e);
  int Mx_e_ord = sum(order_x_e);
}

parameters {
  // Discrete state model
  simplex[K] pi1[N];                 // initial state probabilities
  simplex[K] A[N_, K];              // tvtp coefficient (only initialization)

  // Continuous observation model
  ordered[K] phi;                    // observation AR
  vector[Mx_d_] alpha;                 // continuous observation model coefficients
    // ordering
  vector[K] beta_un[Mx_e_un];
  ordered[K] beta_ord[Mx_e_ord];
  real<lower=0> sigma;               // observation standard deviations
}

transformed parameters {
  matrix[Mx_e_, K] beta;
  vector[K] logalpha[T, N];
  // create beta matrix from components
  {
    int count_ordered = 1;
    int count_unordered = 1;
    for (m in 1:Mx_e_){
      if (order_x_e[m]) {
        beta[m] = beta_ord[count_ordered]'; count_ordered += 1;
      } else {
        beta[m] = beta_un[count_unordered]'; count_unordered += 1;
      }
    }
  }
  { // Forward algorithm log p(z_t = j | x_{1:t})
  // logalpha[t, j] = Pr(z_t = j| x_{1:t})
  // logalpha[t-1, i] = Pr(z_t = i| x_{1:t-1})
  // normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma):  p(x_t|z_t = j) local evidence at t for j

    for (nn in 1:N){
      int nn_ = shared_TP ? 1 : nn;
      for (j in 1:K){
        logalpha[1, nn, j] = log(pi1[nn, j]) + normal_lpdf(y[startstop[nn, 1]] |
                                          x_d_[startstop[nn, 1]] * alpha +
                                          x_e_[startstop[nn, 1]] * beta[:,j], sigma);
      }
      for (t in 2:T){
        for (j in 1:K) {
          real accumulator[K];
          for (i in 1:K) {
            accumulator[i] = logalpha[t - 1, nn, i] +
                              log(A[nn_ , i, j]) +
                              normal_lpdf(y[startstop[nn, 1] + t - 1] |
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
  target += normal_lpdf(phi | 0.5, 0.5);

  // varying
  for (n in 1:N_){
    for (k in 1:K){
      target += dirichlet_lpdf(A[n,k] | rep_vector(1, K));
    }
  }
  target += normal_lpdf(alpha | 0, 2);
  for (i in 1:K){
    target += normal_lpdf(beta_un[:,i] | 0, 2);
    target += normal_lpdf(beta_ord[:,i] | 0, 2);
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
    y_hat[startstop[nn, 1]] = sum(prS[1, nn] .* to_vector(normal_rng(
                                        x_d_[startstop[nn, 1]] * alpha +
                                        to_vector(x_e_[startstop[nn, 1]] * beta), sigma)));
    for (t in 2:T){
      y_hat[startstop[nn, 1] + t - 1] = sum(prS[t, nn] .* to_vector(normal_rng(
                                        phi * y[startstop[nn, 1] + t - 2] +
                                        x_d_[startstop[nn, 1] + t - 1, :] * alpha +
                                        to_vector(x_e_[startstop[nn, 1] + t - 1] * beta)
                                        , sigma)));
    }
    // +
    //                                    to_vector(x_d_[startstop[nn, 1] + t - 1] * alpha) +
    //
  }
}

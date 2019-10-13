data {
  int<lower = 1> N;
  int<lower = 1> T;                   // number of observations (length)
  real NT[N, 2];                         // coverage
  vector[N] y[T];                        // observations
  int z[T];
}

parameters {
  // Discrete state model
  simplex[2] pi1;                   // initial state probabilities
  real gamma[2];                    // tvtp intercept
  real lambda[2];                      // tvtp coefficient

  // Continuous observation model
  ordered[2] mu;                    // observation means
  ordered[2] phi;                   // AR1 parameter
  real<lower=0> sigma;           // observation standard deviations
}

transformed parameters {
  vector[2] logalpha[T];             //
  simplex[2] A[T, 2];                // A[t][i][j] = p(z_t = j | z_{t-1} = i)

  for (t in 1:T){
    A[t, 1, 1] = normal_cdf(gamma[1] + lambda[1] * z[t], 0, 1);
    A[t, 2, 2] = normal_cdf(gamma[2] + lambda[2] * z[t], 0, 1);
    A[t, 1, 2] = 1 - A[t, 1, 1];
    A[t, 2, 1] = 1 - A[t, 2, 2];
  }
  { // Forward algorithm log p(z_t = j | x_{1:t})
  // logalpha[t, j] = Pr(z_t = j| x_{1:t})
  // logalpha[t-1, i] = Pr(z_t = i| x_{1:t-1})
  // normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma):  p(x_t|z_t = j) local evidence at t for j

    real accumulator[2];
    for (j in 1:2){
      logalpha[1, j] = log(pi1[j]) + normal_lpdf(y[1] | mu[j], sigma);
    }

    for (t in 2:T) {
      for (j in 1:2) { // j = current (t)
        for (i in 1:2) { // i = previous (t-1)
                         // Murphy (2012) Eq. 17.48
                         // belief state      + transition prob + local evidence at t
                         // logalpha is the vector of probabilities of being in state j at t -1
                         // log(A[i, j]) is the transition probability
                         // normal_lpdf is the p of the local evidence given s_t = j
                         // then this is mixed in log_sum_exp
          accumulator[i] = logalpha[t-1, i] + log(A[t, i, j]) + normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma);
        } // run over the K states to assess the evidence for j
        logalpha[t, j] = log_sum_exp(accumulator); // then put this through log_sum_exp and put p(s_t = j|x_{1:t}) into logalpha[t, j];
        // HENCE logalpha[t-1, i] + log(A[t, i, j]) == log(pti) and we can introduce multiple observations through normal_lpdf(y[t] | mu[j] + phi[j] * y[t - 1], sigma);
        // shouldn't this just work? y[t] is a vector and mu[j]... become a real
        // y[t - 1] is a real[ ]
      }
    }
  } // Forward
}

model {
  target += normal_lpdf(mu | 4, 2);
  target += normal_lpdf(phi | 0.5, 0.25);
  target += normal_lpdf(gamma | 0.5, 0.25);
  target += normal_lpdf(lambda | 0, 2);
  target += log_sum_exp(logalpha[T]); // Note: update based only on last logalpha
}

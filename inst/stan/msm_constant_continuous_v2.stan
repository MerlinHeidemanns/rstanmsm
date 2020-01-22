data {
  // dimensions and slicing
  int<lower = 1> N;                 // units
  int<lower = 2> T;                 // maximum of time points
  int<lower = N> NT;                // number of observations

  // state process
  int<lower = 1, upper = N> NS;     // N of state processes
  int<lower = 1, upper = NS> NTP;   // N of transition probabilities
  int id_tp[NS];                    // which state process belongs to which transition probability matrix
  int<lower = 0> slicer_T[NS, 2];   // min(N, T):max(N, T)
  int startstop[T, NS * 2];         // slicer
  int<lower = 2> K;                 // N(states)

  // predictor dimensions
  int has_intercept[5];             // 1: alpha, 2: beta; 3: gamma; 4: delta; 5: eta
  int<lower = 0> Mz;
  int<lower = 0> Mx_d;              // N(varying parameters of continuous process)
  int<lower = 0> Mx_e;              // N(varying parameters of continuous process)
  int pp1; int pp2; int pp3;        // N of general, state, and state-state specific predictors
  int<lower = 0, upper = 1> pp_lambda[3, Mz];
  int<lower = 0, upper = 1> pp_gamma[pp2];
  int<lower = 0, upper = 1> pp_eta[pp2];

  // predictors
  matrix[NTP * T, Mz] z;         // matrix of predictors of discrete process
  matrix[NT, Mx_d] x_d;             // matrix of shared predictors of continuous process
  matrix[NT, Mx_e] x_e;             // matrix of varying predictors of continuous process

  // output
  vector[NT] y;                     // vector of output

  // shared
  int state_sigma;                  // 0: general,    1: state-specific

  // time varying transition probabilities
  int tvtp;                         // 0: No,         1: Yes

  // ordering
  int order_x_e[Mx_e + has_intercept[2]]; // 0: unordered, 1: ordered

  // priors
    // dirichlet
  vector<lower = 0>[K] A_prior;
  real priors[7, 4];     // 1: Kind, 2: mean, 3: sd, 4: df
                         // 1: normal, 2: cauchy, 3: student-t

  // missing
  real id_miss[NS * T];  // 0: missing, 1: present / at least one observation
}

transformed data {
  // adjust N of predictors
  int Mz_ = Mz + sum(has_intercept[3:5]); // add one if there is an intercept
  int Mx_d_ = has_intercept[1] ? Mx_d + has_intercept[1] : (Mx_d == 0 ? 1 : Mx_d); // 1st ifelse add 1, 2nd ifelse if empty then one
  int Mx_e_ = has_intercept[2] ? Mx_e + has_intercept[2] : (Mx_e == 0 ? 1 : Mx_e); // 1st ifelse add 1, 2nd ifelse if empty then one

  // adjust matrixes
    // 1st ifelse add intercept, 2nd ifelse if empty then matrix of 0s
  matrix[NTP * T, Mz_] z_ = sum(has_intercept[3:5]) ? append_col(rep_vector(1.0, NS * T), z): z;
  matrix[NT, Mx_d_] x_d_ = has_intercept[1] ? append_col(rep_vector(1.0, NT), x_d): (Mx_d == 0 ? rep_matrix(0.0, NT, 1) : x_d);
  matrix[NT, Mx_e_] x_e_ = has_intercept[2] ? append_col(rep_vector(1.0, NT), x_e): (Mx_e == 0 ? rep_matrix(0.0, NT, 1) : x_e);

  // adjust ordering on discrete process parameters
  int pp1_ = pp1 + has_intercept[3]; // increment by one
  int pp2_ = pp2 + has_intercept[4];
  int pp3_ = pp3 + has_intercept[5];
  matrix[3, Mz_] pp_lambda_ = sum(has_intercept[3:5]) ?
            append_col(to_vector(has_intercept[3:5]), to_matrix(pp_lambda))
            : to_matrix(pp_lambda);
  vector[pp2_] pp_gamma_ = sum(has_intercept[3:5]) ? append_row(has_intercept[4], to_vector(pp_gamma)) : to_vector(pp_gamma);

  // ordering on x_e
  int Mx_e_un = Mx_e_ - sum(order_x_e);
  int Mx_e_ord = sum(order_x_e);

  // specificity of sigma
  int K_sigma = state_sigma ? K : 1;

  // time varying transition probabilities
  int T_ = tvtp ? T : 1;
}

parameters {
  // Discrete state model
  simplex[K] pi1[NS];                       // initial state probabilities
  simplex[K] A[tvtp ? T : 1, NTP, K];       // tvtp coefficient (only initialization)
  real delta[tvtp * Mz_];                   // general
  real gamma[tvtp * pp2_, K * tvtp];        // state
  real eta[tvtp * K, pp3_, K - 1];          // state - state

  // Continuous observation model
  vector[Mx_d_] alpha;                      // continuous observation model coefficients

  // ordering
  vector[K] beta_un[Mx_e_un];
  ordered[K] beta_ord[Mx_e_ord];

  // variance term
  real<lower=0> sigma[K_sigma];             // observation standard deviations
}

transformed parameters {
  matrix[K, K] logA[T_, NTP];
  matrix[Mx_e_, K] beta;
  vector[K] logalpha[NS * T];
  vector[K] llh_tmp[NS * T];
  matrix[Mz_, K] lambda_prime[K];
  // lambda_prime
  {
    for (i in 1:K){
      int count3 = 1;
      for (j in 1:K){
        int count2 = 1;
        int count1 = 1;
        for (q in 1:Mz_){
          if (i == j){
            lambda_prime[i, q, j] = 0.0;
          } else {
            if (pp_lambda_[1, q]) {                            // same persistence/transience
              lambda_prime[i, q, j] = delta[q];
            } else if (pp_lambda_[2, q]) {                     // state specific persistence/transience
              while (pp_gamma_[count2] != 1) {
                count2 += 1;
              }
              lambda_prime[i, q, j] = gamma[count2, i];
              count2 += 1;
            } else if (pp_lambda_[3, q]) {                     // state-state specific persistence/transience
              lambda_prime[i, q, j] = eta[i, count1, count3];
              count1 += 1;
            }
            if (q == Mz_) {
              count3 += 1;
            }
          }
        }

      }
    }
  }
  // fill logalpha and temporary loglikelihood
  for (j in 1:K){
    logalpha[:, j] = rep_array(0.0, NS * T);
    llh_tmp[:, j]  = rep_array(0.0, NS * T);
    if (tvtp == 0){
      for (i in 1:K) for (n in 1:NTP) logA[1, n, i, j] = log(A[1, n, i, j]);
    } else {
      for (n in 1:NTP) for (t in 1:T) logA[t, n, j] = softmax(to_vector(z_[t + T * (n - 1)] * lambda_prime[j]))';
    }
  }

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

  // temporary loglikelihood of p(y | S_t = j)
  for (ns in 1:NS) {
    {
      int n_tp = id_tp[ns];
      int t_ns0 = slicer_T[ns, 1]; int t_ns1 = slicer_T[ns, 2]; // relevant rows in the slicer
      for (t in t_ns0:t_ns1) {
        if (id_miss[t + (ns - 1) * T]){
          llh_tmp[t + (ns - 1) * T] = rep_vector(0.0, K);
        } else {
          {
            int t0 = startstop[t, (ns * 2) - 1]; int t1 = startstop[t, (ns * 2)];
            for (j in 1:K){
              {
                int j_ = state_sigma ? j : 1;
                llh_tmp[t + (ns - 1) * T, j] = normal_lpdf(y[t0:t1] | x_d_[t0:t1] * alpha + x_e_[t0:t1] * beta[:,j], sigma[j_]);
              }
            }
          }
        }
      }
    // forward algorithm
      logalpha[t_ns0 + (ns - 1) * T] = log(pi1[ns]) + llh_tmp[ns, t_ns0];
      for (t in (t_ns0 + 1):t_ns1) {
        for (j in 1:K){
          logalpha[t + (ns - 1) * T, j] = log_sum_exp(logalpha[t + (ns - 1) * T - 1] + logA[1, n_tp, :,j] + rep_vector(llh_tmp[t + (ns - 1) * T, j], K));
        }
      }
    }
  }
}

model {
  // transition probabilities
  if (tvtp == 0){
    for (n in 1:NTP){
      for (k in 1:K){
        target += dirichlet_lpdf(A[1, n, k] | A_prior);
      }
    }
  } else {
    {
      if (priors[1, 1] == 1){
        target += normal_lpdf(delta | priors[1, 2], priors[1, 3]);
      } else if (priors[1, 1] == 2){
        target += cauchy_lpdf(delta | priors[1, 2], priors[1, 3]);
      } else if (priors[1, 1] == 3){
        target += student_t_lpdf(delta | priors[1, 4], priors[1, 2], priors[1, 3]);
      }
      // gamma
      if (priors[2, 1] == 1){
        for (i in 1:pp2_) target += normal_lpdf(gamma[i] | priors[2, 2], priors[2, 3]);
      } else if (priors[2, 1] == 2){
        for (i in 1:pp2_) target += cauchy_lpdf(gamma[i] | priors[2, 2], priors[2, 3]);
      } else if (priors[2, 1] == 3){
        for (i in 1:pp2_) target += student_t_lpdf(gamma[i] | priors[2, 4], priors[2, 2], priors[2, 3]);
      }
      // eta
      if (priors[3, 1] == 1){
        for (i in 1:pp3_) for (j in 1:K) target += normal_lpdf(eta[j, i,:] | priors[3, 2], priors[3, 3]);
      } else if (priors[3, 1] == 2){
        for (i in 1:pp3_) for (j in 1:K) target += cauchy_lpdf(eta[j, i,:] | priors[3, 2], priors[3, 3]);
      } else if (priors[3, 1] == 3){
        for (i in 1:pp3_) for (j in 1:K) target += student_t_lpdf(eta[j, i,:] | priors[3, 4], priors[3, 2], priors[3, 3]);
      }
    }
  }

  // alpha
  if (priors[4, 1] == 1){
    target += normal_lpdf(alpha | priors[4, 2], priors[4, 3]);
  } else if (priors[4, 1] == 2){
    target += cauchy_lpdf(alpha | priors[4, 2], priors[4, 3]);
  } else if (priors[4, 1] == 3){
    target += student_t_lpdf(alpha | priors[4, 4], priors[4, 2], priors[4, 3]);
  }

  // beta unordered
  if (priors[5, 1] == 1){
    for (i in 1:K) target += normal_lpdf(beta_un[:,i] | priors[5, 2], priors[5, 3]);
  } else if (priors[5, 1] == 2){
    for (i in 1:K) target += cauchy_lpdf(beta_un[:,i] | priors[5, 2], priors[5, 3]);
  } else if (priors[5, 1] == 3){
    for (i in 1:K) target += student_t_lpdf(beta_un[:,i] | priors[5, 4], priors[5, 2], priors[5, 3]);
  }

  // beta ordered
  if (priors[6, 1] == 1){
    for (i in 1:K) target += normal_lpdf(beta_ord[:,i] | priors[6, 2], priors[6, 3]);
  } else if (priors[6, 1] == 2){
    for (i in 1:K) target += cauchy_lpdf(beta_ord[:,i] | priors[6, 2], priors[6, 3]);
  } else if (priors[6, 1] == 3){
    for (i in 1:K) target += student_t_lpdf(beta_ord[:,i] | priors[6, 4], priors[6, 2], priors[6, 3]);
  }

  // sigma
  if (priors[7, 1] == 1){
    target += normal_lpdf(sigma | priors[7, 2], priors[7, 3]);
  } else if (priors[7, 1] == 2){
    target += cauchy_lpdf(sigma | priors[7, 2], priors[7, 3]);
  } else if (priors[7, 1] == 3){
    target += student_t_lpdf(sigma | priors[7, 4], priors[7, 2], priors[7, 3]);
  }


  // likelihood
  for (ns in 1:NS){
    {
      int t_ns1 = slicer_T[ns, 2];
      target += log_sum_exp(logalpha[t_ns1 + (ns - 1) * T]);
    }
  }
}

generated quantities {
  vector[K] prS[NS * T];
  vector[NT] y_rep;                     // vector of output

  for (j in 1:K){
    prS[:, j] = rep_array(0.0, NS * T);
  }
  for (n in 1:(NS * T)){
    prS[n] = softmax(logalpha[n]);
  }
  {
    vector[K] sigma_ = state_sigma ? to_vector(sigma) : rep_vector(sigma[1], K);
    for (ns in 1:NS) {
      {
        int t_ns0 = slicer_T[ns, 1]; int t_ns1 = slicer_T[ns, 2]; // relevant rows in the slicer
        for (t in t_ns0:t_ns1) {
          {
            int t0 = startstop[t, (ns * 2) - 1]; int t1 = startstop[t, (ns * 2)];
            int count = categorical_rng(prS[t + (ns - 1) * T]);
            for (t2 in t0:t1){
              y_rep[t2] = normal_rng(x_d_[t2] * alpha + x_e_[t2] * beta[:,count], sigma_[count]);
            }
          }
        }
      }
    }
  }
}









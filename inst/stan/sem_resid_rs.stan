functions {
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
  real generalized_double_pareto_lpdf(vector x, real alpha) {
    // generalized double Pareto
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/
    return(sum(
      -(alpha + 1.0) * log(1.0 + fabs(x) / alpha)
    ));
  }
  real eff(int p, real x) {
    return(
      2 * lmgamma(p, x / 2) - x * p * log(x / 2) + x * p
    );
  }
  real ln_det_spd(matrix S) {
    return(2 * sum(log(diagonal(cholesky_decompose(S)))));
  }
  real gen_matrix_beta_ii_lpdf(matrix S, matrix Omega, real n, real m, real ln_det_S) {
    int p = rows(S);
    real F_1 = eff(p, m) + eff(p, n) - eff(p, m + n);
    real F_2 = -((n - p - 1) * ln_det_S) - (m * ln_det_spd(Omega)) +
      ((m + n) * ln_det_spd((m * Omega + n * S) / (m + n)));
    real ll = (F_1 + F_2) / -2.0;
    return(ll);
  }
}
data {
  int Np;
  int Ni;
  matrix[Ni, Ni] S;
  int Nf;
  int Nce;
  array[Nce, 2] int error_mat; // cor error matrix
  matrix[Ni, Nf] loading_pattern;
  array[Nf] int markers; // markers
  int Nf_corr;
  array[Nf_corr, 2] int F_corr_mat;
  array[Nf, Nf] int coef_pattern;
  real<lower = 0> sc_par;  // sigma_coefficients parameter
  real<lower = 0> sl_par;  // sigma_loading parameter
  real<lower = 0> rs_par;  // residual sd parameter
  real<lower = 1> rc_par;  // residual corr parameter
  real<lower = 1> fc_par; // beta prior shape for phi
  int<lower = 1, upper = 100> method; // which method
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl = 0; // N_non-zero loadings
  // cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  matrix[Ni, Ni] N_S = (Np - 1) * S;  // covariance matrix
  int coef_count = 0;
  int outcome_count = 0;
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  real ln_det_S = ln_det_spd(S);
  int N_Sigma = 1;

  if (method >= 90) {
    Nisqd2 = 0;
  }

  if (method == 91) {
    N_Sigma = Ni;
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  for (i in 1:Nf) {
    if (sum(coef_pattern[i, ]) > 0) outcome_count += 1;
    for (j in 1:Nf) {
      if (coef_pattern[i, j] == 1) coef_count += 1;
    }
  }
  for (i in 1:Ni) {
    for (j in 1:Nf) {
      if (loading_pattern[i, j] == 1) Nl += 1;
    }
  }
}
parameters {
  vector<lower = 0.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;
  vector<lower = 0>[Nl - Nf] loadings; // loadings
  real<lower = 0> sigma_loadings; // sd of loadings, hyperparm
  vector<lower = 0>[Nf] phi_sd;
  vector<lower = 0>[Ni] res_sds; // item residual sds heteroskedastic
  vector<lower = 0, upper = 1>[Nf_corr] phi_cor_01;
  vector<lower = 0, upper = 1>[Nce] res_cor_01; // correlated errors on 01
  vector[coef_count] coefs;
  vector<lower = 0>[outcome_count] sigma_coefs;
  cov_matrix[N_Sigma] Sigma;
}
model {
  rms_src_p ~ std_normal();
  if (method == 1) {
    // normal
    resids ~ std_normal();
  } else if (method == 2) {
    // lasso
    resids ~ double_exponential(0, 1);
  } else if (method == 3) {
    // logistic
    resids ~ logistic(0, 1);
  } else if (method == 4) {
    // https://www.mdpi.com/2075-1680/11/9/462
    gdp_alpha ~ lognormal(1, 1);
    target += generalized_double_pareto_lpdf(
      resids | gdp_alpha[1]);
  }

  loadings ~ normal(0, sigma_loadings);
  sigma_loadings ~ student_t(3, 0, sl_par);
  // TODO: #8 set prior as argument here or use standardized approach
  phi_sd ~ student_t(3, 0, 1);

  res_sds ~ student_t(3, 0, rs_par);

  phi_cor_01 ~ beta(fc_par, fc_par);
  res_cor_01 ~ beta(rc_par, rc_par);

  coefs ~ std_normal();
  sigma_coefs ~ student_t(3, 0, sc_par);

  {
    real m;
    matrix[Ni, Ni] Omega;
    vector[Ni] res_var = square(res_sds);
    vector[Nce] res_cor = res_cor_01 * 2 - 1;
    matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
    matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
    vector[Nf] phi_var = square(phi_sd);
    vector[Nf_corr] phi_cor = phi_cor_01 * 2 - 1;
    matrix[Nf, Nf_corr] F_corr_pe = rep_matrix(0, Nf, Nf_corr);
    matrix[Nf, Nf] F_cov_mat;
    vector[Nf] F_var_resid;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    matrix[Ni, Nf] Lambda_One_min_Beta_inv;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;

    {
      array[3] int pos_3 = rep_array(0, 3);
      for (i in 1:Nf) {
        if (sum(coef_pattern[i, ]) > 0) pos_3[1] += 1;
        for (j in 1:Nf) {
          if (coef_pattern[i, j] == 1) {
            pos_3[2] += 1;
            Coef_mat[i, j] = sigma_coefs[pos_3[1]] * coefs[pos_3[2]];
          }
        }
      }
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            if (i == markers[j]) Load_mat[i, j] = 1;
            else {
              pos_3[3] += 1;
              Load_mat[i, j] = loadings[pos_3[3]];
            }
          }
        }
      }
    }

    Lambda_One_min_Beta_inv = Load_mat * inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);

    for (i in 1:Nf_corr) {
      F_corr_pe[F_corr_mat[i, 1], i] = sqrt(
        abs(phi_cor[i]) * phi_var[F_corr_mat[i, 1]]);
      F_corr_pe[F_corr_mat[i, 2], i] = sign(phi_cor[i]) * sqrt(
        abs(phi_cor[i]) * phi_var[F_corr_mat[i, 2]]);
    }

    F_cov_mat = tcrossprod(F_corr_pe);
    F_var_resid = phi_var - diagonal(F_cov_mat);

    for (i in 1:Nce) {
      loading_par_exp[error_mat[i, 1], i] = sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 1]]);
      loading_par_exp[error_mat[i, 2], i] = sign(res_cor[i]) * sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 2]]);
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    Omega = add_diag(
      quad_form(
        add_diag(F_cov_mat, F_var_resid),
        Lambda_One_min_Beta_inv') + loading_par_exp_2,
      delta_mat_ast);

    total_var = diagonal(Omega);

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          pos += 1;
          Omega[i, j] += resids[pos] * rms_src_p[1] * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    if (method != 91) {
      Sigma ~ inv_wishart(1000, identity_matrix(1));
    }

    if (method >= 90 && method <= 99) {
      m = 1.0 / square(rms_src_p[1]) + Ni - 1;
      if (method == 90) {
        target += gen_matrix_beta_ii_lpdf(S | Omega, Np - 1.0, m, ln_det_S);
      } else if (method == 91) {
        Sigma ~ inv_wishart(m, m * Omega);
        target += wishart_lpdf(N_S | Np - 1, Sigma);
      }
    } else {
      target += wishart_lpdf(N_S | Np - 1, Omega);
    }
  }
}
generated quantities {
  real D_obs;
  real D_rep;
  real<lower = 0, upper = 1> ppp;
  real<lower = 0> rms_src;  // RMSE of residuals
  matrix[Ni, Nf] Load_mat_u = rep_matrix(0, Ni, Nf);
  matrix[Nf, Nf] Coef_mat_u = rep_matrix(0, Nf, Nf);
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
  vector[Nf] phi_var = square(phi_sd);
  vector[Nf] r_square;
  vector[Nf_corr] phi_cor = phi_cor_01 * 2 - 1;
  vector[Ni] res_var = square(res_sds);
  vector[Nce] res_cor = res_cor_01 * 2 - 1;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);

  if (method != 100) {
    rms_src = rms_src_p[1];
  } else {
    rms_src = 0.0;
  }
  if (method == 2) {
    rms_src *= sqrt_two;
  } else if (method == 3) {
    rms_src *= pi_sqrt_three;
  } else if (method == 4) {
    rms_src *= sqrt_two * gdp_alpha[1] / sqrt(
      (gdp_alpha[1] - 1.0) * (gdp_alpha[1] - 2.0)
    );
  }

  if (method < 90) {
    int pos = 0;
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        pos += 1;
        Resid[i, j] = resids[pos] * rms_src_p[1];
        Resid[j, i] = Resid[i, j];
      }
    }
  }

  {
    real m;
    matrix[Ni, Ni] Omega;
    matrix[Ni, Ni] Sigma_p;
    matrix[Ni, Ni] S_sim;
    matrix[Nf, Nf_corr] F_corr_pe = rep_matrix(0, Nf, Nf_corr);
    matrix[Nf, Nf] F_cov_mat;
    vector[Nf] F_var_resid;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    matrix[Nf, Nf] One_min_Beta_inv;
    matrix[Ni, Nf] Lambda_One_min_Beta_inv;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    vector[Nf] d_rt_c_hat;

    {
      array[3] int pos_3 = rep_array(0, 3);
      for (i in 1:Nf) {
        if (sum(coef_pattern[i, ]) > 0) pos_3[1] += 1;
        for (j in 1:Nf) {
          if (coef_pattern[i, j] == 1) {
            pos_3[2] += 1;
            Coef_mat[i, j] = sigma_coefs[pos_3[1]] * coefs[pos_3[2]];
          }
        }
      }
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            if (i == markers[j]) Load_mat[i, j] = 1;
            else {
              pos_3[3] += 1;
              Load_mat[i, j] = loadings[pos_3[3]];
            }
          }
        }
      }
    }
    Coef_mat_u = Coef_mat;
    Load_mat_u = Load_mat;

    One_min_Beta_inv = inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);
    Lambda_One_min_Beta_inv = Load_mat * One_min_Beta_inv;

    for (i in 1:Nf_corr) {
      F_corr_pe[F_corr_mat[i, 1], i] = sqrt(
        abs(phi_cor[i]) * phi_var[F_corr_mat[i, 1]]);
      F_corr_pe[F_corr_mat[i, 2], i] = sign(phi_cor[i]) * sqrt(
        abs(phi_cor[i]) * phi_var[F_corr_mat[i, 2]]);
    }

    F_cov_mat = tcrossprod(F_corr_pe);
    F_var_resid = phi_var - diagonal(F_cov_mat);

    for (i in 1:Nce) {
      loading_par_exp[error_mat[i, 1], i] = sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 1]]);
      loading_par_exp[error_mat[i, 2], i] = sign(res_cor[i]) * sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 2]]);
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    Omega = add_diag(
      quad_form(
        add_diag(F_cov_mat, F_var_resid),
        Lambda_One_min_Beta_inv') + loading_par_exp_2,
      delta_mat_ast);

    total_var = diagonal(Omega);

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          pos += 1;
          Omega[i, j] += resids[pos] * rms_src_p[1] * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    if (method >= 90 && method <= 99) {
      if (method == 90) {
        m = 1.0 / square(rms_src_p[1]) + Ni - 1;
        Sigma_p = inv_wishart_rng(m, m * Omega);
        S_sim = wishart_rng(Np - 1.0, Sigma_p / (Np - 1.0));
        D_obs = -2.0 * gen_matrix_beta_ii_lpdf(S | Omega, Np - 1.0, m, ln_det_S);
        D_rep = -2.0 * gen_matrix_beta_ii_lpdf(S_sim | Omega, Np - 1.0, m, ln_det_S);
      } else if (method == 91) {
        S_sim = wishart_rng(Np - 1.0, Sigma / (Np - 1.0));
        D_obs = -2.0 * wishart_lpdf(S | Np - 1.0, Sigma);
        D_rep = -2.0 * wishart_lpdf(S_sim | Np - 1.0, Sigma);
      }
    } else {
      S_sim = wishart_rng(Np - 1.0, Omega / (Np - 1.0));
      D_obs = -2.0 * wishart_lpdf(S | Np - 1.0, Omega / (Np - 1.0));
      D_rep = -2.0 * wishart_lpdf(S_sim | Np - 1.0, Omega / (Np - 1.0));
    }
    ppp = D_rep > D_obs ? 1.0 : 0.0;

    // Bollen (1989), pg 350
    d_rt_c_hat = sqrt(diagonal(quad_form(
      add_diag(F_cov_mat, F_var_resid),
      One_min_Beta_inv'
    )));

    for (j in 1:Nf) {
      Load_mat[, j] *= d_rt_c_hat[j];
      Coef_mat[, j] *= d_rt_c_hat[j];
      Coef_mat[j, ] /= d_rt_c_hat[j];
      r_square[j] = 1.0 - phi_var[j] / square(d_rt_c_hat[j]);
    }
    for (i in 1:Ni) {
      Load_mat[i, ] /= sqrt(total_var[i]);
    }
  }

  for (i in 1:Nce) {
    res_cov[i] = res_cor[i] * prod(res_sds[error_mat[i, ]]);
  }
}

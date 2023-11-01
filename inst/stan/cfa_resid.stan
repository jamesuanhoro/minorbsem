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
      -(alpha + 1.0) * log(1.0 + abs(x) / alpha)
    ));
  }
  real eff(int p, real x) {
    return(
      2 * lmgamma(p, x / 2) - x * p * log(x / 2) + x * p
    );
  }
  real gen_matrix_beta_ii_lpdf(matrix S, matrix Omega, real n, real m, real ln_det_S) {
    int p = rows(S);
    real F_1 = eff(p, m) + eff(p, n) - eff(p, m + n);
    real F_2 = -((n - p - 1) * ln_det_S) - (m * log_determinant_spd(Omega)) +
      ((m + n) * log_determinant_spd((m * Omega + n * S) / (m + n)));
    real ll = (F_1 + F_2) / -2.0;
    return(ll);
  }
}
data {
  int<lower = 0> Np;  // number persons
  int<lower = 0> Ni;  // number items
  matrix[Ni, Ni] S;  // covariance matrix
  int<lower = 0> Nf;  // number factors
  int<lower = 0> Nce;  // number correlated errors
  array[Nce, 2] int error_mat;  // correlated error matrix
  matrix[Ni, Nf] loading_pattern;  // loading pattern
  array[Nf] int markers;  // marker variables
  int<lower = 0, upper = 1> corr_fac;  // 1 for correlated factors, 0 otherwise
  real<lower = 1> shape_phi_c;  // lkj prior shape for phi
  real<lower = 0> sl_par;  // sigma_loading parameter
  real<lower = 0> rs_par;  // residual sd parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  int<lower = 0, upper = 1> complex_struc;
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  int Nf_corr = corr_fac == 1 ? Nf : 1;
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int N_complex = 0;
  real ln_det_S = log_determinant_spd(S);
  int N_Sigma = 1;

  if (method >= 90) {
    Nisqd2 = 0;
  }

  if (method == 91 || method == 92) {
    N_Sigma = Ni;
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  for (i in 1:Ni) {
    for (j in 1:Nf) {
      if (loading_pattern[i, j] == 1) Nl += 1;
    }
  }

  if (complex_struc == 1) {
    N_complex = Ni * Nf - Nl;
  }
}
parameters {
  vector<lower = 0.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;  // residual vector
  vector[Nl] loadings;  // loadings
  real<lower = 0> sigma_loadings;  // sd of loadings, hyperparm
  vector<lower = 0>[Ni] res_sds;  // item residual sds heteroskedastic
  cholesky_factor_corr[Nf_corr] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce] res_cor_01;  // correlated errors on 01
  vector[N_complex] loadings_complex;
  vector<lower = 0>[complex_struc] sigma_loadings_complex;
  vector<lower = 2.0>[complex_struc] gdp_loadings_complex;
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

  if (complex_struc == 1) {
    sigma_loadings_complex ~ std_normal();
    gdp_loadings_complex ~ lognormal(1, 1);
    target += generalized_double_pareto_lpdf(
      loadings_complex | gdp_loadings_complex[1]);
  }

  loadings ~ normal(0, sigma_loadings);
  sigma_loadings ~ student_t(3, 0, sl_par);
  res_sds ~ student_t(3, 0, rs_par);
  phi_mat_chol ~ lkj_corr_cholesky(shape_phi_c);
  res_cor_01 ~ beta(rc_par, rc_par);

  {
    real m;
    matrix[Ni, Ni] Omega;
    vector[Ni] res_var = square(res_sds);
    matrix[Nf_corr, Nf_corr] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol);
    vector[Nce] res_cor = res_cor_01 * 2 - 1;
    matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
    matrix[Ni, Ni] lamb_phi_lamb;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;

    {
      int pos = 0;
      int pos_complex = 0;
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            pos += 1;
            Load_mat[i, j] = loadings[pos];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = sigma_loadings_complex[1] * loadings_complex[pos_complex];
          }
        }
      }
    }

    if (corr_fac == 1) lamb_phi_lamb = quad_form_sym(phi_mat, Load_mat');
    else lamb_phi_lamb = tcrossprod(Load_mat);

    for (i in 1:Nce) {
      loading_par_exp[error_mat[i, 1], i] = sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 1]]);
      loading_par_exp[error_mat[i, 2], i] = sign(res_cor[i]) * sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 2]]);
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    Omega = add_diag(lamb_phi_lamb + loading_par_exp_2, delta_mat_ast);

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

    if (method != 91 && method != 92) {
      Sigma ~ inv_wishart(1000, identity_matrix(1));
    }

    if (method >= 90 && method <= 99) {
      m = 1.0 / square(rms_src_p[1]) + Ni - 1;
      if (method == 90) {
        target += gen_matrix_beta_ii_lpdf(S | Omega, Np - 1.0, m, ln_det_S);
      } else if (method == 91) {
        Sigma ~ inv_wishart(m, m * Omega);
        target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Sigma));
      } else if (method == 92) {
        Sigma ~ wishart(m, Omega / m);
        target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Sigma));
      }
    } else {
      target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Omega));
    }
  }
}
generated quantities {
  real D_obs;
  real D_rep;
  real<lower = 0, upper = 1> ppp;
  real<lower = 0> rms_src;  // RMSE of residuals
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf_corr, Nf_corr] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol);
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
    matrix[Ni, Ni] lamb_phi_lamb;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;

    {
      int pos = 0;
      int pos_complex = 0;
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            pos += 1;
            Load_mat[i, j] = loadings[pos];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = sigma_loadings_complex[1] * loadings_complex[pos_complex];
          }
        }
      }
    }

    if (corr_fac == 1) lamb_phi_lamb = quad_form_sym(phi_mat, Load_mat');
    else lamb_phi_lamb = tcrossprod(Load_mat);

    for (i in 1:Nce) {
      loading_par_exp[error_mat[i, 1], i] = sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 1]]);
      loading_par_exp[error_mat[i, 2], i] = sign(res_cor[i]) * sqrt(
        abs(res_cor[i]) * res_var[error_mat[i, 2]]);
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    Omega = add_diag(lamb_phi_lamb + loading_par_exp_2, delta_mat_ast);

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
      } else if (method == 91 || method == 92) {
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
  }

  for (j in 1:Nf) {
    if (Load_mat[markers[j], j] < 0) {
      Load_mat[, j] *= -1.0;
      if (corr_fac == 1) {
        phi_mat[, j] *= -1.0;
        phi_mat[j, ] *= -1.0;
      }
      for (i in Load_mat[, j]) {
        if (i != 0) {
          for (k in 1:Nce) {
            if (i == error_mat[k, 1]) res_cor[k] *= -1.0;
            if (i == error_mat[k, 2]) res_cor[k] *= -1.0;
          }
        }
      }
    }
  }

  for (i in 1:Nce) {
    res_cov[i] = res_cor[i] * prod(res_sds[error_mat[i, ]]);
  }
}

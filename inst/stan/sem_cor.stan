functions {
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
  real generalized_double_pareto_lpdf(vector x, real alpha, real scale) {
    // generalized double Pareto
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/
    return(sum(
      -log(2) - log(scale) - (alpha + 1.0) * log1p(abs(x) / (scale * alpha))
    ));
  }
  real eff(int p, real x) {
    return(
      2 * lmgamma(p, x / 2) - x * p * log(x / 2) + x * p
    );
  }
  vector matrix_log_vech(matrix A) {
    int p = rows(A);
    int p_ast = (p * (p - 1)) %/% 2;
    vector[p] lambda_tmp = eigenvalues_sym(A);
    matrix[p, p] q_mat_tmp = eigenvectors_sym(A);
    vector[p] ln_lambda;
    matrix[p, p] q_mat;
    matrix[p, p] ln_A;
    vector[p_ast] ln_A_lower_tri;

    for (i in 1:p) {
      ln_lambda[i] = log(lambda_tmp[p - i + 1]);
      q_mat[, i] = q_mat_tmp[, p - i + 1];
    }

    ln_A = quad_form(diag_matrix(ln_lambda), q_mat');

    {
      int pos = 0;
      for (j in 1:(p - 1)) {
        for (i in (j + 1):p) {
          pos += 1;
          ln_A_lower_tri[pos] = ln_A[i, j];
        }
      }
    }

    return(ln_A_lower_tri);
  }
}
data {
  int<lower = 0> Np;  // number persons
  int<lower = 0> Ni;  // number items
  vector[(Ni * (Ni - 1)) %/% 2] r_obs_vec;  // correlation matrix
  matrix[(Ni * (Ni - 1)) %/% 2, (Ni * (Ni - 1)) %/% 2] r_obs_vec_cov;  // ACOV(correlation matrix)
  int<lower = 0> Nf;  // number factors
  int<lower = 0> Nce;  // number correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni, Nf] int loading_pattern;  // loading pattern
  array[Ni, Nf] real loading_fixed;  // loading pattern
  array[Ni] int res_var_pattern;  // res_var pattern
  array[Ni] real<lower = 0.0> res_var_fixed;  // res_var pattern
  array[Nf, Nf] int coef_pattern;  // coef pattern
  array[Nf, Nf] real coef_fixed;  // coef fixed
  array[Nf] int markers;  // marker variables
  matrix[Nf, Nf] corr_mask;  // 1 for correlated factors, 0 otherwise
  real<lower = 1> shape_phi_c;  // lkj prior shape for phi
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 0> rs_par;  // residual sd parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  int<lower = 0, upper = 1> complex_struc;
  matrix[Ni, Nf] load_est;
  matrix[Ni, Nf] load_se;
  matrix[Nf, Nf] coef_est;
  matrix[Nf, Nf] coef_se;
  int <lower = 0, upper = 1> ret_ll; // return log-likelihood?
  int <lower = 0, upper = 1> has_data; // has data?
  matrix[Np * has_data, Ni] Y; // data
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl_uniq = 0;  // N_non-zero loadings unique
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  int<lower = 0> Nl_fixed = 0;  // N_non-zero loadings unique
  int<lower = 0> Nrv_uniq = 0;  // N_non-zero res_var unique
  int<lower = 0> Nrv = 0;  // N_non-zero res_var
  int<lower = 0> Nrv_fixed = 0;  // N_non-zero res_var unique
  int<lower = 0> Nco_uniq = 0;  // N_non-zero coef unique
  int<lower = 0> Nco = 0;  // N_non-zero coef
  int<lower = 0> Nco_fixed = 0;  // N_non-zero coef unique
  cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int N_complex = 0;
  real ln_det_S = log_determinant_spd(S);
  int N_Sigma = 1;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique
  int Np_ll = Np * ret_ll * has_data;

  if (method >= 90) {
    Nisqd2 = 0;
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  for (i in 1:Ni) {
    for (j in 1:Nf) {
      if (loading_pattern[i, j] != 0) {
        Nl += 1;
        if (loading_pattern[i, j] > Nl_uniq) Nl_uniq = loading_pattern[i, j];
      } else if (loading_fixed[i, j] != -999) {
        Nl_fixed += 1;
      }
    }
  }

  if (complex_struc == 1) {
    N_complex = Ni * Nf - Nl - Nl_fixed;
  }

  for (i in 1:Nf) {
    for (j in 1:Nf) {
      if (coef_pattern[i, j] != 0) {
        Nco += 1;
        if (coef_pattern[i, j] > Nco_uniq) Nco_uniq = coef_pattern[i, j];
      } else if (coef_fixed[i, j] != -999) {
        Nco_fixed += 1;
      }
    }
  }

  for (j in 1:(Ni - 1)) {
    for (i in (j + 1):Ni) {
      if (error_pattern[i, j] > Nce_uniq) {
        Nce_uniq = error_pattern[i, j];
      }
    }
  }

  for (i in 1:Ni) {
    if (res_var_pattern[i] != 0) {
      Nrv += 1;
      if (res_var_pattern[i] > Nrv_uniq) Nrv_uniq = res_var_pattern[i];
    } else if (res_var_fixed[i] != 999) {
      Nrv_fixed += 1;
    }
  }
}
parameters {
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;  // residual vector
  vector[Nl_uniq] loadings;  // loadings
  vector<lower = 0>[Nrv_uniq] res_sds_u;  // item residual sds heteroskedastic
  cholesky_factor_corr[Nf] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector[N_complex] loadings_complex;
  vector<lower = -1, upper = 1>[Nco_uniq] coefs;  // may need to be limited to (-1, 1)?
  vector<lower = 0>[complex_struc] sigma_loadings_complex;
  vector<lower = 2.0>[complex_struc] gdp_loadings_complex;
  cov_matrix[N_Sigma] Sigma;
}
transformed parameters {
  real rms_src_tmp = 0.0;

  if (method != 100) rms_src_tmp = rms_src_p[1];
  if (method == 2) {
    rms_src_tmp /= sqrt_two;
  } else if (method == 3) {
    rms_src_tmp /= pi_sqrt_three;
  } else if (method == 4) {
    rms_src_tmp /= sqrt_two * gdp_alpha[1] / sqrt(
      (gdp_alpha[1] - 1.0) * (gdp_alpha[1] - 2.0)
    );
  }
}
model {
  rms_src_p ~ normal(0, rm_par);
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
      resids | gdp_alpha[1], 1.0);
  }

  if (complex_struc == 1) {
    sigma_loadings_complex ~ normal(0, 0.25);
    gdp_loadings_complex ~ lognormal(1, 1);
    target += generalized_double_pareto_lpdf(
      loadings_complex | gdp_loadings_complex[1], sigma_loadings_complex[1]);
  }

  res_sds_u ~ student_t(3, 0, rs_par);
  phi_mat_chol ~ lkj_corr_cholesky(shape_phi_c);
  res_cor_01 ~ beta(rc_par, rc_par);

  {
    real m;
    matrix[Ni, Ni] Omega;
    vector[Ni] res_var;
    matrix[Nf, Nf] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol) .* corr_mask;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
    matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
    matrix[Nf, Nf] One_min_Beta_inv;
    matrix[Nf, Nf] imp_phi;
    matrix[Ni, Ni] lamb_phi_lamb;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    vector[Nf] phi_sd;

    for (i in 1:Nf) {
      for (j in 1:Nf) {
        if (coef_pattern[i, j] != 0) {
          coefs[coef_pattern[i, j]] ~ normal(coef_est[i, j], coef_se[i, j]);
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] != -999) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    {
      int pos_complex = 0;
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            loadings[loading_pattern[i, j]] ~ normal(load_est[i, j], load_se[i, j]);
            Load_mat[i, j] = loadings[loading_pattern[i, j]];
          } else if (loading_fixed[i, j] != -999) {
            Load_mat[i, j] = loading_fixed[i, j];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = loadings_complex[pos_complex];
          }
        }
      }
    }

    for (i in 1:Nf) {
      phi_sd[i] = fmax(0.0, sqrt(1 - quad_form_sym(phi_mat, Coef_mat[i, ]')));
    }
    One_min_Beta_inv = inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);
    imp_phi = quad_form_sym(quad_form_diag(phi_mat, phi_sd), One_min_Beta_inv');
    lamb_phi_lamb = quad_form_sym(imp_phi, Load_mat');

    for (i in 1:Ni) {
      if (res_var_pattern[i] != 0) {
        res_var[i] = square(res_sds_u[res_var_pattern[i]]);
      } else if (res_var_fixed[i] != 999) {
        res_var[i] = res_var_fixed[i];
      }
    }

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            loading_par_exp[i, pos_err] = sqrt(
              abs(res_cor_u[error_pattern[i, j]]) * res_var[i]
            );
            loading_par_exp[j, pos_err] =
              sign(res_cor_u[error_pattern[i, j]]) * sqrt(
                abs(res_cor_u[error_pattern[i, j]]) * res_var[j]
              );
          }
        }
      }
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
          Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    if (method != 91 && method != 92) {
      Sigma ~ inv_wishart(1000, identity_matrix(1));
    }

    if (method >= 90 && method <= 99) {
      m = 1.0 / square(rms_src_tmp) + Ni - 1;
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
  real<lower = 0> rms_src = 0.0;  // RMSE of residuals
  matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  vector[Nf] r_square = rep_vector(0, Nf);
  matrix[Nf, Nf] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol) .* corr_mask;
  vector[Ni] res_sds;
  vector[Ni] res_var;
  vector[Nce] res_cor;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);
  matrix[Ni, Ni] Omega;
  vector[Np_ll] log_lik;

  if (method != 100) rms_src = rms_src_p[1];

  if (method < 90) {
    int pos = 0;
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        pos += 1;
        Resid[i, j] = resids[pos] * rms_src_tmp;
        Resid[j, i] = Resid[i, j];
      }
    }
  }

  {
    real m;
    matrix[Ni, Ni] Sigma_p;
    matrix[Ni, Ni] S_sim;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Nf, Nf] One_min_Beta_inv;
    matrix[Nf, Nf] imp_phi;
    matrix[Ni, Ni] lamb_phi_lamb;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    vector[Nf] phi_sd;

    for (i in 1:Nf) {
      for (j in 1:Nf) {
        if (coef_pattern[i, j] != 0) {
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] != -999) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    {
      int pos_complex = 0;
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            Load_mat[i, j] = loadings[loading_pattern[i, j]];
          } else if (loading_fixed[i, j] != -999) {
            Load_mat[i, j] = loading_fixed[i, j];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = loadings_complex[pos_complex];
          }
        }
      }
    }

    for (i in 1:Nf) {
      phi_sd[i] = fmax(0.0, sqrt(1 - quad_form_sym(phi_mat, Coef_mat[i, ]')));
    }
    One_min_Beta_inv = inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);
    imp_phi = quad_form_sym(quad_form_diag(phi_mat, phi_sd), One_min_Beta_inv');
    lamb_phi_lamb = quad_form_sym(imp_phi, Load_mat');

    for (i in 1:Ni) {
      if (res_var_pattern[i] != 0) {
        res_var[i] = square(res_sds_u[res_var_pattern[i]]);
      } else if (res_var_fixed[i] != 999) {
        res_var[i] = res_var_fixed[i];
      }
    }

    res_sds = sqrt(res_var);

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            res_cor[pos_err] = res_cor_u[error_pattern[i, j]];
            res_cov[pos_err] = res_cor[pos_err] * res_sds[i] * res_sds[j];
            loading_par_exp[i, pos_err] = sqrt(
              abs(res_cor_u[error_pattern[i, j]]) * res_var[i]
            );
            loading_par_exp[j, pos_err] =
              sign(res_cor_u[error_pattern[i, j]]) * sqrt(
                abs(res_cor_u[error_pattern[i, j]]) * res_var[j]
              );
          }
        }
      }
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
          Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    if (ret_ll == 1) {
      vector[Ni] zero_vec = rep_vector(0.0, Ni);
      if (method == 91 || method == 92) {
        for (i in 1:Np_ll) {
          log_lik[i] = multi_normal_lpdf(Y[i, ] | zero_vec, Sigma);
        }
      } else {
        for (i in 1:Np_ll) {
          log_lik[i] = multi_normal_lpdf(Y[i, ] | zero_vec, Omega);
        }
      }
    }

    if (method >= 90 && method <= 99) {
      if (method == 90) {
        m = 1.0 / square(rms_src_tmp) + Ni - 1;
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

    r_square = 1.0 - square(phi_sd);
  }

  for (j in 1:Nf) {
    if (markers[j] != 0) {
      if (Load_mat[markers[j], j] < 0) {
        Load_mat[, j] *= -1.0;
        phi_mat[, j] *= -1.0;
        phi_mat[j, ] *= -1.0;
        Coef_mat[, j] *= -1.0;
        Coef_mat[j, ] *= -1.0;
      }
    }
  }
}

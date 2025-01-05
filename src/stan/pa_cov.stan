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
  int<lower = 0> Nce;  // number correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni] int res_var_pattern;  // res_var pattern
  array[Ni] real<lower = 0.0> res_var_fixed;  // res_var pattern
  array[Ni, Ni] int coef_pattern;  // coef pattern
  array[Ni, Ni] real coef_fixed;  // coef fixed
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 0> rs_par;  // residual sd parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  array[Ni, Ni] int cond_ind_mat; // conditional independence locations
  matrix[Ni, Ni] coef_est;
  matrix[Ni, Ni] coef_se;
  int <lower = 0, upper = 1> ret_ll; // return log-likelihood?
  int <lower = 0, upper = 1> has_data; // has data?
  matrix[Np * has_data, Ni] Y; // data
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nrv_uniq = 0;  // N_non-zero res_var unique
  int<lower = 0> Nrv = 0;  // N_non-zero res_var
  int<lower = 0> Nrv_fixed = 0;  // N_non-zero res_var unique
  int<lower = 0> Nco_uniq = 0;  // N_non-zero coef unique
  int<lower = 0> Nco = 0;  // N_non-zero coef
  int<lower = 0> Nco_fixed = 0;  // N_non-zero coef unique
  cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  int Ncond_ind = 0;
  int N_rms = 1;
  int N_alpha = 0;
  real ln_det_S = log_determinant_spd(S);
  int is_wb_adj = 0;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique
  int Np_ll = Np * ret_ll * has_data;

  if (method < 90) {
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        if (cond_ind_mat[i, j] == 1) {
          Ncond_ind += 1;
        }
      }
    }
  }

  if (method == 91 || method == 92) {
    is_wb_adj = 1;
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  for (i in 1:Ni) {
    for (j in 1:Ni) {
      if (coef_pattern[i, j] != 0) {
        Nco += 1;
        if (coef_pattern[i, j] > Nco_uniq) Nco_uniq = coef_pattern[i, j];
      } else if (coef_fixed[i, j] > -990) {
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
    } else if (res_var_fixed[i] < 990) {
      Nrv_fixed += 1;
    }
  }
}
parameters {
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Ncond_ind] resids;  // residual vector
  vector<lower = 0>[Nrv_uniq] res_sds_u;  // item residual sds heteroskedastic
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector<lower = -1, upper = 1>[Nco_uniq] coefs;  // may need to be limited to (-1, 1)?
  array[is_wb_adj] cov_matrix[Ni] Sigma;
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

  res_sds_u ~ student_t(3, 0, rs_par);
  res_cor_01 ~ beta(rc_par, rc_par);

  {
    real m;
    matrix[Ni, Ni] Omega;
    vector[Ni] res_var;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] Coef_mat = rep_matrix(0, Ni, Ni);
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;

    for (i in 1:Ni) {
      for (j in 1:Ni) {
        if (coef_pattern[i, j] != 0) {
          coefs[coef_pattern[i, j]] ~ normal(coef_est[i, j], coef_se[i, j]);
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] > -990) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    One_min_Beta_inv = inverse(identity_matrix(Ni) - Coef_mat);

    for (i in 1:Ni) {
      if (res_var_pattern[i] != 0) {
        res_var[i] = square(res_sds_u[res_var_pattern[i]]);
      } else if (res_var_fixed[i] < 990) {
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

    Omega = quad_form_sym(
      add_diag(loading_par_exp_2, delta_mat_ast),
      One_min_Beta_inv'
    );

    total_var = diagonal(Omega);

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          if (cond_ind_mat[i, j] == 1) {
             pos += 1;
             Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
             Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    if (method >= 90 && method <= 99) {
      m = 1.0 / square(rms_src_tmp) + Ni - 1;
      if (method == 90) {
        target += gen_matrix_beta_ii_lpdf(S | Omega, Np - 1.0, m, ln_det_S);
      } else if (method == 91) {
        Sigma[1] ~ inv_wishart(m, m * Omega);
        target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Sigma[1]));
      } else if (method == 92) {
        Sigma[1] ~ wishart(m, Omega / m);
        target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Sigma[1]));
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
  matrix[Ni, Ni] Coef_mat = rep_matrix(0, Ni, Ni);
  vector[Ni] r_square = rep_vector(0, Ni);
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
        if (cond_ind_mat[i, j] == 1) {
          pos += 1;
          Resid[i, j] = resids[pos] * rms_src_tmp;
          Resid[j, i] = Resid[i, j];
        }
      }
    }
    if (Ncond_ind == 0) rms_src = 0;
  }

  {
    real m;
    matrix[Ni, Ni] Sigma_p;
    matrix[Ni, Ni] S_sim;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;

    for (i in 1:Ni) {
      for (j in 1:Ni) {
        if (coef_pattern[i, j] != 0) {
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] > -990) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    One_min_Beta_inv = inverse(identity_matrix(Ni) - Coef_mat);

    for (i in 1:Ni) {
      if (res_var_pattern[i] != 0) {
        res_var[i] = square(res_sds_u[res_var_pattern[i]]);
      } else if (res_var_fixed[i] < 990) {
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

    Omega = quad_form_sym(
      add_diag(loading_par_exp_2, delta_mat_ast),
      One_min_Beta_inv'
    );

    total_var = diagonal(Omega);
    r_square = 1.0 - res_var ./ total_var;

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          if (cond_ind_mat[i, j] == 1) {
            pos += 1;
            Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
            Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    if (ret_ll == 1) {
      vector[Ni] zero_vec = rep_vector(0.0, Ni);
      if (method == 91 || method == 92) {
        for (i in 1:Np_ll) {
          log_lik[i] = multi_normal_lpdf(Y[i, ] | zero_vec, Sigma[1]);
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
        S_sim = wishart_rng(Np - 1.0, Sigma[1] / (Np - 1.0));
        D_obs = -2.0 * wishart_lpdf(S | Np - 1.0, Sigma[1]);
        D_rep = -2.0 * wishart_lpdf(S_sim | Np - 1.0, Sigma[1]);
      }
    } else {
      S_sim = wishart_rng(Np - 1.0, Omega / (Np - 1.0));
      D_obs = -2.0 * wishart_lpdf(S | Np - 1.0, Omega / (Np - 1.0));
      D_rep = -2.0 * wishart_lpdf(S_sim | Np - 1.0, Omega / (Np - 1.0));
    }
    ppp = D_rep > D_obs ? 1.0 : 0.0;
  }
}

functions {
  array[] int find_row_zero(matrix sq_mat, array[] int is_row_fixed) {
    int ni = rows(sq_mat);
    int des_zero = ni - sum(is_row_fixed);
    array[ni] int tmp_res = rep_array(0, ni);
    int idx = 0;
    for (i in 1:ni) {
      if (is_row_fixed[i] == 0) {
        int temp = 0;
        for (j in 1:ni) {
          if (is_row_fixed[j] == 0) {
            if (sq_mat[i, j] == 0) temp += 1;
          }
        }
        if (temp == des_zero) {
          idx += 1;
          tmp_res[idx] = i;
        }
      }
    }
    array[idx] int result = tmp_res[1:idx];
    return(result);
  }
  array[,] int find_recursive_set(matrix coef_mat) {
    int ni = rows(coef_mat);
    array[ni, ni] int result = rep_array(0, ni, ni);
    array[ni] int fix_variable = rep_array(0, ni);
    int ni_sofar = 0;
    int i = 0;
    while (ni_sofar < ni) {
      array[size(find_row_zero(coef_mat, fix_variable))] int temp = find_row_zero(coef_mat, fix_variable);
      fix_variable[temp] = rep_array(1, size(temp));
      i += 1;
      result[i, 1:size(temp)] = temp;
      ni_sofar += size(temp);
    }
    return(result);
  }
  vector find_factor_res_var(matrix coef_mat, matrix cor_psi) {
    int ni = rows(coef_mat);
    array[ni, ni] int full_set = find_recursive_set(coef_mat);
    vector[ni] error_var = rep_vector(1.0, ni);
    vector[ni] total_var_psi = rep_vector(1.0, ni);
    array[ni] int iv = rep_array(0, ni);

    int n_set_1 = 0;
    for (i in 1:ni) {
      if (full_set[1, i] != 0) n_set_1 += 1;
    }
    array[n_set_1] int idx_set_1 = rep_array(0, n_set_1);
    {
      int counter = 0;
      for (i in 1:ni) {
        if (full_set[1, i] != 0) {
          counter += 1;
          idx_set_1[counter] = full_set[1, i];
        }
      }
    }
    error_var[idx_set_1] = total_var_psi[idx_set_1];

    matrix[n_set_1, n_set_1] ic_cor = cor_psi[idx_set_1, idx_set_1];
    vector[n_set_1] start_var = total_var_psi[idx_set_1];
    matrix[ni, ni] iv_cov;
    iv_cov[1:n_set_1, 1:n_set_1] = quad_form_diag(ic_cor, sqrt(start_var));

    int n_sets = 0;
    for (i in 1:ni) {
      if (sum(full_set[i, ]) > 0) n_sets += 1;
    }
    array[n_sets, ni] int set = full_set[1:n_sets, ];

    int n_set_start;
    int n_set_end = 0;
    for (i in 1:(n_sets - 1)) {
      n_set_start = n_set_end + 1;
      int n_set_i = 0;
      int n_set_ip1 = 0;
      for (j in 1:ni) {
        if (set[i, j] != 0) n_set_i += 1;
        if (set[i + 1, j] != 0) n_set_ip1 += 1;
      }
      n_set_end += n_set_i;
      array[n_set_i] int idx_set_i = rep_array(0, n_set_i);
      array[n_set_ip1] int idx_set_ip1 = rep_array(0, n_set_ip1);
      {
        int counter_i = 0;
        int counter_ip1 = 0;
        for (j in 1:ni) {
          if (set[i, j] != 0) {
            counter_i += 1;
            idx_set_i[counter_i] = set[i, j];
          }
          if (set[i + 1, j] != 0) {
            counter_ip1 += 1;
            idx_set_ip1[counter_ip1] = set[i + 1, j];
          }
        }
      }
      iv[n_set_start:n_set_end] = idx_set_i;

      matrix[n_set_ip1, n_set_end] tmp_beta = coef_mat[idx_set_ip1, iv[1:n_set_end]];

      matrix[n_set_ip1, n_set_ip1] var_reg = quad_form_sym(
        iv_cov[1:n_set_end, 1:n_set_end], tmp_beta'
      );

      matrix[n_set_ip1, n_set_ip1] temp_psi = cor_psi[idx_set_ip1, idx_set_ip1];
      vector[n_set_ip1] temp_psi_sd = rep_vector(0.0, n_set_ip1);

      for (j in 1:n_set_ip1) {
        error_var[idx_set_ip1[j]] = total_var_psi[idx_set_ip1[j]] - var_reg[j, j];
        temp_psi_sd[j] = fmax(0.0, sqrt(error_var[idx_set_ip1[j]]));
      }

      if (i < (n_sets - 1)) {
        temp_psi = quad_form_diag(temp_psi, temp_psi_sd);
        int n_agg = n_set_end + n_set_ip1;
        matrix[n_agg, n_agg] real_temp_psi = rep_matrix(0, n_agg, n_agg);
        real_temp_psi[1:n_set_end, 1:n_set_end] = iv_cov[1:n_set_end, 1:n_set_end];
        real_temp_psi[(n_set_end + 1):n_agg, (n_set_end + 1):n_agg] = temp_psi;
        array[n_agg] int agg;
        agg[1:n_set_end] = iv[1:n_set_end];
        agg[(n_set_end + 1):n_agg] = idx_set_ip1;
        matrix[n_agg, n_agg] temp_path2 = rep_matrix(0, n_agg, n_agg);
        temp_path2[(n_set_end + 1):n_agg, ] = coef_mat[idx_set_ip1, agg];
        matrix[n_agg, n_agg] id_mat = identity_matrix(n_agg);
        matrix[n_agg, n_agg] id_tmp2_inv = inverse(id_mat - temp_path2);
        iv_cov[1:n_agg, 1:n_agg] = quad_form(real_temp_psi, id_tmp2_inv');
      }
    }

    return(error_var);
  }
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
  int is_wb_adj = 0;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique
  int Np_ll = Np * ret_ll * has_data;

  if (method >= 90) {
    Nisqd2 = 0;
  }

  if (method == 91 || method == 92) {
    is_wb_adj = 1;
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
      } else if (loading_fixed[i, j] > -990) {
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
  vector[Nisqd2] resids;  // residual vector
  vector[Nl_uniq] loadings;  // loadings
  vector<lower = 0>[Nrv_uniq] res_sds_u;  // item residual sds heteroskedastic
  cholesky_factor_corr[Nf] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector[N_complex] loadings_complex;
  vector<lower = -1, upper = 1>[Nco_uniq] coefs;  // may need to be limited to (-1, 1)?
  vector<lower = 0>[complex_struc] sigma_loadings_complex;
  vector<lower = 2.0>[complex_struc] gdp_loadings_complex;
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
        } else if (coef_fixed[i, j] > -990) {
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
          } else if (loading_fixed[i, j] > -990) {
            Load_mat[i, j] = loading_fixed[i, j];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = loadings_complex[pos_complex];
          }
        }
      }
    }

    phi_sd = fmax(0.0, sqrt(find_factor_res_var(Coef_mat, phi_mat)));
    One_min_Beta_inv = inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);
    imp_phi = quad_form_sym(quad_form_diag(phi_mat, phi_sd), One_min_Beta_inv');
    lamb_phi_lamb = quad_form_sym(imp_phi, Load_mat');

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
        } else if (coef_fixed[i, j] > -990) {
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
          } else if (loading_fixed[i, j] > -990) {
            Load_mat[i, j] = loading_fixed[i, j];
          } else if (complex_struc == 1) {
            pos_complex += 1;
            Load_mat[i, j] = loadings_complex[pos_complex];
          }
        }
      }
    }

    r_square = 1.0 - find_factor_res_var(Coef_mat, phi_mat);
    phi_sd = sqrt(fmax(0.0, 1.0 - r_square));
    One_min_Beta_inv = inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);
    imp_phi = quad_form_sym(quad_form_diag(phi_mat, phi_sd), One_min_Beta_inv');
    lamb_phi_lamb = quad_form_sym(imp_phi, Load_mat');

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

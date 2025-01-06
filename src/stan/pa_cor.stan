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
  int<lower = 0> Nce;  // number correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni, Ni] int coef_pattern;  // coef pattern
  array[Ni, Ni] real coef_fixed;  // coef fixed
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  array[Ni, Ni] int cond_ind_mat; // conditional independence locations
  matrix[Ni, Ni] coef_est;
  matrix[Ni, Ni] coef_se;
  int<lower = 0, upper = 1> centered;
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nco_uniq = 0;  // N_non-zero coef unique
  int<lower = 0> Nco = 0;  // N_non-zero coef
  int<lower = 0> Nco_fixed = 0;  // N_non-zero coef unique
  int Ncond_ind = 0;
  int Nisqd2_vec = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int is_wb_adj = 0;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique
  matrix[Nisqd2_vec, Nisqd2_vec] L_vec_cov = cholesky_decompose(r_obs_vec_cov);

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
}
parameters {
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Ncond_ind] resids;  // residual vector
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector<lower = -1, upper = 1>[Nco_uniq] coefs;  // may need to be limited to (-1, 1)?
  array[is_wb_adj] vector[Nisqd2_vec] r_vec_ri;
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

  res_cor_01 ~ beta(rc_par, rc_par);

  {
    matrix[Ni, Ni] Omega;
    vector[Nisqd2_vec] r_vec;
    vector[Ni] res_var;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] Coef_mat = rep_matrix(0, Ni, Ni);
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    matrix[Ni, Ni] cor_psi = identity_matrix(Ni);

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

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            cor_psi[i, j] = res_cor_u[error_pattern[i, j]];
            cor_psi[j, i] = cor_psi[i, j];
          }
        }
      }
    }

    // fix this, error correlation matrix is NaN
    res_var = find_factor_res_var(Coef_mat, cor_psi);

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

    r_vec = matrix_log_vech(Omega);

    {
      vector[Nisqd2_vec] tmp_loc = r_vec;
      matrix[Nisqd2_vec, Nisqd2_vec] tmp_cov = L_vec_cov;

      if (method >= 90 && method <= 99) {
        if (method == 90) {
          tmp_cov = cholesky_decompose(add_diag(r_obs_vec_cov, square(rms_src_tmp)));
        } else if (method == 91 || method == 92) {
          if (centered == 0) {
            r_vec_ri[1] ~ std_normal();
            tmp_loc = r_vec + rms_src_tmp * r_vec_ri[1];
          } else if (centered == 1) {
            r_vec_ri[1] ~ normal(r_vec, rms_src_tmp);
            tmp_loc = r_vec_ri[1];
          }
        }
      }

      target += multi_normal_cholesky_lupdf(r_obs_vec | tmp_loc, tmp_cov);
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
    vector[Nisqd2_vec] r_vec;
    vector[Nisqd2_vec] r_vec_sim;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    matrix[Ni, Ni] cor_psi = identity_matrix(Ni);

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

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            cor_psi[i, j] = res_cor_u[error_pattern[i, j]];
            cor_psi[j, i] = cor_psi[i, j];
          }
        }
      }
    }

    res_var = find_factor_res_var(Coef_mat, cor_psi);
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

    r_vec = matrix_log_vech(Omega);

    {
      vector[Nisqd2_vec] tmp_loc = r_vec;
      matrix[Nisqd2_vec, Nisqd2_vec] tmp_cov = L_vec_cov;

      if (method >= 90 && method <= 99) {
        if (method == 90) {
          tmp_cov = cholesky_decompose(add_diag(r_obs_vec_cov, square(rms_src_tmp)));
        } else if (method == 91 || method == 92) {
          if (centered == 0) {
            tmp_loc = r_vec + rms_src_tmp * r_vec_ri[1];
          } else if (centered == 1) {
            tmp_loc = r_vec_ri[1];
          }
        }
      }

      r_vec_sim = multi_normal_cholesky_rng(tmp_loc, tmp_cov);
      D_obs = -2.0 * multi_normal_cholesky_lpdf(r_obs_vec | tmp_loc, tmp_cov);
      D_rep = -2.0 * multi_normal_cholesky_lpdf(r_vec_sim | tmp_loc, tmp_cov);
    }
    ppp = D_rep > D_obs ? 1.0 : 0.0;
  }
}

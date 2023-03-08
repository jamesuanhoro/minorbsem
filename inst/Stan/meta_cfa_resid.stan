// Each S is dist according to some E
// And each E follows some structured matrix
// And that structured matrix is IW (see Wu & Browne)
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
  real gen_matrix_beta_ii_lpdf(matrix S, matrix Omega, real n, real m) {
    int p = rows(S);
    real F_1 = eff(p, m) + eff(p, n) - eff(p, m + n);
    real F_2 = -((n - p - 1) * log_determinant_spd(S)) - (m * log_determinant_spd(Omega)) +
      ((m + n) * log_determinant_spd((m * Omega + n * S) / (m + n)));
    real ll = (F_1 + F_2) / -2.0;
    return(ll);
  }
}
data {
  int<lower = 0> Ng;  // number of groups
  array[Ng] int<lower = 0> Np;  // number persons by matrix
  int<lower = 0> Ni;  // number items
  array[Ng] matrix[Ni, Ni] S;  // covariance matrices
  int<lower = 0> Nf; // N_factors
  int<lower = 0> Nce; // N_correlated errors
  array[Nce, 2] int error_mat; // cor error matrix
  matrix[Ni, Nf] loading_pattern;
  array[Nf] int markers; // markers
  int<lower = 0, upper = 1> corr_fac;  // 1 for correlated factors, 0 otherwise
  real<lower = 1> shape_phi_c; // lkj prior shape for phi
  real<lower = 0> sl_par;  // sigma_loading parameter
  real<lower = 0> rs_par;  // residual sd parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  int<lower = 0, upper = 1> complex_struc;
  int Nmiss;  // number missing correlations
  int Nitem_miss;  // number items with missing correlations across groups
  array[Ni, Ng] int valid_var;  // indicator of valid items within group
  array[Ni, Ng] int miss_ind;  // items with missing correlations
  int p;  // number of moderators
  matrix[Ng, p] X;  // moderator matrix
  real<lower = 0> mln_par;  // meta-reg int hyper-parameter
  real<lower = 0> mlb_par;  // meta-reg beta hyper-parameter
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  int Nf_corr = corr_fac == 1 ? Nf : 1;
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int N_complex = 0;
  // int N_Sigma = 1;

  if (method >= 90) {
    Nisqd2 = 0;
  }

  // if (method == 91) {
  //   N_Sigma = Ni;
  // }

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
  // cov_matrix[N_Sigma] Sigma;
  real m_ln_int;
  vector[p] m_ln_beta;
  vector<lower = 0, upper = 1>[Nmiss] miss_cor_01;
  vector<lower = 0>[Nitem_miss] var_shifts;
}
model {
  vector[Ng] m_s = exp(m_ln_int + X * m_ln_beta) + Ni - 1;

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

  m_ln_int ~ student_t(3, 0, mln_par);
  m_ln_beta ~ student_t(3, 0, mlb_par);
  var_shifts ~ std_normal();

  {
    matrix[Ni, Ni] S_impute;
    int pos_miss_c = 0;
    int pos_miss_i = 0;
    int pos_valid;
    matrix[Ni, Ni] Omega;
    vector[Ni] total_var;

    {
      vector[Ni] res_var = square(res_sds);
      matrix[Nf_corr, Nf_corr] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol);
      vector[Nce] res_cor = res_cor_01 * 2 - 1;
      matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
      matrix[Ni, Ni] lamb_phi_lamb;
      matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
      matrix[Ni, Ni] loading_par_exp_2;
      vector[Ni] delta_mat_ast;

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
    }

    for (i in 1:Ng) {
      array[sum(valid_var[, i])] int idxs;

      pos_valid = 0;
      for (j in 1:Ni) {
        if (valid_var[j, i] == 1) {
          pos_valid += 1;
          idxs[pos_valid] = j;
        }
      }

      S_impute = S[i];
      for (j in 2:Ni) {
        if (valid_var[j, i] == 1) {
          for (k in 1:(j - 1)) {
            if (valid_var[k, i] == 1) {
              if (S_impute[j, k] == 999) {
                pos_miss_c += 1;
                S_impute[j, k] = (miss_cor_01[pos_miss_c] * 2 - 1) *
                  sqrt(S_impute[j, j]) * sqrt(S_impute[k, k]);
                S_impute[k, j] = S_impute[j, k];
              }
            }
          }
        }
      }

      for (j in 1:Ni) {
        if (miss_ind[j, i] == 1) {
          pos_miss_i += 1;
          S_impute[j, j] += (var_shifts[pos_miss_i]);
        }
      }

      target += gen_matrix_beta_ii_lpdf(
        S_impute[idxs, idxs] | Omega[idxs, idxs], Np[i] - 1.0, m_s[i]);
    }
  }
}
generated quantities {
  real<lower = 0> rms_src;  // RMSE of residuals
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf_corr, Nf_corr] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol);
  vector[Ni] res_var = square(res_sds);
  vector[Nce] res_cor = res_cor_01 * 2 - 1;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);
  vector[Nmiss] miss_cor = miss_cor_01 * 2 - 1;
  vector[Ng] m_s = exp(m_ln_int + X * m_ln_beta) + Ni - 1;
  real v_mn = mean(1.0 ./ m_s);
  real rmsea_mn = sqrt(v_mn);
  vector[p] rmsea_beta;

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

  {
    vector[Ng] ebx = exp(m_ln_int + X * m_ln_beta);

    for (i in 1:p) {
      rmsea_beta[i] = -m_ln_beta[i] * mean(ebx ./ (2 * (ebx + p - 1) ^ (3.0 / 2)));
    }
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

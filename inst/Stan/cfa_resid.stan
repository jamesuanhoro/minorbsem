functions {
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
  real generalized_double_pareto_lpdf(vector x, real alpha) {
    return(sum(
      -(alpha + 1.0) * log(1.0 + abs(x) / alpha)
    ));
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
  int<lower = 1, upper = 4> method; // which method
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  int Nf_corr = corr_fac == 1 ? Nf : 1;
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_alpha = 0;

  if (method == 4) N_alpha = 1;

  for (i in 1:Ni) {
    for (j in 1:Nf) {
      if (loading_pattern[i, j] == 1) Nl += 1;
    }
  }
}
parameters {
  real<lower = 0> rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;  // residual vector
  vector[Nl] loadings;  // loadings
  real<lower = 0> sigma_loadings;  // sd of loadings, hyperparm
  vector<lower = 0>[Ni] res_sds;  // item residual sds heteroskedastic
  cholesky_factor_corr[Nf_corr] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce] res_cor_01;  // correlated errors on 01
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
    // generalized double Pareto
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/
    gdp_alpha ~ cauchy(0, 1);
    target += generalized_double_pareto_lpdf(
      resids | gdp_alpha[1]);
  }

  loadings ~ normal(0, sigma_loadings);
  sigma_loadings ~ student_t(3, 0, sl_par);
  res_sds ~ student_t(3, 0, rs_par);
  phi_mat_chol ~ lkj_corr_cholesky(shape_phi_c);
  res_cor_01 ~ beta(rc_par, rc_par);

  {
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
      for (i in 1:Ni) {
        for (j in 1:Nf) {
          if (loading_pattern[i, j] != 0) {
            pos += 1;
            Load_mat[i, j] = loadings[pos];
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

    {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          pos += 1;
          Omega[i, j] += resids[pos] * rms_src_p * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Omega));
  }
}
generated quantities {
  real<lower = 0> rms_src = rms_src_p;  // RMSE of residuals
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf_corr, Nf_corr] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol);
  vector[Ni] res_var = square(res_sds);
  vector[Nce] res_cor = res_cor_01 * 2 - 1;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);

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
    int pos = 0;
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        pos += 1;
        Resid[i, j] = resids[pos] * rms_src_p;
        Resid[j, i] = Resid[i, j];
      }
    }
  }

  {
    int pos = 0;
    for (i in 1:Ni) {
      for (j in 1:Nf) {
        if (loading_pattern[i, j] != 0) {
          pos += 1;
          Load_mat[i, j] = loadings[pos];
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

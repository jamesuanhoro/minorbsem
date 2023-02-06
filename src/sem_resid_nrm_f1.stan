functions {
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
}
data {
  int Np;
  int Ni;
  matrix[Ni, Ni] S;
  real<lower = 1> shape_beta; // beta prior shape for phi
  int Nf;
  int Nce;
  array[Nce, 2] int error_mat; // cor error matrix
  matrix[Ni, Nf] loading_pattern;
  array[Nf] int markers; // markers
  int Nf_corr;
  array[Nf_corr, 2] int F_corr_mat;
  array[Nf, Nf] int coef_pattern;
  real sc_par;  // sigma_coefficients parameter
  real sl_par;  // sigma_loading parameter
  real rs_par;  // residual sd parameter
  real rc_par;  // residual corr parameter
}
transformed data {
  int<lower = 0> Nl = 0; // N_non-zero loadings
  cholesky_factor_cov[Ni] NL_S = sqrt(Np - 1) * cholesky_decompose(S);  // covariance matrix-chol
  int coef_count = 0;
  int outcome_count = 0;
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;

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
  real<lower = 0> rms_src;
  vector[Nisqd2] resids;
  vector<lower = 0>[Nl] loadings; // loadings
  real<lower = 0> sigma_loadings; // sd of loadings, hyperparm
  // vector<lower = 0>[Nf] phi_sd;
  vector<lower = 0>[Ni] res_sds; // item residual sds heteroskedastic
  vector<lower = 0, upper = 1>[Nf_corr] phi_cor_01;
  vector<lower = 0, upper = 1>[Nce] res_cor_01; // correlated errors on 01
  vector[coef_count] coefs;
  vector<lower = 0>[outcome_count] sigma_coefs;
}
model {
  rms_src ~ std_normal();
  resids ~ std_normal();

  loadings ~ normal(0, sigma_loadings);
  sigma_loadings ~ student_t(3, 0, sl_par);
  // TODO: set prior as argument here or use standardized approach
  // phi_sd ~ student_t(3, 0, 1);

  res_sds ~ student_t(3, 0, rs_par);

  phi_cor_01 ~ beta(shape_beta, shape_beta);
  res_cor_01 ~ beta(rc_par, rc_par);

  coefs ~ std_normal();
  sigma_coefs ~ student_t(3, 0, sc_par);

  {
    matrix[Ni, Ni] Omega;
    vector[Ni] res_var = square(res_sds);
    vector[Nce] res_cor = res_cor_01 * 2 - 1;
    matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
    matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
    // vector[Nf] phi_var = square(phi_sd);
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
            pos_3[3] += 1;
            Load_mat[i, j] = loadings[pos_3[3]];
            // if (i == markers[j]) Load_mat[i, j] = 1;
            // else {
            //   pos_3[3] += 1;
            //   Load_mat[i, j] = loadings[pos_3[3]];
            // }
          }
        }
      }
    }

    Lambda_One_min_Beta_inv = Load_mat * inverse(diag_matrix(rep_vector(1, Nf)) - Coef_mat);

    for (i in 1:Nf_corr) {
      // F_corr_pe[F_corr_mat[i, 1], i] = sqrt(
      //   abs(phi_cor[i]) * phi_var[F_corr_mat[i, 1]]);
      // F_corr_pe[F_corr_mat[i, 2], i] = sign(phi_cor[i]) * sqrt(
      //   abs(phi_cor[i]) * phi_var[F_corr_mat[i, 2]]);
      F_corr_pe[F_corr_mat[i, 1], i] = sqrt(abs(phi_cor[i]));
      F_corr_pe[F_corr_mat[i, 2], i] = sign(phi_cor[i]) * sqrt(abs(phi_cor[i]));
    }

    F_cov_mat = tcrossprod(F_corr_pe);
    F_var_resid = 1.0 - diagonal(F_cov_mat);

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

    total_var = diagonal(
      quad_form(
        add_diag(F_cov_mat, F_var_resid),
        Lambda_One_min_Beta_inv')) + res_var;

    {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          pos += 1;
          Omega[i, j] += resids[pos] * rms_src * sqrt(total_var[i] * total_var[j]);
          Omega[j, i] = Omega[i, j];
        }
      }
    }

    target += wishart_cholesky_lupdf(NL_S | Np - 1, cholesky_decompose(Omega));
  }
}
generated quantities {
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf, Nf] Coef_mat = rep_matrix(0, Nf, Nf);
  // vector[Nf] phi_var = square(phi_sd);
  vector[Nf_corr] phi_cor = phi_cor_01 * 2 - 1;
  vector[Ni] res_var = square(res_sds);
  vector[Nce] res_cor = res_cor_01 * 2 - 1;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);

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
          pos_3[3] += 1;
          Load_mat[i, j] = loadings[pos_3[3]];
          // if (i == markers[j]) Load_mat[i, j] = 1;
          // else {
          //   pos_3[3] += 1;
          //   Load_mat[i, j] = loadings[pos_3[3]];
          // }
        }
      }
    }
  }

  for (j in 1:Nf) {
    if (Load_mat[markers[j], j] < 0) {
      Load_mat[, j] *= -1.0;
      Coef_mat[, j] *= -1.0;
      Coef_mat[j, ] *= -1.0;
      for (k in 1:Nf_corr) {
        if (j == F_corr_mat[k, 1]) phi_cor[k] *= -1.0;
        if (j == F_corr_mat[k, 2]) phi_cor[k] *= -1.0;
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

  {
    int pos = 0;
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        pos += 1;
        Resid[i, j] = resids[pos] * rms_src;
        Resid[j, i] = Resid[i, j];
      }
    }
  }

  for (i in 1:Nce) {
    res_cov[i] = res_cor[i] * prod(res_sds[error_mat[i, ]]);
  }
}

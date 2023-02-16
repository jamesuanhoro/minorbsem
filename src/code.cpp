// #include <Rcpp.h>
// using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]
/* https://gallery.rcpp.org/articles/dmvnorm_arma/ */
#include <RcppArmadillo.h>

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::colvec &x, arma::mat const &trimat)
{
    arma::uword const n = trimat.n_cols;

    for (unsigned j = n; j-- > 0;)
    {
        double tmp(0.);
        for (unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

// [[Rcpp::export]]
arma::mat ldmvnrm_list_arma_fast(arma::mat const &x,
                                 arma::colvec const &mean,
                                 arma::mat const &sigma_list)
{
    using arma::uword;
    uword const n = x.n_cols,
                xdim = x.n_rows,
                iter_len = sigma_list.n_cols;
    arma::mat out(iter_len, n);
    arma::mat sigma(xdim, xdim);
    arma::mat rooti(xdim, xdim);
    arma::colvec z(xdim);
    double const constants = -(double)xdim / 2.0 * log2pi;
    double rootisum,
        other_terms;

    for (uword i = 0; i < iter_len; i++)
    {
        sigma = reshape(sigma_list.col(i), xdim, xdim);
        rooti = arma::inv(trimatu(arma::chol(sigma)));
        rootisum = arma::sum(log(rooti.diag()));
        other_terms = rootisum + constants;

        for (uword j = 0; j < n; j++)
        {
            z = (x.col(j) - mean);
            inplace_tri_mat_mult(z, rooti);
            out(i, j) = other_terms - 0.5 * arma::dot(z, z);
        }
    }

    return out;
}

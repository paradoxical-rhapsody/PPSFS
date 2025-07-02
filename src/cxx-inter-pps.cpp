// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


//' @rdname interaction-effect-pps
//' @description `cxx_inter_pps`.
//' @noRd
// [[Rcpp::export]]
arma::mat cxx_inter_pps (arma::vec W, arma::vec r, const arma::mat x) {
	unsigned int p0 = x.n_cols ;
	arma::mat out(p0, p0) ; out.fill(NA_REAL) ;

	arma::vec pps(3) ;
	arma::mat V(3, 3) ;
	unsigned int iR = 0, iC = 0 ;
	for (iR = 0; iR < p0; iR++) { 
		pps(0) = arma::sum(r % x.col(iR)) ;
		V(0, 0) = arma::sum(W % x.col(iR) % x.col(iR)) ;

		for (iC = iR+1; iC < p0; iC++) {
			pps(1) = arma::sum(r % x.col(iC)) ;
			pps(2) = arma::sum(r % x.col(iR) % x.col(iC)) ;

			V(0, 1) = arma::sum(W % x.col(iR) % x.col(iC)) ;
			V(0, 2) = arma::sum(W % x.col(iR) % x.col(iR) % x.col(iC)) ;

			V(1, 1) = arma::sum(W % x.col(iC) % x.col(iC)) ;
			V(1, 2) = arma::sum(W % x.col(iR) % x.col(iC) % x.col(iC)) ;
			
			V(2, 2) = arma::sum(W % x.col(iR) % x.col(iR) % x.col(iC) % x.col(iC)) ;
			
			V = arma::symmatu(V) ;
			out(iR, iC) = arma::sum(pps % solve(V, pps)) ;
		}
	}

	return out ;
}

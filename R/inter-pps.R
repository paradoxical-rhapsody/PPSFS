#' @name interaction-effect-pps
#' @title Partial Profile Score for Composite Interaction Feature
#' 
#' @param W Vector.
#' @param r Vector.
#' @param x Matrix.
#' @return Matrix.
#' 
#' @examples
#' set.seed(2022)
#' n <- 100
#' p <- 30
#' x <- matrix(rnorm(n*p), n)
#' W <- rnorm(n)
#' r <- rnorm(n)
#' system.time( out.r <- r_inter_pps(W, r, x) )
#' system.time( out.c <- cxx_inter_pps(W, r, x) )
#' print( all.equal(out.r, out.c) )
#' 
#' @noRd
NULL
#> NULL

#' @rdname interaction-effect-pps
#' @description `r_inter_pps`.
#' @noRd
r_inter_pps <- function (W, r, x) {
	p0 <- NCOL(x)
	out <- matrix(NA_real_, p0, p0)
	for (iR in seq(1, p0-1)) { for (iC in seq(iR+1, p0)) {
		pps <- c( sum(r*x[, iR]), sum(r*x[, iC]), sum(r*x[, iR]*x[, iC]) )
		
		V <- matrix(c(
			sum(W*x[, iR]^2),	sum(W*x[, iR]*x[, iC]),	sum(W*x[, iR]^2 * x[, iC]),
			NA_real_,			sum(W*x[, iC]^2),		sum(W*x[, iR] * x[, iC]^2),
			NA_real_,			NA_real_,				sum(W*x[, iR]^2 * x[, iC]^2)
		), nrow=3, byrow=TRUE)
		V[lower.tri(V)] <- V[upper.tri(V)]
		stopifnot( isSymmetric(V) )

		out[iR, iC] <- sum(pps * solve(V, pps))
	}}

	return(out)
}

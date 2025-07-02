#' @name PPS-method
#' @title Partial Profile Score Feature Selection for GLMs
#' 
#' @param x Matrix.
#' @param y Vector.
#' @param family See [glm] and [family].
#' @param keep Initial set of features that are included in model fitting.
#' @param I0 Index set of interaction effects to be identified. 
#' @param maxK Maximum number of identified features.
#' @param ebicFlag The procedure stops when the EBIC increases after `ebicFlag` times.
#' @param ... Additional parameters for [glm.fit].
#' @param verbose Print the procedure path?
#' 
#' @return Index set of identified features.
NULL
#> NULL


#' @rdname PPS-method
#' @description `ppsfs`: PPSFS for **main-effects**.
#' @details 
#' That `ppsfs(x, y, family="gaussian")` is an implementation to 
#' *sequential lasso* method proposed by Luo and Chen(2014, <\doi{10/f6kfr6}>).
#' 
#' @examples
#' ## ***************************************************
#' ## Identify main-effect features
#' ## ***************************************************
#' set.seed(2022)
#' n <- 300
#' p <- 1000
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- eta + rnorm(n, sd=sd(eta)/5)
#' print( A <- ppsfs(x, y, 'gaussian', verbose=TRUE) )
#' 
#' @references 
#' Z. Xu, S. Luo and Z. Chen (2022). Partial profile score feature selection in 
#' high-dimensional generalized linear interaction models. Statistics and Its Interface.
#' \doi{10.4310/21-SII706}
#' 
#' @export 
ppsfs <- function (x, y, family, 
                    keep=NULL, I0=NULL, 
                    ..., 
                    ebicFlag=1,
                    maxK=min( NROW(x)-1, NCOL(x)+length(I0) ), 
                    verbose=FALSE) {
    n0 <- NROW(x)
    p0 <- NCOL(x)
    stopifnot(length(y) == n0)
    
    y  <- na.fail(y)
	x  <- na.fail(x)
    colnames(x) <- as.character(seq_len(p0))

    I0 <- unique( c(I0, grep(':', keep, value=TRUE)) )
    if (length(I0)) {
        x <- cbind(x, get_inter_features(I0, FALSE))
        p0 <- NCOL(x)
    }
    x.names <- colnames(x)
    
    if (is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))
		family <- family()
	if (is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}

    maxK <- min( maxK+length(keep), n0-1, p0 )
    rEBIC <- max( 0.0, 1.0-log(n0)/(2.0*log(p0)) )
    fitFun <- brglm2::brglmFit # ifelse(biasReduction, brglm2::brglmFit, glm.fit)

    idx <- keep
    egg <- fitFun(x=cbind(intercept=1.0, x[, idx]), y=y, family=family, ...)
    class(egg) <- 'glm'
    # stopifnot(egg$converged)
    # ebic <- -2*logLik(egg) + attr(logLik(egg), 'df') * log(attr(logLik(egg), 'nobs'))
    ebic <- BIC(egg)
    if (verbose) 
        message("Initial model (EBIC=", round(ebic, 4), ") ...")
    ebicCount <- 0
    while (TRUE) {
        D0 <- egg$family$mu.eta(egg$linear.predictor)
        S0 <- egg$family$variance(egg$fitted.values)
        pps <- drop( crossprod(x, D0/S0*(y-fitted(egg))) )
        pps.sd <- drop( sqrt(crossprod(x*x, D0^2/S0)) )

        m0 <- names( which.max( abs(pps)/pps.sd ) )
        tmp.idx <- unique( c(idx, m0) )
        # stopifnot( !(m0 %in% idx) )
        if (length(tmp.idx) > maxK) 
            break

        egg <- fitFun(x=cbind(intercept=1.0, x[, tmp.idx]), y=y, family=family, ...)
        class(egg) <- 'glm'
        if (!egg$converged) 
            break
        ebic.penalty <- 2*rEBIC*lchoose(p0-length(keep), length(tmp.idx)-length(keep))
        # tmp.ebic <- -2*logLik(egg) + 
        #   attr(logLik(egg), 'df') * log(attr(logLik(egg), 'nobs')) + 
        #   ebic.penalty
        tmp.ebic <- BIC(egg) + ebic.penalty

        if (tmp.ebic >= ebic)
            ebicCount <- ebicCount + 1
        if (ebicCount >= ebicFlag)
            break
        if (verbose)
            message("Add ", m0, " (EBIC=", round(tmp.ebic, 4), ") ...")
            
        ebic <- tmp.ebic
        idx  <- tmp.idx
    }

    # names(idx) <- colnames(x)[idx]
    # return(idx)
    as.character(idx)
}

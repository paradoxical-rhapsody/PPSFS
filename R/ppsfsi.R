#' @rdname PPS-method
#' @description `ppsfsi`: PPSFS for **interaction effects**.
#' @examples
#' ## ***************************************************
#' ## Identify interaction effects
#' ## ***************************************************
#' set.seed(2022)
#' n <- 300
#' p <- 150
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop( cbind(x[, 1:3], x[, 4:6]*x[, 7:9]) %*% runif(6, 1.0, 1.5) )
#' y <- eta + rnorm(n, sd=sd(eta)/5)
#' print( group <- ppsfsi(x, y, 'gaussian', verbose=TRUE) )
#' print( A <- ppsfs(x, y, "gaussian", I0=group, verbose=TRUE) )
#' 
#' print( A <- ppsfs(x, y, "gaussian", keep=c(1, "5:8"), 
#'                   I0=group, verbose=TRUE) )
#' 
#' @export
ppsfsi <- function (x, y, family, 
                    keep=NULL, 
                    ..., 
                    ebicFlag=1, 
                    maxK=min( NROW(x)-1, choose(NCOL(x), 2) ), 
                    verbose=FALSE) {
	n0 <- NROW(x)
    p0 <- choose(NCOL(x), 2) # NCOL(x)*(NCOL(x)-1)/2
    stopifnot(length(y) == n0)

    y <- na.fail(y)
    x <- na.fail(x)

    if (is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))
		family <- family()
	if (is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}

    rEBIC <- max( 0.0, 1.0-log(n0)/(2.0*log(p0)) )
    fitFun <- brglm2::brglmFit # ifelse (biasReduction, brglm2::brglmFit, glm.fit)
    
    group <- keep
    egg <- fitFun(x=cbind(intercept=1.0, get_inter_features(group, TRUE)), 
                y=y, family=family, ...)
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
        r0 <- D0 / S0 * (y - fitted(egg))
        W0 <- D0^2 / S0
        pps <- cxx_inter_pps(W0, r0, x)
        
        g0 <- paste(arrayInd(which.max(pps), dim(pps)), collapse=":")
        tmp.group <- unique(c(group, g0))
        # stopifnot( !(g0 %in% group) )
        if (length(tmp.group) > maxK) 
            break

        egg <- fitFun(x=cbind(intercept=1.0, get_inter_features(tmp.group, TRUE)), 
                    y=y, family=family, ...)
		class(egg) <- 'glm'
        if (!egg$converged) 
            break
        ebic.penalty <- 2*rEBIC*lchoose(p0-length(keep), length(tmp.group)-length(keep))
        # tmp.ebic <- -2*logLik(egg) + 
        #   attr(logLik(egg), 'df') * log(attr(logLik(egg), 'nobs')) + 
        #   ebic.penalty
        tmp.ebic <- BIC(egg) + ebic.penalty

        if (tmp.ebic >= ebic) 
            ebicCount <- ebicCount + 1
        if (ebicCount >= ebicFlag) 
            break
        if (verbose) 
            message("Add ", g0, " (EBIC=", round(tmp.ebic, 4), ") ...")

        ebic  <- tmp.ebic
        group <- tmp.group
    }

    group
}

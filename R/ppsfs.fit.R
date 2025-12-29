#' @title Partial Profile Score Feature Selection for GLMs
#' @description
#' This is a simplified but more stable implementation for main-effects selection
#' compared to [ppsfs].
#' It removes the standaridization, but adds the weights in the partial profile score.
#' 
#' @param x Matrix.
#' @param y Vector.
#' @param family See [glm] and [family].
#' @param fitFun Fitting function, by default [glm.fit].
#' @param ... Additional parameters for [glm.fit].
#' @param keep Initial index set of features that are included in model fitting.
#' @param maxK Maximum number of identified features.
#' @param verbose Print the feature path?
#' 
#' @details
#' That `ppsfs(x, y, family="gaussian")` is an implementation to
#' sequential lasso method proposed by Luo and Chen(2014, \doi{10/f6kfr6}).
#' 
#' @examples
#' set.seed(2025)
#' n <- 300
#' p <- 1000
#' x <- matrix(rnorm(n*p), n)
#' eta <- drop( x[, 1:3] %*% runif(3, 1.0, 1.5) )
#' y <- eta + rnorm(n, sd=sd(eta)/5)
#' print( A <- ppsfs.fit(x, y, 'gaussian', verbose=TRUE) )
#' 
#' @references
#' Z. Xu, S. Luo and Z. Chen (2022). Partial profile score feature selection in
#' high-dimensional generalized linear interaction models.
#' Statistics and Its Interface. \doi{10.4310/21-SII706}
#' 
#' @export
ppsfs.fit <- function (x, y, family, fitFun = glm.fit, ...,
                   keep = NULL, maxK = NULL, verbose = FALSE) {
    y  <- na.fail(y)
	x  <- na.fail(as.matrix(x))

    n0 <- NROW(x)
    p0 <- NCOL(x)
    stopifnot(length(y) == n0)

    p.keep <- length(keep)
    ebic.r <- max( 0.0, 1.0 - log(n0) / (2.0*log(p0)) )
    
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    if (!is.null(keep))
        stopifnot( all(is.numeric(keep)) )

    if (is.null(colnames(x)))
        colnames(x) <- as.character(seq_len(p0))
    x.names <- colnames(x)

    if (!is.null(maxK)) {
        stopifnot(maxK < NCOL(x))
        maxK <- min( maxK + length(keep), n0 - 1, p0 - 1 )
    }

    showiter <- function(verbose, k=length(idx), level=obj.level)
        if (verbose) message(sprintf("Step %i: model level = %.3f", k, level))

    scoreFun <- function(obj) {
        w0 <- obj[["prior.weights"]]
        D0 <- obj$family$mu.eta(obj$linear.predictor)
        S0 <- obj$family$variance(obj$fitted.values)
        return(drop( crossprod(x, w0 * D0 / S0 * residuals(obj, type="response"))) )
    }

    EBIC <- function(obj) {
        dof <- attr(logLik(obj), "df")
        BIC(obj) + 2.0 * ebic.r * lchoose(p0 - p.keep, dof - p.keep)
    }

    idx <- keep
    obj <- fitFun(x = x[, idx, drop = FALSE], y = y, family = family, ...)
    class(obj) <- "glm"
    obj.level <- EBIC(obj)
    while (TRUE) {
        showiter(verbose)

        m0 <- which.max(abs(scoreFun(obj)))
        idx.tmp <- unique( c(idx, m0) )
        obj <- fitFun(x = x[, idx.tmp, drop = FALSE], y = y, family = family, ...)
        class(obj) <- "glm"
        obj.level.tmp <- EBIC(obj)

        if (!is.null(maxK)) {
            if (length(idx.tmp) > maxK) break
        } else {
            if (obj.level.tmp > obj.level) {
                showiter(verbose, length(idx.tmp), obj.level.tmp)
                break
            }
        }

        obj.level <- obj.level.tmp
        idx  <- idx.tmp
    }

    return(idx)
}

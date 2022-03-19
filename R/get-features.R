#' @title Expand Interaction Features
#' @param group Index set of interaction effects.
#' @param include.main Whether to include the parent features of interaction effects?
#' 
#' @return Data matrix.
#' @noRd
get_inter_features <- function (group, include.main) { 
    stopifnot( exists('x') )
    # stopifnot( is.null(group) || all(grepl(':', group)) )

    if (is.null(group)) 
        return( get('x')[, numeric(0)] )

    idx <- sapply(strsplit(group, ':'), sort)
    mode(idx) <- 'numeric'
    j1 <- idx[1, ]
    j2 <- idx[2, ]

    j0 <- as.numeric( grep(':', group, value=TRUE, invert=TRUE) )
    if (include.main) 
        j0 <- sort(unique(c(j0, idx)))
    features <- cbind(get('x')[, j0], get('x')[, j1]*get('x')[, j2])
    colnames(features) <- c(j0, group)
    return(features)
}

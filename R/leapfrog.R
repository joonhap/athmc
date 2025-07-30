#' Leapfrog function
#'
#' carries out a specified number of leapfrog steps
#'
#' @param x a vector of initial position
#' @param p a vector of initial momentum
#' @param gd an R function that returns the gradient of log target density
#' @param jsize leapfrog step size
#' @param njumps number of leapfrog jumps
#' @param massInv mass-inverse matrix (currently, it is allowed to be only a scalar or a vector of the same length as v. When scalar, this means that the mass-inverse matrix is that constant times the identity matrix, and if `massInv` is a vector, then the mass-inverse matrix is the diagonal matrix having the diagonal entries given by `massInv`.)
#' @param temperature temperature value. The gradient of the potential function is divided by `temperature`.
#' @param bounds upper or lower bounds for each variable can be given. Should be a 2-by-length(x.init) matrix (first row: lower bound, second row: upper bound) or NULL (no bounds). If given, the constructed trajectory will bounce on the boundary, so that it stays within the specified bounds. Defaults to NULL (no bounds specified).
#' @return A concatenated vector `c(x,p)` of the final position and velocity
#' @export
lf <- function(x, p, gd, jsize, njumps, massInv=1, temperature=1, bounds=NULL) {
    ## to ensure time reversibility, the momentum reflection step is inserted between the two consecutive half-step updates of p (when `bounds` is specified).
    x.d <- length(x)
    if (is.null(bounds)) {
        bounds <- rbind(rep(-Inf, x.d), rep(Inf, x.d))
    }
    x <- x + .5*jsize*massInv*p
    j <- 0
    repeat {
        j <- j+1
        if (j >= njumps) break
        p <- p + .5*jsize*gd(x)/temperature
        p[which(x<bounds[1,])] <- abs(p[which(x<bounds[1,])]) # velocity reflection at the boundaries
        p[which(x>bounds[2,])] <- -abs(p[which(x>bounds[2,])])
        p <- p + .5*jsize*gd(x)/temperature
        x <- x + jsize*massInv*p
    }
    p <- p + .5*jsize*gd(x)/temperature
    p[which(x<bounds[1,])] <- abs(p[which(x<bounds[1,])]) # velocity reflection at the boundaries
    p[which(x>bounds[2,])] <- -abs(p[which(x>bounds[2,])])
    p <- p + .5*jsize*gd(x)/temperature
    x <- x + .5*jsize*massInv*p
    return(c(x,p))
}

## NOTE: A better approach could be as follows. Do the updates in the order p-x-x-p, and insert the momentum reflection step in between the two x-updates. This modification may circumvent the issue that the gradient of log target density might not be defined outside the interval constraints. Moreover, negating the momentum in between the two x-updates ensures that at the end of each full leapfrog step, all x-coordinates are within the boundaries. (NOW IMPLEMENTED BELOW AS lf_alternative())
lf_alternative <- function(x, p, gd, jsize, njumps, massInv=1, temperature=1, bounds=NULL) {
    ## x, p: initial position and momentum
    ## gd: gradient of log target density (R function)
    ## jsize: leapfrog step size, njumps: number of leapfrog jumps
    ## massInv: mass-inverse matrix (currently assumes it is a scalar or a vector of diagonal entries)
    ## temperature: temperature value alpha. The potential function U is replaced by (1/alpha)*U
    ## ... diagonal vector and defined as the vector of diagonal entries)
    ## bounds: each column gives the lower (first row) and the upper (second row) bounds for each variable.
    x.d <- length(x)
    if (is.null(bounds)) {
        bounds <- rbind(rep(-Inf, x.d), rep(Inf, x.d))
    }
    p <- p + .5*jsize*gd(x)/temperature
    j <- 0
    repeat {
        j <- j+1
        if (j >= njumps) break
        x <- x + .5*jsize*massInv*p
        p[which(x<bounds[1,])] <- abs(p[which(x<bounds[1,])]) # velocity reflection at the boundaries
        p[which(x>bounds[2,])] <- -abs(p[which(x>bounds[2,])])
        x <- x + .5*jsize*massInv*p
        p <- p + jsize*gd(x)/temperature
    }
    x <- x + .5*jsize*massInv*p
    p[which(x<bounds[1,])] <- abs(p[which(x<bounds[1,])]) # velocity reflection at the boundaries
    p[which(x>bounds[2,])] <- -abs(p[which(x>bounds[2,])])
    x <- x + .5*jsize*massInv*p
    p <- p + .5*jsize*gd(x)/temperature
    return(c(x,p))
}


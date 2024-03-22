#' Leapfrog function
#'
#' carries out a certain number of leapfrog steps
#'
#' @param x a vector of initial position
#' @param v a vector of initial velocity
#' @param gd an R function that returns the gradient of log target density
#' @param jsize leapfrog step size
#' @param njumps number of leapfrog jumps
#' @param massInv mass-inverse matrix (currently, it is allowed to be only a scalar or a vector of the same length as v. When scalar, this means that the mass-inverse matrix is that constant times the identity matrix, and if `massInv` is a vector, then the mass-inverse matrix is the diagonal matrix having the diagonal entries given by `massInv`.)
#' @return A concatenated vector `c(x,v)` of the final position and velocity
#' @export
lf <- function(x, v, gd, jsize, njumps, massInv=1) {
  x <- x + .5*v*jsize
  j <- 0
  repeat {
    j <- j+1
    if (j >= njumps) break
    v <- v + massInv*gd(x)*jsize
    x <- x + v*jsize
  }
  v <- v + massInv*gd(x)*jsize
  x <- x + .5*v*jsize
  return(c(x,v))
}

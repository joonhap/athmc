## leapfrog function
lf <- function(x, v, gd, jsize, njumps, massInv=1) {
  ## x, v: initial position and velocity
  ## gd: gradient of log target density (R function)
  ## jsize: leapfrog step size, njumps: number of leapfrog jumps
  ## massInv: mass-inverse matrix (currently assumes it is a ...
  ## ... diagonal vector and defined as the vector of diagonal entries)
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

#' Automatically-tuned, tempered Hamiltonian Monte Carlo
#'
#' runs tempered Hamiltonian Monte Carlo that facilitates sampling from high-dimensional, strongly multimodal distributions with the capability of automatic tuning. See Park (2024) for details.
#' 
#' @param x.init initial condition
#' @param logtarget log target density function (R function)
#' @param gradlt gradient of log target density (R function)
#' @param sumstat an R function that is applied to the vector `x` to create a summary statistic for each MCMC iteration for output. Here `x` is the current state of the Markov chain. The return value should be a vector. The default is the identity function. This argument is useful for avoiding high memory usage when `x` is high dimensional.
#' @param niter number of MCMC iterations. Note that the length of the constructed chain is `niter+1`, including the initial state.
#' @param massInv mass-inverse matrix (currently, it is allowed to be only a scalar or a vector of the same length as v. When scalar, this means that the mass-inverse matrix is that constant times the identity matrix, and if `massInv` is a vector, then the mass-inverse matrix is the diagonal matrix having the diagonal entries given by `massInv`.)
#' @param power starting value for the estimated polynomial degree of the (local) growth rate of the potential function. The leapfrog step size for each step is given by jsize*exp(4*eta/(power+2))
#' @param jsize starting value for the baseline leapfrog step size
#' @param maxEta starting value for the maximum of the eta sequence, which is one half times the log of the mass scaling factor (alpha)
#' @param lenEta starting value for the length of the eta sequence.
#' @param etaType either "piecewiselinear" or "sinusoidal". If "piecewiselinear", eta(k) = maxEta*2/lenEta*min(k,lenEta-k), and if "sinusoidal", eta(k) = maxEta/2*(1-cos(2*pi*k/lenEta)).
#' @param tune logical. should `jsize`, `maxEta`, `lenEta`, `power` be tuned? If so, the provided values for these arguments are used as the starting point of tuning.
#' @param MaxTuningIter the maximum number of tuning iterations before stopping further tuning, for each MCMC iteration
#' @param maxEta_tuning_method tuning method for maxEta, one of "rectangular", "ellipsoidal", "Uthreshold", or "none". "rectangular" requires the simulated path to reach outside a certain interval for every coordinate axis at least once. "ellipsoidal" requires the simulated path to have a point that is at least a certain (scaled) Euclidean distance away from a specified center point. "Uthreshold" requires the simulated trajectory to attain a potential energy level higher than the provided value. If "none", no tuning of `maxEta` is carried out.
#' @param search_scale a scalar or a vector of length `length(x.init)`, indicating the scale of search for isolated modes for each coordiate direction. If a scalar, every coordinate direction has the same search scale. For the "rectangular" tuning method, this gives the half-width of the interval for each coordinate direction with center at the corresponding entry in `center_point`. For "ellipsoidal" tuning method, this gives the coordinate-scale for the ellipsoid, scaled by sqrt(d) where d is the dimension of the target space, i.e., the criterion is distance^2 = sum((x-center_point)^2/search_scale^2) > d. Only used when `maxEta_tuning_method` is "rectangular" or "ellipsoidal".
#' @param center_point the center of search for isolated modes, used when `maxEta_tuning_method` is "rectangular" or "ellipsoidal". A vector of equal length as `x.init`. When missing, The default is set to `x.init`.
#' @param Uthreshold If `maxEta_tuning_method` is "Uthreshold", then this argument gives the target value that the potential function along the simulated trajectory should exceed at least once.
#' @param leapfrog optional parameter. If this argument is not given, the default leapfrog function will be used (see R/leapfrog.R). Otherwise, a custom R function for the leapfrog function can be used (see, for example, R/leapfrog_sensor.R).
#' @return A list consisting of several named entries. Note that `ATHMC` constructs a chain of length `niter+1`. The named entries are as follows.
#' \itemize{
#' \item{sumstat_MCMC: a matrix with `niter+1` rows, where each row gives the `sumstat` function applied to a state in the constructed chain.}
#' \item{lenEta_chrono: a vector of length `niter+1` recording the tuned values for the length of the eta (mass-scaling) schedule at the end of each iteration. The first value gives the initial, supplied value.}
#' \item{jsize_chrono: a vector of length `niter+1` recording the tuned values for the baseline leapfrog step size.}
#' \item{power_chrono: a vector of length `niter+1` recording the estimated polynomial degree of the potential function growth (gamma-hat).}
#' \item{maxEta_chrono: a vector of length `niter+1` recording the tuned maximum value of the eta sequence.}
#' \item{tuningIter_chrono: a vector of length `niter` recording the number of tuning cycles for each MCMC iteration.}
#' \item{Hinc_chrono: a vector of length `niter` recording the net increase in Hamiltonian over the trajectory constructed after tuning for that iteration.}
#' }
#' 
#' @references Park, J. (2024). Sampling from high-dimensional, multimodal distributions using automatically-tuned, tempered Hamiltonian Monte Carlo <https://doi.org/10.48550/arXiv.2111.06871>
#' @export
ATHMC <- function(x.init, logtarget, gradlt, sumstat=identity, massInv=1, niter, power=2, jsize=0.1, maxEta=1, lenEta=200, etaType="piecewiselinear", tune=TRUE, MaxTuningIter=70, maxEta_tuning_method="rectangular", search_scale, center_point, Uthreshold, leapfrog) {
    if (!missing(leapfrog)) {
        lf <- leapfrog # use custom leapfrog function 
    }
    x.d <- length(x.init) # space dimension
    if (tune && maxEta_tuning_method %in% c("rectangular", "ellipsoidal")) {
        if (missing(search_scale)) {
            stop("If `maxEta_tuning_method` is either 'rectangular' or 'ellipsoidal', `search_scale` should be provided.")
        }
        if (missing(center_point)) {
            center_point <- x.init
        }
        MAXETA_TUNE <- TRUE
    } else if (tune && maxEta_tuning_method=="Uthreshold") {
        if (missing(Uthreshold)) {
            stop("If `maxEta_tuning_method` is 'Uthreshold', then the `Uthreshold` arugment should be provided.")
        }
        MAXETA_TUNE <- TRUE
    } else if (tune && maxEta_tuning_method=="none") {
        MAXETA_TUNE <- FALSE
    } else if (tune) {
        stop("`maxEta_tuning_method` should be one of 'rectangular', 'ellipsoidal', 'Uthreshold', or 'none'.")
    } else { # if tune == FALSE
        MAXETA_TUNE <- FALSE
    }

    i <- 0 # mcmc iteration count
    x <- x.init # the current state of the Markov chain
    sumstat_length <- length(sumstat(x.init)) # length of the summary statistic
    sumstat_MCMC <- matrix(NA, niter+1, sumstat_length)
    sumstat_MCMC[1,] <- sumstat(x.init)
    lenEta_chrono <- numeric(niter+1); lenEta_chrono[1] <- lenEta # chronology of lenEta
    jsize_chrono <- numeric(niter+1); jsize_chrono[1] <- jsize
    power_chrono <- numeric(niter+1); power_chrono[1] <- power
    maxEta_chrono <- numeric(niter+1); maxEta_chrono[1] <- maxEta
    tuningIter_chrono <- numeric(niter) # chronology of number of tuning iterations
    Hinc_chrono <- numeric(niter) # chronology of increment in H over the (final) simulated trajectory

    repeat { ## start MCMC
        Lambda <- runif(1) # to accept/reject a proposed candidate
        i <- i+1; if (i > niter) { break }
        
        iter <- 0 # iteration number for adaptive tuning cycle
        END_ITER_NEXT_TIME <- FALSE

        if (MAXETA_TUNE) {
            maxEta <- ifelse(i==1, maxEta, max(maxEta-1,0.5)) # starting value for tuning maxEta. At each MCMC step, restart from the tuned value of maxEta minus 1.
            REACHED <- FALSE # indicates if the simulated path met the corresponding search criteria.
        } else {
            REACHED <- TRUE
        }

        v <- rnorm(x.d, 0, sqrt(massInv)) # this is actually v-tilde in the manuscript
        repeat{ # start iterative tuning
            iter <- iter + 1

            if (tune && maxEta_tuning_method == "rectangular") {
                REACHED_COORD <- rep(FALSE, x.d) # each entry becomes TRUE if the simulated trajectory cross the interval boundary for the corresponding coordinate.
            }

            if (etaType=="piecewiselinear") {
                eta <- maxEta*pmin(0:(2*lenEta), 2*lenEta-0:(2*lenEta))/lenEta # piecewise linear eta sequence, including both integer and half-integer k (even and odd indices, respectively)
            } else if (etaType=="sinusoidal") {
                eta <- maxEta*1/2*(1-cos(0:(2*lenEta)*2*pi/(2*lenEta))) # sinusoidal eta sequence, including both integer and half-integer k (even and odd indices, respectively) 
            } else {
                stop("etaType should be either 'piecewiselinear' or 'sinusoidal'.")
            }
            alpha <- exp(2*eta) # mass scaling factor
            simlen <- lenEta+1 # number of intermediate points along the trajectory, including the first and the last

            xv <- matrix(NA, simlen, 2*x.d)
            Utrace <- numeric(simlen)
            Ktrace <- numeric(simlen)
            alphatrace <- numeric(simlen)
            Utrace[1] <- U1 <- -logtarget(x)
            k <- 0#sample(c(0:4, lenEta-4:1), 1)#round(lenEta*.4) # initial mass scaling factor
            Ktrace[1] <- .5*sum(v*v/massInv)
            xv[1,] <- c(x,v)
            alphatrace[1] <- alpha[2*k+1]
            n <- 1; xold <- x; vold <- v

            ## simulate the trajectory for tilde-t, tilde-v, and x
            repeat{
                n <- n+1
                xvnew <- lf(xold, vold, gd=gradlt, jsize=jsize*(alpha[2*k+2])^(2/(power+2)),
                    njumps=1, massInv=massInv/alpha[2*k+2])
                xv[n,] <- xvnew
                xnew <- xvnew[1:x.d]; vnew <- xvnew[x.d+1:x.d]
                Utrace[n] <- Unew <- -logtarget(xnew)
                Ktrace[n] <- .5*sum(vnew*vnew/massInv)
                if (MAXETA_TUNE && !REACHED) { # if the search criterion is not yet met,
                    if (maxEta_tuning_method == "rectangular") {
                        REACHED_COORD <- pmax(REACHED_COORD, abs(xnew-center_point)>search_scale)
                        REACHED <- all(REACHED_COORD)
                    } else if (maxEta_tuning_method == "ellipsoidal") {
                        REACHED <- max(REACHED, sum(((xnew-center_point)/search_scale)^2)>x.d)
                    } else if (maxEta_tuning_method == "Uthreshold") {
                        REACHED <- max(REACHED, Unew > Uthreshold)
                    }
                }
                k <- k+1
                alphatrace[n] <- alpha[2*k+1]
                xold <- xnew; vold <- vnew
                if (n >= simlen) { break }
            }

            vbar <- apply(xv[,x.d+1:x.d], 2, function(vseq) {
                vseq*alphatrace^(1/(power+2))
            })
            Kbar <- Ktrace*alphatrace^(2/(power+2))

            ## assessment of the trajectory
            stabl_lim <- max(which(Kbar<1e30)) # approximate limit where the simulation remains stable
            if (stabl_lim < 3) { 
                stop("Numerical simulation is unstable at the third leapfrog step. Consider decreasing the starting value of jsize (leapfrog step size).")
            }
            sign_inc <- (Kbar[2:stabl_lim] > Kbar[1:(stabl_lim-1)]) # is the increment positive?
            cycles_len <- numeric(0) # the length of cycles
            cycles_amp <- numeric(0) # the amplitude of cycles
            for (n in 1:(stabl_lim-1)) {
                if (n==1) {
                    cumulator <- 0 # start counting the length of the current cycle
                    cycle_index <- 1 # the number of cycles counted at the moment
                    cycle_min <- min(Kbar[1], Kbar[2])
                    cycle_max <- max(Kbar[1], Kbar[2])
                    next
                }
                cumulator <- cumulator + 1
                if (sign_inc[n]) { # if the current increment is positive
                    if (sign_inc[n-1]) { # if the previous increment was also positive,
                        cycle_max <- max(cycle_max, Kbar[n+1]) # update the maximum value in the current cycle
                    } else { 
                        cycles_len[cycle_index] <- cumulator # record the current cycle length
                        cycles_amp[cycle_index] <- cycle_max - cycle_min # record the amplitude
                        cycle_index <- cycle_index + 1 # consider a new cycle
                        cumulator <- 0 # refresh the cumulator
                        cycle_min <- cycle_max <- Kbar[n+1] # refresh the min and max for the new cycle
                    }
                } else { # if the current increment is negative,
                    cycle_min <- min(cycle_min, Kbar[n+1])
                }
                if (n==stabl_lim-1) { # finally, add the remainder to the last cycle
                    cycles_len[cycle_index] <- cumulator 
                    cycles_amp[cycle_index] <- cycle_max - cycle_min 
                }
            }
            ncycles <- length(cycles_len) # number of cycles
            Hinc <- Utrace[length(Utrace)]+Ktrace[length(Ktrace)]-Utrace[1]-Ktrace[1]

            if (!tune) { break }
            if (END_ITER_NEXT_TIME) { break } # if tuning was ended in the previous round, exit loop here.

### ASSESSMENTS
            
            ## ASSESSMENT 1. Is the number of cycles adequate? Specifically, are there at least five completed cycles and no more than forty completed cycles in the trajectory?
            TOO_FEW_CYCLES <- (ncycles < 15)
            TOO_MANY_CYCLES <- (ncycles > 100)
            ETA_LEN_ADJUST_FACTOR <- 25/ncycles # adjust the length of the eta sequence by this factor, target=10
            ## ASSESSMENT 2. Is jsize (=epsilon) adequate? Rule of thumb: at least about 10 leapfrog steps per cycle of K(vbar). Count the number of consecutive increases in K(vbar) and the number of consecutive decreases that follow. The sum of the two should be about 10~15, and should be consistent throughout the trajectory.
            MEDIAN_CYCLE_LEN <- median(cycles_len[1:ncycles])
            EPS_TOO_LARGE <- (MEDIAN_CYCLE_LEN < 10)
            EPS_TOO_SMALL <- (MEDIAN_CYCLE_LEN > 100)
            EPS_ADJUST_FACTOR <- MEDIAN_CYCLE_LEN/20 # adjust jsize by this factor, target=20

            ## ASSESSMENT 3. Is `power' (=gamma-hat in the manuscript) adequate? Specifically, is the median amplitude for cycles in the middle third of the trajectory within [1/2, 2] times median amplitude of the cycles in first third?
            VBAR_FIRST_PART <- apply(vbar, 2, function(vbseq) {
                max(vbseq[1:ceiling(stabl_lim*1/8)]^2)
            }) 
            VBAR_SECOND_PART <- apply(vbar, 2, function(vbseq) {
                max(vbseq[ceiling(stabl_lim*3/8):ceiling(stabl_lim*4/8)]^2)
            }) 
            OSC_RATIO <- median(VBAR_SECOND_PART / VBAR_FIRST_PART)
            POWER_TOO_LARGE <-  !is.na(OSC_RATIO) && (log(OSC_RATIO) < -0.2)
            POWER_TOO_SMALL <- !is.na(OSC_RATIO) && log(OSC_RATIO) > 0.2
            EST_A_NEXT <- 2/(power+2) - 0.3*log(OSC_RATIO)/(eta[2*stabl_lim*3.5/8+1]-eta[2*stabl_lim*1/8+1])
            EST_POWER <- 2/EST_A_NEXT - 2

            if (abs(power-EST_POWER) < .3) { 
                power <- EST_POWER
            } else if (power < EST_POWER) { # if the difference between the current power and EST_POWER is more than 0.3, then change power by at most 0.3 (in the direction to EST_POWER) 
                power <- power + .3
            } else {
                power <- power - .3
            }

            ## if both eta sequence length and epsilon need change, then adjust it by the square roots of ETA_LEN_ADJUST_FACTOR and EPS_ADJUST_FACTOR
            lenEta <- min(round(lenEta * sqrt(ETA_LEN_ADJUST_FACTOR)), 10000)
            jsize <- jsize * sqrt(EPS_ADJUST_FACTOR) * runif(1,0.95,1.05)

            ## if the simulated path didn't meet the specified search criterion, increase maxEta by 0.4; if it did, decrease maxEta by 0.1.
            if (MAXETA_TUNE && !REACHED) {
                maxEta <- maxEta + 0.4
            } 

            if (!(POWER_TOO_LARGE|POWER_TOO_SMALL) && !(TOO_FEW_CYCLES|TOO_MANY_CYCLES) && !(EPS_TOO_SMALL|EPS_TOO_LARGE) && (!MAXETA_TUNE || REACHED) && !is.nan(Hinc)){ 
                END_ITER_NEXT_TIME <- TRUE # tuning is sufficient. simulate with the tuned parameter and end the loop
            }

            if (iter > MaxTuningIter) { break }

        } # END OF AUTOMATIC TUNING LOOP
        
        if (Hinc < -log(Lambda)) {
            x <- xnew; v <- vnew
        }
        sumstat_MCMC[i+1,] <- sumstat(x)
        lenEta_chrono[i+1] <- lenEta
        jsize_chrono[i+1] <- jsize    
        power_chrono[i+1] <- power
        maxEta_chrono[i+1] <- maxEta
        tuningIter_chrono[i] <- iter
        Hinc_chrono[i] <- Hinc

        if (i >= 5 && mean(tuningIter_chrono[i-(4:0)]) < 4) { # if the average number of tuning iterations in the past five MCMC iterations is less than 3,
            tune <- FALSE # freeze the tuned parameters
        }
    } # END OF MCMC LOOP

    return(list(sumstat_MCMC=sumstat_MCMC, lenEta_chrono=lenEta_chrono, jsize_chrono=jsize_chrono, power_chrono=power_chrono, maxEta_chrono=maxEta_chrono, tuningIter_chrono=tuningIter_chrono, Hinc_chrono=Hinc_chrono))
}

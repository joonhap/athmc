#' Automatically-tuned, tempered Hamiltonian Monte Carlo
#'
#' runs tempered Hamiltonian Monte Carlo that facilitates sampling from high-dimensional, strongly multimodal distributions with the capability of automatic tuning. See Park (2025) for details.
#' @param x.init initial condition
#' @param logtarget log target density function (R function)
#' @param gradlt gradient of log target density (R function)
#' @param sumstat an R function that is applied to the vector `x` to create a summary statistic for each MCMC iteration for output. Here `x` is the current state of the Markov chain. The return value should be a vector. The default is the identity function. This argument is useful for avoiding high memory usage when `x` is high dimensional.
#' @param niter number of MCMC iterations. Note that the length of the constructed chain is `niter+1`, including the initial state.
#' @param massInv inverse of the mass. a scalar, or a vector giving the diagonal entries of a diagonal matrix. A future upgrade of this package may support massInv that is a symmetric, positive definite matrix
#' @param power starting value for the estimated polynomial degree of the (local) growth rate of the potential function. The leapfrog step size for each step is given by jsize*exp(4*eta/(power+2))
#' @param power.fixed logical. If TRUE, the provided value of `power` is used without change, even if `tune`=TRUE.
#' @param jsize starting value for the baseline leapfrog step size
#' @param jsize.fixed: logical. If TRUE (default), jsize is not adaptively tuned when `tune`=TRUE. If FALSE, jsize is adaptively tuned when `tune`=TRUE.
#' @param bounds upper or lower bounds for each variable can be given. Should be a 2-by-length(x.init) matrix (first row: lower bound, second row: upper bound) or NULL (no bounds).
#' @param pilot_jsize_tune logical. If `TRUE`, standard HMC will be run to tune `jsize` (leapfrog step size), where maxEta is fixed (default value is 0). Useful for selecting the reference leapfrog step size (denoted epsilon-bar in Park (2025)) before running tempered HMC. `jsize` is adjusted so that the average acceptance rate is close to `pilot_target_avg_accep`. If `pilot_jsize_tune` is TRUE, `tune` is forced to be FALSE.
#' @param pilot_target_avg_accep target acceptance probability when tuning jsize in a pliot run with maxEta=0
#' @param maxEta starting value for the maximum of the eta sequence, which is one half times the log of the mass scaling factor (alpha)
#'     ## maxEta: starting value for the maximum of the eta sequence, which is one half times the log of the mass scaling factor (alpha)
#' @param maxEta.fixed logical. If TRUE, the provided value of `maxEta` is used without change, even if `tune`=TRUE.
#' @param nsteps number of leapfrog steps if tune=FALSE or the starting number of leapfrog steps (while adaptive tuning) if tune=TRUE.
#' @param nsteps.fixed logical. If TRUE, the provided value of `nsteps` is used without change, even if `tune`=TRUE.
#' @param target_avg_accep target acceptance probability when tuning the rate at which eta changes. Used for adaptively tuning `nsteps`. `nsteps` is determined as the integer ceiling of 2*maxEta/(etaRate*jsize).
#' @param etaType either "piecewiselinear" or "sinusoidal". If "piecewiselinear", eta(k) = maxEta*2/lenEta*min(k,lenEta-k), and if "sinusoidal", eta(k) = maxEta/2*(1-cos(2*pi*k/lenEta)).
#' @param tune logical. Should at least one of `maxEta`, `nsteps`, or `power` be tuned? If so, the provided values for these arguments are used as the starting point of tuning.
#' @param MaxTuningIter the maximum number of tuning iterations before stopping further tuning, for each MCMC iteration
#' @param maxEta_tuning_method tuning method for maxEta, one of "rectangular", "ellipsoidal", "Uthreshold", or "none". "rectangular" requires that for at least half of the coordinates, the simulated trajectory reaches both sides of the complement of the interval (`search_center` +- `search_scale`). "ellipsoidal" requires the simulated path to have a point that is at least a certain (scaled) Euclidean distance away from a specified center point. "Uthreshold" requires the simulated trajectory to attain a potential energy level higher than the provided value. If "none", no tuning of `maxEta` is carried out.
#' @param search_center the center of search for isolated modes, used when `maxEta_tuning_method` is "rectangular" or "ellipsoidal". A vector of equal length as `x.init`. When missing, The default is set to `x.init`.
#' @param search_scale a scalar or a vector of length `length(x.init)`, indicating the scale of search for isolated modes for each coordiate direction. If a scalar, every coordinate direction has the same search scale. For the "rectangular" tuning method, this gives the half-width of the interval for each coordinate direction with center at the corresponding entry in `center_point`. For "ellipsoidal" tuning method, this gives the coordinate-scale for the ellipsoid, scaled by sqrt(d) where d is the dimension of the target space, i.e., the criterion is distance^2 = sum((x-center_point)^2/search_scale^2) > d. Only used when `maxEta_tuning_method` is "rectangular" or "ellipsoidal".
#' @param Uthreshold If `maxEta_tuning_method` is "Uthreshold", then this argument gives the target value that the potential function along the simulated trajectory should exceed at least once.
#' @param stay_at_maxEta The proportion of time (in t-bar, which is approximately proportional to the number of leapfrog steps) at which the temperature (eta) stays at the maximum value. Default value is 0.
#' @param pbar.plot The coordinates of the sequence of pbar (scaled momentum) to be plotted for each trajectory constructed. If NULL, no plot is generated.
#' @param verbose logical. If TRUE, diagnostic messages will be printed for the ongoing tuning process.
#' @return A list consisting of several named entries. Note that `ATHMC` constructs a chain of length `niter+1`. The named entries are as follows.
#' \itemize{
#' \item{sumstat: a matrix with `niter+1` rows, where each row gives the `sumstat` function applied to a state in the constructed chain.}
#' \item{jsize_chrono: If `pilot_jsize_tune` is TRUE, `jsize_chrono` gives the sequence of `jsize` (leapfrog step size) values obtained during the tuning process. If 'pilot_jsize_tune` is FALSE, 'jsize_chrono` returns the scalar value of `jsize`.}
#' \item{nsteps_chrono: If `nsteps` is tuned, a vector of length `niter+1` recording the tuned values for the length of the eta (temperature) schedule at the end of each iteration. The first value gives the initial, supplied value. If nsteps is not tuned, return the supplied scalar value of `nsteps`.}
#' \item{power_chrono: If `power` is tuned, a vector of length `niter+1` recording the estimated polynomial degree of the potential function growth (gamma-hat). Otherwise, the provided value of `power' (scalar).}
#' \item{maxEta_chrono: If `maxEta` is tuned, a vector of length `niter+1` recording the tuned maximum value of the eta sequence. Otherwise, the scalar value of the provided `maxEta`.}
#' \item{Hinc_chrono: a vector of length `niter` recording the net increase in Hamiltonian over the trajectory.}
#' }
#'
#' @references Park, J. (2025). Sampling from high-dimensional, multimodal distributions using automatically-tuned, tempered Hamiltonian Monte Carlo <https://doi.org/10.48550/arXiv.2111.06871>
#' @export
ATHMC <- function(x.init, logtarget, gradlt, sumstat=identity, massInv=1, niter, power=2, power.fixed=FALSE, jsize=0.1, jsize.fixed=TRUE, bounds=NULL, pilot_jsize_tune=FALSE, pilot_target_avg_accep=0.9, maxEta=0, maxEta.fixed=FALSE, nsteps=50, nsteps.fixed=FALSE, target_avg_accep=0.4, etaType="piecewiselinear", tune=FALSE, maxEta_tuning_method="rectangular", search_center, search_scale, Uthreshold, stay_at_maxEta=0, pbar.plot=NULL, verbose=FALSE) {
    x.d <- length(x.init)

    if (pilot_jsize_tune) {
        tune <- FALSE
    }

    if (tune && !maxEta.fixed && maxEta_tuning_method %in% c("rectangular", "ellipsoidal")) {
        if (missing(search_scale)) {
            stop("If `maxEta_tuning_method` is either 'rectangular' or 'ellipsoidal', `search_scale` should be provided.")
        }
        if (missing(search_center)) {
            search_center <- x.init
        }
        MAXETA_TUNE <- TRUE
    } else if (tune && !maxEta.fixed && maxEta_tuning_method=="Uthreshold") {
        if (missing(Uthreshold)) {
            stop("If `maxEta_tuning_method` is 'Uthreshold', then the `Uthreshold` arugment should be provided.")
        }
        MAXETA_TUNE <- TRUE
    } else if (tune && maxEta.fixed) {
        MAXETA_TUNE <- FALSE
    } else if (tune) { # if tune == TRUE and maxEta.fixed == FALSE and `maxEta_tuning_method` is NOT provided
        stop("`maxEta_tuning_method` should be one of 'rectangular', 'ellipsoidal', 'Uthreshold', or 'none'.")
    } else { # if tune == FALSE
        MAXETA_TUNE <- FALSE
    }

    i <- 0 # mcmc iteration count
    x <- x.init # the current state of the Markov chain
    sumstat_length <- length(sumstat(x.init)) # length of the summary statistic
    sumstat_MCMC <- matrix(NA, niter+1, sumstat_length)
    sumstat_MCMC[1,] <- sumstat(x.init)
    jsize_chrono <- numeric(niter+1); jsize_chrono[1] <- jsize # chronology of nsteps
    nsteps_chrono <- numeric(niter+1); nsteps_chrono[1] <- nsteps # chronology of nsteps
    power_chrono <- numeric(niter+1); power_chrono[1] <- power
    maxEta_chrono <- numeric(niter+1); maxEta_chrono[1] <- maxEta
    tuningIter_chrono <- numeric(niter) # chronology of number of tuning iterations
    Hinc_chrono <- numeric(niter) # chronology of increment in H over the (final) simulated trajectory

    etaRate <- 2*maxEta/nsteps/jsize # initial value of the rate of change of eta

    repeat { ## start MCMC
        Lambda <- runif(1) # to accept/reject a proposed candidate
        i <- i+1; if (i > niter) { break }

        iter <- 0 # iteration number for adaptive tuning cycle

        if (MAXETA_TUNE) {
            REACHED <- FALSE # indicates if the simulated path met the corresponding search criteria.
        } else {
            REACHED <- NA
        }

        p <- rnorm(x.d, 0, sqrt(1/massInv))
        jsize_perturb <- runif(1, 0.9, 1.1) # NEW in Revision1

        if (MAXETA_TUNE && maxEta_tuning_method == "rectangular") {
            REACHED_COORD_LO <- rep(FALSE, x.d) # each entry becomes TRUE if the simulated trajectory cross the lower boundary of the search criterion
            REACHED_COORD_UP <- rep(FALSE, x.d) # becomes TRUE if the trajectory cross the upper boundary
        }

        if(!nsteps.fixed && maxEta!=0) {
            nsteps <- ceiling(2*maxEta/etaRate/jsize) # number of leapfrog steps
        }

        if (etaType=="piecewiselinear") {
            eta <- maxEta*pmin(0:(2*nsteps), 2*nsteps-0:(2*nsteps))/nsteps # piecewise linear eta sequence, including both integer and half-integer k (even and odd indices, respectively)
        } else if (etaType=="sinusoidal") {
            eta <- maxEta*1/2*(1-cos(0:(2*nsteps)*2*pi/(2*nsteps))) # sinusoidal eta sequence, including both integer and half-integer k (even and odd indices, respectively)
        } else {
            stop("etaType should be either 'piecewiselinear' or 'sinusoidal'.")
        }
        nstaysteps <- round(nsteps*stay_at_maxEta/2) # the number of leapfrog steps that stays at maxEta
        eta <- c(eta[0:(nsteps-1)], rep(maxEta, 2*nstaysteps+1), eta[(nsteps+1):(2*nsteps)])
        alpha <- exp(2*eta) # mass scaling factor
        simlen <- nsteps+nstaysteps+1 # number of intermediate points along the trajectory, including the first and the last

        xp <- matrix(NA, simlen, 2*x.d)
        Utrace <- numeric(simlen)
        Ktrace <- numeric(simlen)
        alphatrace <- numeric(simlen)
        Utrace[1] <- -logtarget(x)
        k <- 0
        Ktrace[1] <- .5*sum(p*p*massInv)
        xp[1,] <- c(x,p)
        alphatrace[1] <- alpha[2*k+1]
        n <- 1; xold <- x; pold <- p

        ## simulate the trajectory
        repeat{
            n <- n+1
            xpnew <- lf(xold, pold, gd=gradlt, jsize=jsize*jsize_perturb*(alpha[2*k+2])^(2/(power+2)), njumps=1, massInv=massInv, temperature=alpha[2*k+2], bounds=bounds)
            xp[n,] <- xpnew
            xnew <- xpnew[1:x.d]; pnew <- xpnew[x.d+1:x.d]
            Utrace[n] <- Unew <- -logtarget(xnew)
            Ktrace[n] <- .5*sum(pnew*pnew*massInv)
            if (MAXETA_TUNE && !REACHED) { # if the search criterion is not yet met,
                if (maxEta_tuning_method == "rectangular") {
                    REACHED_COORD_LO <- pmax(REACHED_COORD_LO, xnew-search_center < -search_scale)
                    REACHED_COORD_UP <- pmax(REACHED_COORD_UP, xnew-search_center > search_scale)
                    REACHED_COORD <- REACHED_COORD_LO & REACHED_COORD_UP
                    REACHED <- (mean(REACHED_COORD) > 0.5)
                } else if (maxEta_tuning_method == "ellipsoidal") {
                    REACHED <- max(REACHED, sum(((xnew-search_center)/search_scale)^2)>x.d)
                } else if (maxEta_tuning_method == "Uthreshold") {
                    REACHED <- max(REACHED, Unew > Uthreshold)
                }
            }
            NUMERICAL_INSTABILITY <- FALSE
            if (MAXETA_TUNE && is.na(REACHED)) { # REACHED is NA due to numerical instability of the constructed path
                NUMERICAL_INSTABILITY <- TRUE
                break
            }
            k <- k+1
            alphatrace[n] <- alpha[2*k+1]
            xold <- xnew; pold <- pnew
            #if (runif(1) < .02) { newdirec <- rnorm(x.d); pold <- newdirec / sqrt(sum(newdirec^2)) * sqrt(sum(pold^2)) } # random scattering
            #if (runif(1) < .02) { pold <- pold*exp(runif(1,-.5,.5)) } # random scattering
            if (n >= simlen) { break }
        }

        pbar <- apply(xp[,x.d+1:x.d,drop=FALSE], 2, function(pseq) {
            pseq*alphatrace^(1/(power+2))
        })
        Kbar <- Ktrace*alphatrace^(2/(power+2))
        if (!is.null(pbar.plot)) {
            no_row <- ceiling(sqrt(length(pbar.plot)))
            no_col <- ceiling(length(pbar.plot)/no_row)
            par(mfrow=c(no_row, no_col), mar=c(2,2,0,0))
            for (coor in pbar.plot) {
                lplot(pbar[,coor]); points(pbar[,coor], col='green', pch='+')
            }
        }
        Hinc <- ifelse(NUMERICAL_INSTABILITY, NA, Utrace[length(Utrace)]+Ktrace[length(Ktrace)]-Utrace[1]-Ktrace[1])

        if (verbose) {
            if (pilot_jsize_tune) { cat("## PILOT RUN, TUNING jsize ") }
            cat("## i=",i, "maxEta=", maxEta, "nsteps=", nsteps, ", etaRate=", etaRate, ", jsize=", jsize, ", power=", power, ", Hinc=", Hinc, "REACHED=", REACHED, "U(x)=", ifelse(NUMERICAL_INSTABILITY, NaN, Utrace[length(Utrace)]), "\n")
        }

        ### If this is a pilot run for selecting the reference leapfrog step size (i.e., `pilot_jsize_tune`=TRUE), adjust `jsize`.
        if (pilot_jsize_tune) {
            if (is.na(Hinc)) {
                jsize <- exp(log(jsize) - (1/(1+i)^0.6))
            } else {
                jsize <- exp(log(jsize) + (1/(1+i)^0.6) * (min(1, exp(-Hinc)) - target_avg_accep))
            }
            jsize_chrono[i+1] <- jsize
        }

        ### TUNING
        if (tune) {
            ## ADJUSTMENT 1. Is eta varying at an appropriate rate?
            if (!nsteps.fixed) {
                if (is.na(Hinc) || NUMERICAL_INSTABILITY) {
                    etaRate <- exp(log(etaRate) - (1/(1+i)^0.6))
                } else {
                    etaRate <- exp(log(etaRate) - (1/(1+i)^0.6) * (target_avg_accep - min(1, exp(-Hinc))) )
                }
            }

            ## (Optional.) ADJUSTMENT. Is jsize small enough?
            if (!jsize.fixed) {
                if (is.na(Hinc) || NUMERICAL_INSTABILITY) {
                    jsize <- exp(log(jsize) - (0.3/(1+i)^0.6))
                } else {
                    jsize <- exp(log(jsize) - (0.3/(1+i)^0.6) * (target_avg_accep - min(1, exp(-Hinc))) )
                }
            }

            ## ADJUSTEMTN 2. Is `power' (=gamma-hat in the manuscript) adequate? Specifically, is the median amplitude for cycles in the middle third of the trajectory within [1/2, 2] times median amplitude of the cycles in first third?
            if (!power.fixed) { # if power is specified, skip tuning of power
                stabl_lim <- max(which(Kbar<1e30)) # approximate limit where the simulation remains stable
                if (stabl_lim < 3) {
                    stop("Numerical simulation is unstable at the third leapfrog step. Consider decreasing the starting value of jsize (leapfrog step size).")
                }
                PBAR_FIRST_PART <- apply(pbar, 2, function(pbseq) {
                    max(pbseq[1:ceiling(stabl_lim*1/8)]^2)
                })
                PBAR_SECOND_PART <- apply(pbar, 2, function(pbseq) {
                    max(pbseq[ceiling(stabl_lim*3/8):ceiling(stabl_lim*4/8)]^2)
                })
                KBAR_FIRST_PART <- max(Kbar[1:ceiling(stabl_lim*1/8)])
                KBAR_SECOND_PART <- max(Kbar[ceiling(stabl_lim*3/8):ceiling(stabl_lim*4/8)])
                OSC_RATIO <- median(PBAR_SECOND_PART / PBAR_FIRST_PART)
                ##OSC_RATIO <- KBAR_FIRST_PART / KBAR_SECOND_PART
                POWER_TOO_LARGE <-  !is.na(OSC_RATIO) && log(OSC_RATIO) < -0.2
                POWER_TOO_SMALL <- !is.na(OSC_RATIO) && log(OSC_RATIO) > 0.2
                EST_A_NEXT <- 2/(power+2) - 0.3*log(OSC_RATIO)/(eta[2*stabl_lim*3.5/8+1]-eta[2*stabl_lim*0.5/8+1]) # correction: the second index for eta was corrected from 2*stabl_lim*1/8+1 to 2*stabl_lim*0.5/8+1
                EST_POWER <- 2/EST_A_NEXT - 2
                if (verbose) { cat("## Tuning power (gamma-hat)"," power:", power,"EST_POWER:",EST_POWER,"EST_A_NEXT:",EST_A_NEXT,"OSC_RATIO:",OSC_RATIO,"\n") }
                ## Adjust power
                if (abs(power-EST_POWER) < .3) {
                    power <- EST_POWER
                } else if (power < EST_POWER) { # if the difference between the current power and EST_POWER is more than 0.3, then change power by at most 0.3 (in the direction to EST_POWER)
                    power <- power + (.2 + runif(1,0,0.1)) # NOTE: The change in power was randomized to avoid endlessly alternating power in Revision1.
                } else {
                    power <- power - (.2 - runif(1,0,0.1)) # NOTE: The change in power was randomized to avoid endlessly alternating power in Revision1.
                }
            }

            ## ADJUSTEMENT 3. If the simulated path didn't meet the specified search criterion, increase maxEta by 0.4; if it did, decrease maxEta by 0.1.
            if (MAXETA_TUNE && !is.na(REACHED)) {
                if (!REACHED) {
                    maxEta <- maxEta + 1/(i+1)^0.6 * 2
                } else {
                    maxEta <- max(maxEta - 1/(i+1)^0.6 * 1, 0)
                }
            }
        }

        if (!NUMERICAL_INSTABILITY && !is.na(Hinc) && Hinc < -log(Lambda)) {
            x <- xnew; p <- pnew
        }
        sumstat_MCMC[i+1,] <- sumstat(x)
        nsteps_chrono[i+1] <- nsteps
        power_chrono[i+1] <- power
        maxEta_chrono[i+1] <- maxEta
        Hinc_chrono[i] <- Hinc

    } # END OF MCMC LOOP

    returnval <- list(sumstat=sumstat_MCMC)

    if (pilot_jsize_tune) {
        returnval$jsize_chrono <- jsize_chrono
        returnval$nsteps_chrono <- nsteps
        returnval$power_chrono <- power
        returnval$maxEta_chrono <- maxEta
    }

    if (tune) {
        returnval$jsize_chrono <- jsize

        if (!nsteps.fixed) {
            returnval$nsteps_chrono <- nsteps_chrono
        } else {
            returnval$nsteps_chrono <- nsteps
        }

        if (!power.fixed) {
            returnval$power_chrono <- power_chrono
        } else {
            returnval$power_chrono <- power
        }

        if (!maxEta.fixed) {
            returnval$maxEta_chrono <- maxEta_chrono
        } else {
            returnval$maxEta_chrono <- maxEta
        }
    }

    returnval$Hinc_chrono <- Hinc_chrono

    return(returnval)
}

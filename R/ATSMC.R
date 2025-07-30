#' Tempered sequential Monte Carlo (TSMC) with adaptive tuning
#'
#' Tempered sequential Monte Carlo is also known as \emph{anealed importance sampling}. The adaptive tuning procedure mostly follow the methods proposed by Buchholz et al. (2021).
#' @param logtarget log target density function (R function)
#' @param gradlt gradient of log target density (R function)
#' @param sumstat an R function that is applied to the vector `x` to create a summary statistic for each MCMC iteration for output. Here `x` is the current state of the Markov chain. The return value should be a vector. The default is the identity function. This argument is useful for avoiding high memory usage when `x` is high dimensional.
#' @param base_gen An R function that generates a single draw from the base distribution. This function does not take any argument.
#' @param base_lpdf The log probability density function of the base distribution, p_0(x). Intermediate tempered densities are given by p_0(x)^(1-beta) * p(x)^beta, where p(x) is the target density, and beta is the inverse temperature.
#' @param base_glpdf The gradient of the log probability density function of the base distribution, log(p_0(x)) (see the description for `base_lpdf').
#' @param nparticle number of particles to be used.
#' @param adaptune logical. Should `ntemp' and `tempschedule' be tuned adaptively? If TRUE, both parameters will be tuned such that the effective sample size (ESS) for the importance sampling step for adjacent temperature levels is close to the specified `targetESSratio' times `nparticle'.
#' @param targetESSratio The target ratio of the effective sample size to number of particles for adaptively tuning `invtemp_schedule'.
#' @param invtemp_schedule `invtemp_schedule' gives a numeric vector of inverse temperature (beta) levels. It is expected that the vector forms an increasing sequence, where the first entry equals 0 and the last entry equals 1. If adaptune=TRUE, this argument will be ignored.
#' @param ntemp The number of temperature levels, including the base (beta=0) and the target (beta=1). If adaptune=FALSE, `invtemp_schedule' will be a uniform sequence between 0 and 1 of length `ntemp'. If 'invtemp_schedule' is provided, 'ntemp' will be ignored.
#' @param lfstepsize Starting value for tuning the leapfrog step size for the HMC kernel used in between importance sampling steps to perturb parameters (and replenish particle diversity). The leapfrog step size is tuned so that the acceptance rate is close to 0.9.
#' @param n_lfsteps Base number of leapfrog steps for the HMC kernel. Repeated until the correlation is close to zero (see Buchholz, Chopin, Jacob (2021) for details.)
#' @param MCMCsteps Number of HMC perturbations carried out in sequence for each particle at each temperature level
#' @param bounds the bounds for the variables (used for constructing HMC trajectories). See '?leapfrog' for more details.
#' @param return_all: logical. If TRUE, return (the summary statistics of) all particles at all temperature levels
#' @param verbose: logical. If TRUE, diagnostic messages are printed.
#' @param parallel: logical. Is particle rejuvenation step carried out in parallel using multiple cores? Defaults to TRUE.
#' @param ncore: optional. The number of computing cores to use for parallel particle rejuvenation. Automatically capped at the number of available cores minus one.
#' @return A return list. If `return_all` is FALSE, the list consists of:
#' \itemize{
#' \item{particles: The ensemble of particles at the final temperature level. These particles provide a Monte Carlo approximation to the target distribution. Note that the return values are the values of the function `sumstat` applied to the particles. Thus `particles` do not return the actual particles unless `sumstat` is the identity function (default).}
#' \item{weights: The normalized weights for the particles at the final step. Note: The `weights` vector does not correspond to the `particles` vector, in the sense that `particles` give the values after importance resampling with `weights`. Thus `weights` are provided mainly for diagnostic purposes.}
#' \item{invtemp_schedule: The sequence of inverse temperature levels. If `adaptune` is TRUE, the tuned sequence of values will be returned. Otherwise, the provided sequence will be returned.}
#' }
#' If `return_all` is TRUE, the list additionally consists of
#' \itemize{
#' \item{particles_book: A list consisting of particle ensembles _before_ the importance resampling steps at all temperature levels. The first list entry is a collection of iid draw from the base distribution (`base_gen`).}
#' \item{weights_book: A list of normalized weights at all temperature levels.}
#' }
#'
#' @references
#' \itemize{
#' \item{Buchholz, A., Chopin, N., and Jacob, P. E. (2021). Adaptive tuning of Hamiltonian Monte Carlo within sequential Monte Carlo. Bayesian Analysis, 16(3):745–771}
#' \item{Neal, R. M. (2001). Annealed importance sampling. Statistics and Computing, 11:125–139.}
#' \item{Park, J. (2025). Sampling from high-dimensional, multimodal distributions using automatically-tuned, tempered Hamiltonian Monte Carlo <https://doi.org/10.48550/arXiv.2111.06871>}
#' }
#'
#' @export
ATSMC <- function(logtarget, gradlt, sumstat=identity, base_gen, base_lpdf, base_glpdf, nparticle, adaptune=TRUE, targetESSratio=0.5, invtemp_schedule, ntemp=15, lfstepsize=0.1, n_lfsteps=10, MCMCsteps=5, bounds=NULL, return_all=FALSE, verbose=FALSE, parallel=FALSE, ncore) {

    if (!adaptune) {
        if (!missing(invtemp_schedule)) {
            ntemp <- length(invtemp_schedule)
            if (invtemp_schedule[1]!=0 || invtemp_schedule[ntemp]!=1) {
                warning("The first and last entries of invtemp_schedule should be 0 and 1, respectively. Proceeding as if they are.")
            }
        } else { ## if invtemp_schedule has not been provided
            invtemp_schedule <- seq(from=0, to=1, length.out=ntemp)
        }
    } else { # if adaptune=TRUE
        invtemp_schedule <- 0
        temp_saturated <- FALSE # have we found all necessary intermediate temperature levels?
    }

    particles <- rbind(replicate(nparticle, base_gen())) # draw particles from the base distribution, each column is a draw.
    if (return_all) {
        particles_book <- list() # empty list to keep the particles at all temperature levels
        particles_book[[1]] <- rbind(apply(particles, 2, sumstat))
        weights_book <- list() # empty list to keep the particle weights at all temperature levels
    }
    x.d <- dim(particles)[1] # dimension of each particle
    itlvl <- 1 # inverse temperature level (index)
    jsize <- lfstepsize # initialize the leapfrog step size
    repeat {
        itlvl <- itlvl + 1
        evalESS <- function(beta) {
            logweights <- (beta - invtemp_schedule[itlvl-1]) * (apply(particles, 2, logtarget) - apply(particles, 2, base_lpdf))
            adj_lw <- logweights - max(logweights)
            ESS <- sum(exp(adj_lw))^2 / sum(exp(2*adj_lw)) # effective sample size
            return(ESS)
        }
        if (adaptune) {
            ## find next beta (invtemp) that gives ESS = targetESSratio*nparticle
            beta_gap_ub <- ifelse(itlvl==2, 0.2, invtemp_schedule[itlvl-1] - invtemp_schedule[itlvl-2]) # upper bound for gap in beta
            repeat { # increase the upper bound for temperature gap until the ratio ESS/nparticle falls below targetESSratio
                if (evalESS(invtemp_schedule[itlvl-1]+beta_gap_ub) > targetESSratio*nparticle) {
                    beta_gap_ub <- beta_gap_ub * 2
                } else {
                    break
                }
            }
            beta_lower <- invtemp_schedule[itlvl-1] # lower bound for the beta search
            beta_upper <- beta_lower + beta_gap_ub
            repeat { # bisection method
                beta_mid <- (beta_upper + beta_lower) / 2
                ESS_beta_mid <- evalESS(beta_mid)
                if (ESS_beta_mid > (0.95*targetESSratio+0.05)*nparticle) {
                    beta_lower <- beta_mid
                } else if (ESS_beta_mid < (1.05*targetESSratio-0.05)*nparticle) {
                    beta_upper <- beta_mid
                } else { # if the temperature at temp_mid is between (targetESSratio +- .05*(1-targetESSratio)) * nparticle, stop.
                    break
                }
            }
            if (beta_mid > 1) { # no more temperature level is needed
                ntemp <- itlvl
                temp_saturated <- TRUE
                invtemp_schedule[itlvl] <- 1
            } else {
                invtemp_schedule[itlvl] <- beta_mid
            }
        }
        ## weighting and resampling
        if (verbose) { cat("## Inverse temperature: ", invtemp_schedule[itlvl], " (", itlvl, "-th level). Resampling the particles.\n") }
        logweights <- (invtemp_schedule[itlvl] - invtemp_schedule[itlvl-1]) * (apply(particles, 2, logtarget) - apply(particles, 2, base_lpdf))
        adj_lw <- logweights - max(logweights)
        ESS <- sum(exp(adj_lw))^2 / sum(exp(2*adj_lw)) # effective sample size
        norm_w <- exp(adj_lw) / sum(exp(adj_lw)) # normalized weights
        if (return_all) {
            weights_book[[itlvl-1]] <- norm_w
        }
        newparticles <- matrix(NA, x.d, nparticle)
        cw <- 0 # cumulative weight
        offset <- runif(1,0,1)/nparticle # offset for the thresholds for systematic resampling
        ii <- 0; jj <- 0
        repeat { # systematic resampling
            ii <- ii + 1 # increase the index for the new particle
            if (ii > nparticle) {
                break
            }
            if (offset + (ii-1)/nparticle > cw) {
                repeat {
                    jj <- jj + 1 # increase the index for the original particle
                    cw <- cw + norm_w[jj]
                    if (cw >= offset + (ii-1)/nparticle) { # if the cumulative weight exceeds the next threshold
                        newparticles[,ii] <- particles[,jj]
                        break
                    }
                }
            } else {
                newparticles[,ii] <- particles[,jj]
            }
        }
        ## at the final invtemp level, return the particles and the normalized weights
        if ((!adaptune && itlvl == ntemp) || (adaptune && temp_saturated)) {
            if (return_all) {
                return(list(particles=rbind(apply(newparticles, 2, sumstat)), weights=norm_w, invtemp_schedule=invtemp_schedule, particles_book=particles_book, weights_book=weights_book))
            } else {
                return(list(particles=rbind(apply(newparticles, 2, sumstat)), weights=norm_w, invtemp_schedule=invtemp_schedule))
            }
        }
        ## HMC particle rejuvenation
        logtempered <- function(x) {
            invtemp_schedule[itlvl] * logtarget(x) + (1-invtemp_schedule[itlvl]) * base_lpdf(x)
        }
        grad_logtempered <- function(x) {
            invtemp_schedule[itlvl] * gradlt(x) + (1-invtemp_schedule[itlvl]) * base_glpdf(x)
        }
        pilot_iter <- 30 # number of iterations for a pilot run to determine the leapfrog step size
        pilot <- ATHMC(newparticles[,sample(1:nparticle, 1)], logtempered, grad_logtempered, sumstat=identity, massInv=1, niter=pilot_iter, jsize=jsize, pilot_jsize_tune=TRUE, nsteps=n_lfsteps, target_avg_accep=0.9, bounds=bounds, verbose=verbose) # pilot run to find a mode and tune the reference leapfrog step size
        jsize <- pilot$jsize_chrono[pilot_iter+1]
        HMCmove <- function(x) { # apply HMC kernel
            p <- rnorm(x.d)
            H1 <- -logtempered(x) + 0.5*sum(p*p) # initial Hamiltonian
            xpnew <- lf(x, p, gd=grad_logtempered, jsize=jsize, njumps=n_lfsteps, bounds=bounds)
            xnew <- xpnew[1:x.d]
            pnew <- xpnew[x.d+1:x.d]
            Hinc <- -logtempered(xnew) + 0.5*sum(pnew*pnew) - H1 # increment in Hamiltonian
            if (-log(runif(1,0,1)) > Hinc) {
                x <- xnew # proposal is accepted
            }
            return(x)
        }
        if (verbose) { cat("## Inverse temperature: ", invtemp_schedule[itlvl], " (", itlvl, "-th level). Rejuvenating particles using HMC.\n") }
        particles <- newparticles
        MCMCiter <- 0
        if (parallel) {
            require(parallel)
            if (missing(ncore)) {
                ncore <- Inf
            }
            ncore <- min(ncore, detectCores()-1)
            repeat {
                MCMCiter <- MCMCiter + 1
                particles <- rbind(simplify2array(mclapply(1:nparticle, function(pid) { HMCmove(particles[,pid]) }, mc.cores=ncore)))
                if (MCMCiter >= MCMCsteps) {
                    break
                }
            }
        } else { # if not running in parallel
            repeat{
                MCMCiter <- MCMCiter + 1
                particles <- rbind(apply(particles, 2, HMCmove))
                if (MCMCiter >= MCMCsteps) {
                    break
                }
            }
        }
        if (return_all) {
            particles_book[[itlvl]] <- rbind(apply(particles, 2, sumstat))
        }
    }
}


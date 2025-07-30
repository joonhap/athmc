#' Parallel tempering using a HMC kernel with an adaptive tuning procedure
#'
#' Parallel tempering, in which each chain is constructed using HMC. Adopts the temperature level adjustment strategy proposed by Miasojedow, Moulines, and Vihola (2013), JCGS. Additionally, the number of parallel chains is adaptively tuned, such that the top level chain satisfies a specified search criterion.
#' @param x.init starting point for the Markov chain
#' @param target an R function that evaluates the log target density
#' @param gd.target an R function that evaluates the gradient of the log target density
#' @param maxiter the maximum number of MCMC iterations (the length of each parallel chain)
#' @param maxtemp starting value for the maximum log temperature
#' @param ntemp starting value for the number of parallel chains minus one (number of auxiliary chains besides the chain at temperature = 1
#' @param njumps number of leapfrog steps for HMC at each iteration of parallel chains
#' @param sumstat_fun the function to apply to compute the summary statistic to be returned. Defaults to the identity function.
#' @param ntemp_tuning_method method for tuning the number of parallel chains. Either "rectangular" (default) or "Uthreshold".
#' @param search_center the reference point for assessing whether the top level chain meets a search criterion. Used if `ntemp_tuning_method` is "rectangular".
#' @param search_scale the vector of length equal to the dimension of the space specifying the search scale for each coordinate. Used if `ntemp_tuning_method` is "rectangular".
#' @param Uthreshold the required value of the potential function U to be reached before stopping the tuning of ntemp. Used if `ntemp_tuning_method` is "Uthreshold".
#' @param span: a scalar value specifying how many past MCMC iterations at the top level should be considered when checking the search criterion.
#' @param verbose: logical. should detailed internal outputs be printed in buffer?
#' @return A list consisting of:
#' \itemize{
#' \item{jsizes: Tuned values of leapfrog step sizes for parallel chains.}
#' \item{tempseq: The (tuned) temperature levels.}
#' \item{acc_prob: The history of acceptance probabilities for the parallel chains.}
#' \item{sumstat: The history of parallel chains: the values are the evaluation of `sumstat_fun` on the chain states.}
#' \item{cumnjump: Cumulative number of leapfrog steps up to the current iteration (recorded for all parallel chains).}
#' }
#'
#' @references
#' \itemize{
#' \item{Park, J. (2025). Sampling from high-dimensional, multimodal distributions using automatically-tuned, tempered Hamiltonian Monte Carlo <https://doi.org/10.48550/arXiv.2111.06871>}
#' \item{Miasojedow, B., Moulines, E., and Vihola, M. (2013). An adaptive parallel tempering algorithm. Journal of Computational and Graphical Statistics, 22(3):649–664.}
#' }
#'
#' @export
apt <- function(x.init, target, gd.target, maxiter,  maxtemp, ntemp, njumps=50, sumstat_fun=identity, ntemp_tuning_method="rectangular", search_center, search_scale, Uthreshold, span=100, verbose=FALSE) {
    x.d <- length(x.init)
    tempseq <- exp(seq(from=0, to=maxtemp, length.out=ntemp)) # starting temperature sequence
    tempseq_history <- list()
    tempseq_history[[1]] <- tempseq
    acc_prob_history <- list()
    rho <- log(tempseq[2:ntemp] - tempseq[1:(ntemp-1)]) # initial value for the adaptive tuning parameter (see "An Adaptive Parallel Tempering Algorithm" (2013) by Błażej Miasojedow, Eric Moulines & Matti Vihola, JCGS)
    ## pilot run to select reference leapfrog step size
    jsizes <- numeric(ntemp)
    X <- matrix(NA, x.d, ntemp) # each column gives the current state of each chain
    for (k in 1:ntemp) {
        tem_target <- function(x) { target(x) / tempseq[k] }
        tem_gd.target <- function(x) { gd.target(x) / tempseq[k] }
        pilot_iter <- 50
        if (verbose) { cat("Tuning the leapfrog step size for the", k, "-th temperature", tempseq[k], ". \n") }
        pilot <- ATHMC_rev1(x.init, tem_target, tem_gd.target, sumstat=identity, massInv=1, niter=pilot_iter, jsize=0.5*tempseq[k]^0.5, pilot_jsize_tune=TRUE, nsteps=njumps, target_avg_accep=0.9, verbose=verbose) # pilot run to find a mode and tune the reference leapfrog step size
        X[,k] <- pilot$sumstat[pilot_iter+1,] # update the starting position
        jsizes[k] <- pilot$jsize_chrono[pilot_iter+1]
    }
    ## initialize parallel chains
    sumstat <- list() # summary statistic to report for each chain
    sumstat[[1]] <- sapply(1:ntemp, function(k) { sumstat_fun(X[,k]) })
    cumnjump <- rep(0, maxiter-1) # cumulative number of leapfrog steps (for all parallel chains)
    topchain <- matrix(NA, x.d, span) # the record of the past `span` MCMC iterations of the top level chain
    REACHED_TIMES <- 0 # the number of times the search criterion as reached. The search criterion is checked every ceiling(span/2) MCMC iterations. Tuning of ntemp is stopped if the search criterion is satisfied twice.
    NTEMP_TUNED <- FALSE # has `ntemp` been tuned? (no further tuning afterwards)
    ## adaptive parallel tempering
    for (i in 1:(maxiter-1)) {
        ## advance chains in parallel
        for (k in 1:ntemp) {
            x0 <- X[,k]
            p0 <- rnorm(x.d)
            H0 <- -target(x0) / tempseq[k] + 0.5*sum(p0^2)
            xpprop <- lf(x0, p0, gd=gd.target, jsize=jsizes[k]*runif(1,0.95,1.05), njumps=njumps, temperature=tempseq[k])
            Hprop <- -target(xpprop[1:x.d]) / tempseq[k] + 0.5*sum(xpprop[x.d+1:x.d]^2)
            if (!is.na(-Hprop+H0) && runif(1) < exp(-Hprop+H0)) {
                X[,k] <- xpprop[1:x.d]
            } else {
                X[,k] <- x0
            }
            ## update leapfrog step size
            if (is.na(Hprop-H0)) {
                jsizes[k] <- exp(log(jsizes[k]) - (1/(1+i)^0.6))
            } else {
                jsizes[k] <- exp(log(jsizes[k]) + (1/(1+i)^0.6) * (min(1, exp(-Hprop+H0)) - 0.9))
            }
        }
        ## exchange states
        target_accprob <- 0.234
        acc_prob <- numeric(ntemp-1)
        nonreversible <- TRUE # non-reversible chain? (if TRUE, use the even-odd temperature swap schedule as in Syed et al.)
        if (nonreversible) {
            if (i%%2==0) {
                lowertemp <- seq(2,ntemp-1,by=2)
            } else {
                lowertemp <- seq(1,ntemp-1,by=2)
            }
        } else {
            evenswap <- sample(c(T,F), 1)
            if (evenswap) {
                lowertemp <- seq(2,ntemp-1,by=2)
            } else {
                lowertemp <- seq(1,ntemp-1,by=2)
            }
        }
        for (k in lowertemp) {
            acc_prob[k] <- min(1, exp((target(X[,k+1])-target(X[,k]))*(1/tempseq[k]-1/tempseq[k+1])))
            if (runif(1) < acc_prob[k]) {
                lowerstate <- X[,k+1]
                X[,k+1] <- X[,k]
                X[,k] <- lowerstate
            }
            rho[k] <- rho[k] + (acc_prob[k] - target_accprob)/(i+1)^0.6 # adaptive update of temperature levels
        }
        for (k in 2:ntemp) {
            tempseq[k] <- tempseq[k-1] + exp(rho[k-1])
        }
        ## record the summary statitsic for each chain & cumulative number of leapfrog steps
        sumstat[[i+1]] <- sapply(1:ntemp, function(k) { sumstat_fun(X[,k]) })
        cumnjump[i] <- ifelse(i==1, 0, cumnjump[i-1]) + ntemp*njumps # cumulative number of leapfrog steps
        topchain[,(i-1)%%span+1] <- X[,ntemp] # replace the current value of the top chain in the record book (for checking the search criterion)
        ## Decide whether to add or remove the top chain
        if (!NTEMP_TUNED) {
            if (ntemp_tuning_method=="rectangular") {
                REACHED <- mean(apply(topchain, 1, function(vec) { any(vec < search_center - search_scale) && any(vec > search_center + search_scale) })) >= 0.5 # check if the search criterion is met for at least half of the coordinates at the top level chain
                if (!is.na(REACHED)) { REACHED_TIMES <- REACHED_TIMES + REACHED }
                if (REACHED_TIMES >= span) { NTEMP_TUNED <- TRUE; REACHED <- NA } # no further tuning if REACHED=TRUE for at least `span` times
                if (i%%ceiling(span/2)==0 && !is.na(REACHED) && !REACHED) { # add another chain
                    X <- cbind(X, X[,ntemp]) # initialize the new top chain to be equal to the previous top chain
                    rho[ntemp] <- rho[ntemp-1]
                    tempseq[ntemp + 1] <- tempseq[ntemp] + exp(rho[ntemp])
                    jsizes[ntemp + 1] <- jsizes[ntemp]
                    ntemp <- length(tempseq) # ntemp increases by 1
                }
            }
            if (ntemp_tuning_method=="Uthreshold") {
                REACHED <- any(apply(topchain, 1, target) < -Uthreshold, na.rm=TRUE) # check if the search criterion is met for at least half of the coordinates at the top level chain
                if (!is.na(REACHED)) { REACHED_TIMES <- REACHED_TIMES + REACHED }
                if (REACHED_TIMES >= span) { NTEMP_TUNED <- TRUE; REACHED <- NA } # no further tuning if REACHED=TRUE for at least `span` times
                if (i%%ceiling(span/2)==0 && !is.na(REACHED) && !REACHED) { # add another chain
                    X <- cbind(X, X[,ntemp]) # initialize the new top chain to be equal to the previous top chain
                    rho[ntemp] <- rho[ntemp-1]
                    tempseq[ntemp + 1] <- tempseq[ntemp] + exp(rho[ntemp])
                    jsizes[ntemp + 1] <- jsizes[ntemp]
                    ntemp <- length(tempseq)
                }
            }
        }
        tempseq_history[[i+1]] <- tempseq
        acc_prob_history[[i]] <- acc_prob
        if (verbose) { cat("i: ", i, " k: ", k, " ntemp: ", ntemp, " maxtemp: ", tempseq[k], " NTEMP_TUNED: ", NTEMP_TUNED, " REACHED: ", REACHED, "\n") }
    }
    return(list(jsizes=jsizes, tempseq=tempseq_history, acc_prob=acc_prob_history, sumstat=sumstat, cumnjump=cumnjump))
}

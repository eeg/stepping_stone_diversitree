### Xie et al. 2011. SystBiol. 10.1093/sysbio/syq085
### Estimation of marginal likelihood via the "stepping stone" algorithm.
###
### Dec 2013. Implemented with diversitree. Emma Goldberg.

require(diversitree)

#--------------------------------------------------
# Functions for wrapping the stepping stone algorithm around diversitree
#--------------------------------------------------

# The "power-likelihood".  Returns log(likelihood ^ beta)
loglik.beta <- function(loglik, b)
{
    llb <- function(pars)
        loglik(pars) * b
    diversitree:::set.info(llb, diversitree:::get.info(loglik))
    attributes(llb) <- attributes(loglik)
    return(llb)
}

# Obtain the power-posterior samples for all beta_k.
# (All we really need are the log-likelihoods, but for now it's instructive to
# return the complete samples.)
sample.pp.all <- function(beta.all, mcmc.args)
{
    beta.all <- sort(beta.all, decreasing=T)
    stopifnot(beta.all[1] == 1 & beta.all[length(beta.all)] == 0)

    burnin <- 100
    mcmc.args$nsteps <- mcmc.args$nsteps + burnin

    loglik <- mcmc.args$lik
    an <- argnames(loglik)

    ans <- list()
    
    # Sample the real, not-power posterior.  (Could omit this---SS does not use
    # the samples from k = K, and one would probably have them already anyway.)
    samp <- do.call("mcmc", mcmc.args)
    ans[[1]] <- samp[-seq(burnin), -1]

    # Get also the log-likelihood values.
    # (This is faster than actually re-computing loglik for each sample.)
    ans[[1]]$ll <- ans[[1]]$p - apply(ans[[1]][,an], 1, mcmc.args$prior)

    # Sample the beta_k power-posteriors.
    # Each beta_k helps to burn in the next.  (Could instead use increasing k,
    # starting from the prior.)
    for (i in seq(2, length(beta.all)))
    {
        mcmc.args$lik <- loglik.beta(loglik, beta.all[i])
        mcmc.args$x.init <- as.numeric(ans[[i-1]][nrow(ans[[i-1]]), an])
        samp <- do.call("mcmc", mcmc.args)
        ans[[i]] <- samp[-seq(burnin), -1]
        ans[[i]]$ll <- (ans[[i]]$p - apply(ans[[i]][,an], 1,
                                           mcmc.args$prior)) / beta.all[i]
    }

    ans <- rev(ans)
    beta.all <- rev(beta.all)
    names(ans) <- as.character(beta.all)

    # Need to repair ll for beta_k = 0.
    ans[[1]]$ll <- apply(ans[[1]], 1, function(x) loglik(as.numeric(x[an])))

    vals <- list(ll = sapply(ans, function(x) x$ll), 
                beta = beta.all, pp.samples = ans)
    return(vals)
}

# Estimate the (log) marginal likelihood from the power-posterior samples.
get.logml.ss <- function(beta.all, ll.all)
{
    stopifnot(beta.all[1] == 0 & beta.all[length(beta.all)] == 1)

    # Compute log(r_SS,k) for a particular k.
    #   beta.diff = beta_k - beta_k-1
    get.logrk <- function(beta.diff, ll.k)
    {
        ll.max <- max(ll.k)
        lrk <- beta.diff * ll.max + log(mean(exp(beta.diff * (ll.k - ll.max))))
        return(lrk)
    }

    # Get all the log(r_SS,k).
    lrk <- sapply(seq(length(beta.all) - 1), 
              function(k) get.logrk(beta.all[k+1] - beta.all[k], ll.all[, k]))

    # Add up all k's to get log(r_SS), the estimate of the marginal likelihood.
    return(sum(lrk))
}

#--------------------------------------------------
# Test data and preparations.
#--------------------------------------------------

# A birth-death tree.
pars <- c(1, 0.3)  # birth rate, death rate
set.seed(1)
phy <- tree.bd(pars, max.t=5)

# A likelihood function.
loglik <- make.bd(phy)

# Rates for exponential priors.
prior.rates <- c("b" = 0.1, "d" = 0.2)

# Stepping stone parameters.
K <- 50
b.unif <- seq(0, 1, length.out = K + 1)
b.beta <- qbeta(b.unif, 0.3, 1)

nsteps <- 2000  # number of MCMC iterations, after burnin
burnin <- 200

#--------------------------------------------------
# Run.
#--------------------------------------------------

mcmc.args <- list(lik=loglik, x.init=pars, nsteps=nsteps,
                  prior=make.prior.exponential(prior.rates), 
                  w=1, lower=0, upper=10, print.every=1000)

set.seed(2)
ans <- sample.pp.all(b.beta, mcmc.args)
# From ans$pp.samples, can check for burnin.  And can see that for the
# parameters, pp.samples[["1"]] looks like a real posterior while
# pp.samples[["0"]] looks like the prior.

ans$lml <- get.logml.ss(ans$beta, ans$ll)
# 250.9406

# todo: should adjust w for each beta_k -- will run faster

#--------------------------------------------------
# Compare with inferior methods.
#--------------------------------------------------

library(LaplacesDemon)

aname <- argnames(loglik)

# Use lots of samples from the true posterior.
mcmc.args$nsteps <- 2000 * (K+1)
set.seed(2)
ansL <- do.call("mcmc", mcmc.args)
ansL$ll <- ansL$p - apply(ansL[,aname], 1, mcmc.args$prior)

# Harmonic mean
LML(theta = as.matrix(ansL[, aname]), LL = ansL$ll, method="HME")$LML
# 254.0632

# Arrogance sampling
LML(theta = as.matrix(ansL[, aname]), LL = ansL$ll, method="NSIS")$LML
# 254.4796

# The answers from these "easy" methods are a bit too large.  This is what Xie
# et al. 2010 report, too, at least for the HME.

# I also compared with the steppingStoneSampling() function in the TESS
# library.  After a bugfix, it yields 250.1622 for this test data -- similar
# to my SS function, less like the LaplacesDemon functions.

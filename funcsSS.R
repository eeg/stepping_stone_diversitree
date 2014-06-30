require(diversitree)

# Rather than load the whole TESS library, we need only the
# tess.steppingStoneSampling() function.  This slightly-modified version
# returns its likelihood samples, and it includes a small correction (I think)
# to the marginal likelihood calculation.
source("tessSSSmod.R")

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
        ans[[i]]$ll <- (ans[[i]]$p - apply(ans[[i]][,an], 1, mcmc.args$prior)) / beta.all[i]
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

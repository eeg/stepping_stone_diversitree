### Xie et al. 2011. SystBiol. 10.1093/sysbio/syq085
### Estimation of marginal likelihood via the "stepping stone" algorithm.
###
### Dec 2013.  Emma Goldberg.

###--------------------------------------------------
### Normal distribution
###--------------------------------------------------
### Their test case.  See page 154 and Table 1.

#--------------------------------------------------
# Test data
#--------------------------------------------------

# "A single data set of size n = 100 was simulated from a normal distribution
# having mean mu = 0.0 and standard deviation tau = 1.0."

p <- list(n = 100, mu = 0, tau = 1)     # parameters
set.seed(6)                             # gets the marginal likelihood close to their value
y <- rnorm(p$n, mean=p$mu, sd=p$tau)    # the "data" itself

# Lower-left of page 154.
p$sig2 <- p$tau^2 / p$n
p$s2 <- (sum(y^2) - p$n * mean(y)^2) / p$n

# "The prior for mu was normal with mean mu_0 = 0.0 and standard deviation
# sigma_0 = 1.0."
p$mu.0 <- 0
p$sig2.0 <- 1

# True log-marginal likelihood, Eq (2)
get.logml <- function(dat, p)
{
    ybar <- mean(dat)
    ans <- (-p$n/2) * log(2 * pi * p$n * p$sig2) +
        (-1/2) * log(p$sig2.0 / p$sig2 + 1) +
        (-1/2) * ((p$s2 + ybar^2) / p$sig2 + 
                  p$mu.0^2 / p$sig2.0 - 
                  (ybar/p$sig2 + p$mu.0/p$sig2.0)^2 / (1/p$sig2 + 1/p$sig2.0))
    return(ans)
}

logml.true <- get.logml(y, p)
logml.true  # -147.0712
# Not exactly their number (middle of page 154 column 2), but that's due to
# random sampling.

# Log-likelihood of the test data.  They come from a normal distribution with
# mean to be estimated.
make.loglik <- function(dat, param)
{
    function(param)
        sum(log(dnorm(dat, mean=param, sd=1)))
}
loglik <- make.loglik(y, param)

# The power posterior.  MCMC is not used -- the distribution can be written and
# sampled from explicitly (very bottom-left of p. 154).
# ["b" is \beta; "beta" is an R function]
sample.pp <- function(b, n, dat, p)
{
    ybar <- mean(dat)
    m.pp <- (b * ybar / p$sig2 + p$mu.0/p$sig2.0) / (b/p$sig2 + 1/p$sig2.0)
    s.pp <- 1 / sqrt(b / p$sig2 + 1 / p$sig2.0)
    rnorm(n, mean=m.pp, sd=s.pp)
}
n.samp <- 2000

#--------------------------------------------------
# K and beta
#--------------------------------------------------
# Mainly for stepping-stone, but also use K for harmonic mean.
# "b" is \beta; the "beta" function happens to be used to produce values of b

K <- 50
b.unif <- seq(0, 1, length.out = K + 1)
b.beta <- qbeta(b.unif, 0.3, 1)             # more samples near b = 0

#--------------------------------------------------
# Harmonic mean
#--------------------------------------------------

# "Harmonic mean analyses were based on 2000(K+1) posterior samples"
set.seed(1)
samp.hm <- sample.pp(1, n.samp * (K + 1), y, p)

ll <- sapply(samp.hm, loglik)

logml.hm <- log(1/(mean(1/exp(ll))))
logml.hm # -146.0855
# Too high by about 1, as in their example.

#--------------------------------------------------
# Stepping stone
#--------------------------------------------------

# Estimate the (log) marginal likelihood via the stepping stone algorithm.
# Top right of page 153.
get.logml.ss <- function(betas, samp.par, loglik)
{
    # Compute log(r_SS,k) for a particular k.
    #   beta.diff = beta_k - beta_k-1
    #   samp = parameter samples from the power posterior for this beta_k
    get.logrk <- function(beta.diff, samp)
    {
        # likelihood of the data given the model (i.e., the parameter samples)
        ll <- sapply(samp, loglik)
        
        # log(r_SS,k)
        ll.max <- max(ll)
        lrk <- beta.diff * ll.max + log(mean(exp(beta.diff * (ll - ll.max))))
        return(lrk)
    }

    # Get all the log(r_SS,k).
    lrk <- sapply(seq(length(betas) - 1), 
                  function(k) get.logrk(betas[k+1] - betas[k], samp.par[, k]))
    print(lrk)

    # Add up all k's to get log(r_SS), the estimate of the marginal likelihood.
    return(sum(lrk))
}

# "We sampled 2000 values of mu directly from each power posterior"
set.seed(2)
samp.unif <- sapply(b.unif, sample.pp, n.samp, y, p)
samp.beta <- sapply(b.beta, sample.pp, n.samp, y, p)

logml.ss.unif <- get.logml.ss(b.unif, samp.unif, loglik)
logml.ss.unif   # -147.0821
# Very close to the true value.

logml.ss.beta <- get.logml.ss(b.beta, samp.beta, loglik)
logml.ss.beta   # -147.0628
# Even closer, as in their example.

#--------------------------------------------------
# Summarize the incorrectness.
#--------------------------------------------------

# This looks like a great package.  Include its MCMC-based methods.
library(LaplacesDemon)

LML(theta=samp.hm, LL=ll, method="HME")$LML # -146.0855
# Identical to logml.hm above

logml.nsis <- LML(theta=matrix(samp.hm), LL=ll, method="NSIS")$LML # -146.159
# Equivalent to margLikArrogance package.

round(c(HME = logml.hm, NSIS = logml.nsis, SS.unif = logml.ss.unif, SS.beta = logml.ss.beta) - logml.true, 4)
#    HME    NSIS SS.unif SS.beta 
# 0.9857  0.9122 -0.0109  0.0085 

# Well, that arrogance-sampling method wasn't very good!
# But stepping stone performs very well, as in their example.


###--------------------------------------------------
### Phylogenetic data
###--------------------------------------------------
### An example with diversitree.

library(diversitree)

#--------------------------------------------------
# Functions
#--------------------------------------------------

# Returns log(likelihood ^ beta)
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
    # the samples from k = K, and would probably have them already anyway.)
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
                beta.all = beta.all, pp.samples = ans)
    return(vals)
}

# Estimate the (log) marginal likelihood via the stepping stone algorithm.
# Top right of page 153.
get.logml.ss <- function(beta.all, ll.all)
{
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
# Try it out
#--------------------------------------------------

set.seed(1)
pars <- c(1, 0.3)
phy <- tree.bd(pars, max.t=5)

loglik <- make.bd(phy)
prior <- make.prior.exponential(c(0.1, 0.2))

K <- 50
b.unif <- seq(0, 1, length.out = K + 1)
b.beta <- qbeta(b.unif, 0.3, 1)

mcmc.args <- list(lik=loglik, x.init=pars, nsteps=2000, w=1, prior=prior,
                  lower=0, upper=10, print.every=1000)
samp.all <- sample.pp.all(b.beta, mcmc.args)

# From samp.all$pp.samples, can check for burnin, and can see that
# pp.samples[["1"]] looks like a real posterior while pp.samples[["0"]] looks
# like the prior.

logml.ss <- get.logml.ss(samp.all$beta.all, samp.all$ll)
logml.ss # 250.9399

#--------------------------------------------------
# Compare with the "easy" methods
#--------------------------------------------------

library(LaplacesDemon)
an <- argnames(loglik)

# Straight from the not-power posterior
ans.1 <- samp.all$pp.samples[["1"]]
LML(theta = as.matrix(ans.1[, an]), LL = ans.1$ll, method="HME")$LML  # 254.2774
LML(theta = as.matrix(ans.1[, an]), LL = ans.1$ll, method="NSIS")$LML # 254.4239

# Or better, use more samples
mcmc.args$nsteps <- 2000 * (K+1)
ans <- do.call("mcmc", mcmc.args)
ans$ll <- ans$p - apply(ans[,an], 1, mcmc.args$prior)

LML(theta = as.matrix(ans[, an]), LL = ans$ll, method="HME")$LML  # 253.8072
LML(theta = as.matrix(ans[, an]), LL = ans$ll, method="NSIS")$LML # 254.3445

# As before, the answers from these easy methods are a bit too large.

#--------------------------------------------------
# Compare with TESS
#--------------------------------------------------

# library(TESS)
# Or source the modified function instead.

priors <- c(function(x) dexp(x,rate=0.1,log=TRUE), function(x) dexp(x,rate=0.2,log=TRUE))
ans.tess <- tess.steppingStoneSampling(loglik, priors, pars, c(TRUE,TRUE), 2000, burnin=100, K=50) # 247.8093

# Pretty close.  But who is right?!

# Modify the tess function to return its samples.
# fixInNamespace("tess.steppingStoneSampling", "TESS")  # see other file

ans <- steppingStoneSampling(loglik, priors, pars, c(TRUE,TRUE), 10000, burnin=100, K=50) # 247.964

plot(ans$samples[,1])
# Browsing through these, the tess sampler definitely wanders off into lower-probability regions more than diversitree's slice sampler does.

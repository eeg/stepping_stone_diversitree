### Xie et al. 2011. SystBiol. 10.1093/sysbio/syq085
### Estimation of marginal likelihood via the "stepping stone" algorithm.
###
### Comparing my implementation with diversitree to the TESS implementation.
### Dec 2013.  Emma Goldberg.

source("funcsSS.R")

#--------------------------------------------------
# Test data and common preparations.
#--------------------------------------------------

# A birth-death tree.
pars <- c(1, 0.3)  # birth rate, death rate
set.seed(1)
phy <- tree.bd(pars, max.t=5)

# A likelihood function.  Used with both the diversitree and tess methods.
loglik <- make.bd(phy)

# For both methods, use exponential priors with these rates.
prior.rates <- c("b" = 0.1, "d" = 0.2)

# Stepping stone parameters.  Using the same as inside the tess function.
K <- 50
b.unif <- seq(0, 1, length.out = K + 1)
b.beta <- qbeta(b.unif, 0.3, 1)

nsteps <- 2000  # number of MCMC iterations, after burnin
burnin <- 200

#--------------------------------------------------
# With diversitree
#--------------------------------------------------

mcmc.args <- list(lik=loglik, x.init=pars, nsteps=nsteps,
                  prior=make.prior.exponential(prior.rates), 
                  w=1, lower=0, upper=10, print.every=1000)

set.seed(2)
ans.divt <- sample.pp.all(b.beta, mcmc.args)
# From ans.divt$pp.samples, can check for burnin.  And can see that for the
# parameters, pp.samples[["1"]] looks like a real posterior while
# pp.samples[["0"]] looks like the prior.

ans.divt$lml <- get.logml.ss(ans.divt$beta, ans.divt$ll)
# 250.9406

# todo: should adjust w for each beta_k -- will run faster

#--------------------------------------------------
# With tess
#--------------------------------------------------

priors.tess <- c(function(x) dexp(x, rate=prior.rates["b"], log=TRUE),
                 function(x) dexp(x, rate=prior.rates["d"], log=TRUE))

set.seed(2)
ans.tess <- steppingStoneSampling(loglik, priors.tess, pars,
                                  logTransforms=c(TRUE, TRUE),
                                  iterations=nsteps, burnin=burnin)
ans.tess$lml
# 250.1622 with my edit to the function
ans.tess$lml0
# 247.7725 with the original version

# Note that the above value includes a small bugfix (I think) in the
# computation of "BF" near the end of the tess function.  I noticed a
# discrepancy when comparing answers from running the tess log-likelihood
# samples through my own function for computing the marginal likelihood:
get.logml.ss(rev(ans.tess$beta), ans.tess$ll[, rev(seq(ncol(ans.tess$ll)))])
# 250.1622

# I haven't checked if the (corrected) tess lml values are systematically
# slightly-lower than the diversitree-based ones.  They might be---a brief
# glance at the samples indicates that the tess sampler may be more prone to
# diversions into lower-probability parameter space.  On the other hand, the
# tess method seems to be significantly faster in cranking out samples.

#--------------------------------------------------
# Inferior methods
#--------------------------------------------------

library(LaplacesDemon)

aname <- argnames(loglik)

# Use lots of samples from the true posterior.
mcmc.args$nsteps <- 2000 * (K+1)
set.seed(2)
ans <- do.call("mcmc", mcmc.args)
ans$ll <- ans$p - apply(ans[,aname], 1, mcmc.args$prior)

# Harmonic mean
LML(theta = as.matrix(ans[, aname]), LL = ans$ll, method="HME")$LML
# 254.0632

# Arrogance sampling
LML(theta = as.matrix(ans[, aname]), LL = ans$ll, method="NSIS")$LML
# 254.4796

# The answers from these "easy" methods are a bit too large.  This is what Xie
# et al. 2010 report, too, at least for the HME.

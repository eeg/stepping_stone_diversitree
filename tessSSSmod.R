message("Loading a slightly-modified form of TESS' steppingStoneSampling()")
steppingStoneSampling <- function(likelihoodFunction,priors,parameters,logTransforms,iterations,burnin=round(iterations/3),K=50) {

  x <- K:0 / K
  beta <- qbeta(x,0.3,1)

  # pre-compute current posterior probability
  pp    <- likelihoodFunction(parameters)
  prior <- 0
  for ( j in 1:length(parameters) ) {
    prior <- prior + priors[[j]](parameters[j])
  }

  # the path values
  BF <- 0
  BFeeg <- 0 # NEW
  ll.samples <- matrix(NA, nrow=iterations, ncol=K+1) # NEW

  for (k in 1:length(beta)) {
    b <- beta[k]

    samples <- c()
    for (i in 1:(iterations+burnin)) {
    
      # propose new values
      for ( j in 1:length(parameters) ) {
        new_prior <- prior - priors[[j]](parameters[j])
        if ( logTransforms[j] == TRUE ) {
          if (parameters[j] == 0) {
            stop("Cannot propose new value for a parameter with value 0.0.")
          }
          eta           <- log(parameters[j]) ### propose a new value for parameter[j]
          new_eta       <- eta + rnorm(1,0,1)
          new_val       <- exp(new_eta)
          hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
          parameters[j] <- new_val
          new_pp        <- likelihoodFunction(parameters)
          new_prior     <- new_prior + priors[[j]](parameters[j])
        
          # compute acceptance ratio for power posterior
          if ( b == 0.0 ) {
            acceptance_ratio <- new_prior-prior+hr
          } else { 
            acceptance_ratio <- new_prior-prior+b*new_pp-b*pp+hr
          }
          # accept / reject
          if (is.finite(new_pp) && is.finite(new_prior) && acceptance_ratio > log(runif(1,0,1)) ) {
            pp <- new_pp
            prior <- new_prior
          } else {
            parameters[j] <- exp(eta)
          }
        } else {
          eta           <- parameters[j] ### propose a new value for parameter[j]
          new_val       <- eta + rnorm(1,0,1)
          hr            <- 0.0 # calculate the Hastings ratio
          parameters[j] <- new_val
          new_pp        <- likelihoodFunction(parameters)
          new_prior     <- new_prior + priors[[j]](parameters[j])
          # compute acceptance ratio for power posterior
          if ( b == 0.0 ) {
            acceptance_ratio <- new_prior-prior+hr
          } else { 
            acceptance_ratio <- new_prior-prior+b*new_pp-b*pp+hr
          }
          # accept / reject
          if (is.finite(new_pp) && is.finite(new_prior) && acceptance_ratio > log(runif(1,0,1)) ) {
            pp <- new_pp
            prior <- new_prior
          } else {
            parameters[j] <- eta
          }
        }

      }

      # sample the likelihood
      if (i > burnin) {
        samples[i-burnin] <- pp
      }
    }
    ll.samples[, k] <- samples

    max <- samples[which.max(samples)]
    if ( k > 1 ) {
       BF <- BF + log(mean( exp(samples-max)^(beta[k-1]-beta[k]) ))+(beta[k-1]-beta[k])*max
       BFeeg <- BFeeg + log(mean( exp((samples-max)*(beta[k-1]-beta[k])) ))+(beta[k-1]-beta[k])*max
       # Referring to the equation for $\log \hat{r}_{SS}$ in the upper-right
       # of page 153, I have made a modification to the grouping of the 
       # "exp{ ... }" part near the end.  Does this look right to you?
       # Actually, they're mathematically the same but BFeeg is numerically more stable.
    }
                              
  }
  return(list(lml = BFeeg, lml0 = BF, ll = ll.samples, beta = beta))
}

# From: Sebastian Höhna <hoehna@math.su.se>
# It makes a lot of sense to store the samples from the power posterior run
# and then to compute both the path-sampling as well as
# stepping-stone-sampling approximations of the marginal likelihood. You
# basically just need to add the line
# BFPS <- BFPS + ( mean(ll.samples[, k-1]) + mean(ll.samples[, k])
# )*(beta[k-1]-beta[k]) / 2.0
# This line should be added just below "BF <- BF + ...". Using the
# path-sampling should also give you another benchmark for the marginal
# likelihood estimator.
# Let me know if you need any help with that but it looks that you already
# adapted my code for this.

#' Estimating a dynamic mixture via the noisy Cross-Entropy method
#'
#' This function estimates a dynamic mixture by means of the noisy Cross-Entropy method.
#' Currently only implemented for the lognormal - generalized Pareto case,
#' with Cauchy or exponential weight. This function
#' is mostly an auxiliary function, not suitable for the final user. Use CeNoisyFitBoot.R instead.
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset. 
#' @param rawdata either a list of vectors or a vector: in the former case, each
#' vector contains a dataset to be used for estimation.
#' @param rho real in (0,1): parameter determining the quantile of the log-likelihood values to be used at each iteration.
#' @param maxiter non-negative integer: maximum number of iterations.
#' @param alpha real in (0,1): smoothing parameter.
#' @param nsim non-negative integer: number of replications used in the normal and lognormal updating.
#' @param nrepsInt non-negative integer: number of replications used in the Monte Carlo estimate of the normalizing constant.
#' @param eps non-negative real: tolerance for the stopping criterion of the noisy Cross-Entropy method.
#' @param r positive integer: length of window to be used in the stopping criterion.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return For each dataset, a list with the following elements is returned:
#'
#' V (nreps x 12) matrix: updated mean and variance of the distributions used in the stochastic program.
#' nit (positive integer): number of iterations needed for convergence.
#' loglik (scalar): maximized log-likelihood.

#' @details See Rubinstein and Kroese (2004, chap. 6).
#' @keywords dynamic mixture; cross-entropy.
#' @export
#' @examples
#' maxiter = 10
#' alpha = .5
#' rho = .05
#' eps = 1e-2
#' nsim = 1000
#' nrepsInt = 1000
#' res <- CENoisyFit(1,Metro2019,rho,maxiter,alpha,nsim,nrepsInt,eps)
#'
#'@seealso CENoisyFitBoot
#' @references{
#'   \insertRef{rubio4}{FitDynMix}
#' }
#'
#'
#' @importFrom Rdpack reprompt

CENoisyFit <- function(x,rawdata,rho,maxiter,alpha,nsim,nrepsInt,eps,r=5,weight)
{
  if (is.list(rawdata) == TRUE)
    yObs = rawdata[[x]]
  else
    yObs = rawdata
  n = length(yObs)
  mediana = median(yObs)
  y1 = yObs[yObs < mediana]
  y2 = yObs[yObs >= mediana]
  mu0 = MASS::fitdistr(y1, "lognormal")$estimate[1]
  sigma0 = MASS::fitdistr(y1, "lognormal")$estimate[2]
  xi0Est = evir::gpd(y2, mediana)$par.ests["xi"]
  xi0 = pmax(xi0Est,.01)
  beta0 = evir::gpd(y2, mediana)$par.ests["beta"]
  
  if (weight == 'cau')
  {
    muc0 = quantile(yObs, 0.5)
    tau0 = abs(log(sd(yObs)/2))
    v0 = c(as.double(muc0), 1, tau0, 1.5, as.double(mu0), 1, as.double(sigma0), 1, as.double(log(xi0)), 2, as.double(log(beta0)), 2)
    
    X = matrix(0,nsim,6) # columns = number of parameters to be estimated
    v = matrix(0,nsim,12) # columns = number of parameters of the instrumental distributions
    v[1,] = v0
    loglik = matrix(0,nsim,1)
    gammma = matrix(0,nsim,1)
    Bt = matrix(0,nsim,1)
    nit = 2
    change = 100
    while (change > eps & nit <= maxiter)             # start iterations
    {
      X[,1] = rnorm(nsim,v[nit-1,1],v[nit-1,2]) # CA1
      X[,2] = rlnorm(nsim,v[nit-1,3],v[nit-1,4])# CA2
      X[,3] = rnorm(nsim,v[nit-1,5],v[nit-1,6]) # mu
      X[,4] = rlnorm(nsim,v[nit-1,7],v[nit-1,8]) # sigma
      X[,5] = rlnorm(nsim,v[nit-1,9],v[nit-1,10]) # xi
      X[,6] = rlnorm(nsim,v[nit-1,11],v[nit-1,12]) # beta
      for (i in 1:nsim)
      {
        loglik[i] = dynloglikMC(X[i,],yObs,nrepsInt,'cau')
      }
      gammma[nit] = quantile(loglik,1-rho,na.rm=TRUE)
      indici = which(loglik>gammma[nit])
      v[nit,1] = mean(X[indici,1])
      v[nit,2] = sd(X[indici,1])
      v[nit,3] = mean(log(X[indici,2]))
      v[nit,4] = sd(log(X[indici,2]))
      v[nit,5] = mean(X[indici,3])
      v[nit,6] = sd(X[indici,3])
      v[nit,7] = mean(log(X[indici,4]))
      v[nit,8] = sd(log(X[indici,4]))
      v[nit,9] = mean(log(X[indici,5]))
      v[nit,10] = sd(log(X[indici,5]))
      v[nit,11] = mean(log(X[indici,6]))
      v[nit,12] = sd(log(X[indici,6]))
      v[nit,1] = alpha * v[nit,1] + (1-alpha) * v[nit-1,1]
      v[nit,2] = alpha * v[nit,2] + (1-alpha) * v[nit-1,2]
      v[nit,3] = alpha * v[nit,3] + (1-alpha) * v[nit-1,3]
      v[nit,4] = alpha * v[nit,4] + (1-alpha) * v[nit-1,4]
      v[nit,5] = alpha * v[nit,5] + (1-alpha) * v[nit-1,5]
      v[nit,6] = alpha * v[nit,6] + (1-alpha) * v[nit-1,6]
      v[nit,7] = alpha * v[nit,7] + (1-alpha) * v[nit-1,7]
      v[nit,8] = alpha * v[nit,8] + (1-alpha) * v[nit-1,8]
      v[nit,9] = alpha * v[nit,9] + (1-alpha) * v[nit-1,9]
      v[nit,10] = alpha * v[nit,10] + (1-alpha) * v[nit-1,10]
      v[nit,11] = alpha * v[nit,11] + (1-alpha) * v[nit-1,11]
      v[nit,12] = alpha * v[nit,12] + (1-alpha) * v[nit-1,12]
      change = max(c(v[nit,2],v[nit,4],v[nit,6],v[nit,8],v[nit,10],v[nit,12]))
      if (nit>10)
      {
        Bt[nit] = mean(abs(gammma[(nit-10):nit]))
      }
      if (nit > 10 + r)
      {
        Bmin = min(Bt[(nit-r):nit])
        Bmax = max(Bt[(nit-r):nit])
        change = (Bmax-Bmin)/Bmin
      }
      nit = nit + 1
      if (nit < maxiter && change < eps)
      {
        results = list(v[1:(nit-1),],(nit-1),sum(loglik))
        break
      }
      if (nit >= maxiter)
      {
        results = list(V=v[1:(maxiter-1),],nit=(nit-1),loglik=sum(loglik))
        break
      }
      
      if (nit>10)
      {
        Bt[nit] = mean(abs(gammma[(nit-10):nit]))
      }
      if (nit > 10 + r)
      {
        Bmin = min(Bt[(nit-r):nit])
        Bmax = max(Bt[(nit-r):nit])
        change = (Bmax-Bmin)/Bmin
      }
      if (nit < maxiter && change < eps)
      {
        results = list(v[1:(nit-1),],(nit-1),sum(loglik))
        break
      }
      if (nit >= maxiter)
      {
        results = list(V=v[1:(maxiter-1),],nit=(nit-1),loglik=sum(loglik))
        break
      }
    }
  }

  if (weight == 'exp')
  {
    lambda0 = abs(log(sd(yObs)/2))
    v0 = c(as.double(lambda0), 1, as.double(mu0), 1, as.double(sigma0), 1, as.double(log(xi0)), 2, as.double(log(beta0)), 2)
    
    X = matrix(0,nsim,5) # columns = number of parameters to be estimated
    v = matrix(0,nsim,10) # columns = number of parameters of the instrumental distributions
    v[1,] = v0
    loglik = matrix(0,nsim,1)
    gammma = matrix(0,nsim,1)
    Bt = matrix(0,nsim,1)
    nit = 2
    change = 100
    while (change > eps & nit <= maxiter)             # start iterations
    {
      X[,1] = rlnorm(nsim,v[nit-1,1],v[nit-1,2]) # lambda
      X[,2] = rnorm(nsim,v[nit-1,3],v[nit-1,4]) # mu
      X[,3] = rlnorm(nsim,v[nit-1,5],v[nit-1,6]) # sigma
      X[,4] = rlnorm(nsim,v[nit-1,7],v[nit-1,8]) # xi
      X[,5] = rlnorm(nsim,v[nit-1,9],v[nit-1,10]) # beta
      for (i in 1:nsim)
      {
        loglik[i] = dynloglikMC(X[i,],yObs,nrepsInt,'exp')
      }
      gammma[nit] = quantile(loglik,1-rho,na.rm=TRUE)
      indici = which(loglik>gammma[nit])
      v[nit,1] = mean(log(X[indici,1]))
      v[nit,2] = sd(log(X[indici,1]))
      v[nit,3] = mean(X[indici,2])
      v[nit,4] = sd(X[indici,2])
      v[nit,5] = mean(log(X[indici,3]))
      v[nit,6] = sd(log(X[indici,3]))
      v[nit,7] = mean(log(X[indici,4]))
      v[nit,8] = sd(log(X[indici,4]))
      v[nit,9] = mean(log(X[indici,5]))
      v[nit,10] = sd(log(X[indici,5]))
      v[nit,1] = alpha * v[nit,1] + (1-alpha) * v[nit-1,1]
      v[nit,2] = alpha * v[nit,2] + (1-alpha) * v[nit-1,2]
      v[nit,3] = alpha * v[nit,3] + (1-alpha) * v[nit-1,3]
      v[nit,4] = alpha * v[nit,4] + (1-alpha) * v[nit-1,4]
      v[nit,5] = alpha * v[nit,5] + (1-alpha) * v[nit-1,5]
      v[nit,6] = alpha * v[nit,6] + (1-alpha) * v[nit-1,6]
      v[nit,7] = alpha * v[nit,7] + (1-alpha) * v[nit-1,7]
      v[nit,8] = alpha * v[nit,8] + (1-alpha) * v[nit-1,8]
      v[nit,9] = alpha * v[nit,9] + (1-alpha) * v[nit-1,9]
      v[nit,10] = alpha * v[nit,10] + (1-alpha) * v[nit-1,10]
      change = max(c(v[nit,2],v[nit,4],v[nit,6],v[nit,8],v[nit,10]))
    if (nit>10)
    {
      Bt[nit] = mean(abs(gammma[(nit-10):nit]))
    }
    if (nit > 10 + r)
    {
      Bmin = min(Bt[(nit-r):nit])
      Bmax = max(Bt[(nit-r):nit])
      change = (Bmax-Bmin)/Bmin
    }
    nit = nit + 1
    if (nit < maxiter && change < eps)
    {
      results = list(v[1:(nit-1),],(nit-1),sum(loglik))
      break
    }
    if (nit >= maxiter)
    {
      results = list(V=v[1:(maxiter-1),],nit=(nit-1),loglik=sum(loglik))
      break
    }

  if (nit>10)
  {
    Bt[nit] = mean(abs(gammma[(nit-10):nit]))
  }
  if (nit > 10 + r)
  {
    Bmin = min(Bt[(nit-r):nit])
    Bmax = max(Bt[(nit-r):nit])
    change = (Bmax-Bmin)/Bmin
  }
  if (nit < maxiter && change < eps)
  {
    results = list(v[1:(nit-1),],(nit-1),sum(loglik))
    break
  }
  if (nit >= maxiter)
  {
    results = list(V=v[1:(maxiter-1),],nit=(nit-1),loglik=sum(loglik))
    break
  }
  }
  }
  return(results)
}


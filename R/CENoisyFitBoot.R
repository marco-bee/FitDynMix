#' Cross-Entropy estimation and bootstrap standard errors
#'
#' This function estimates a dynamic mixture by means of the noisy Cross-Entropy
#' method and computes bootstrap standard errors.
#' Currently only implemented for the lognormal - generalized Pareto,
#' with Cauchy or exponential weight. Bootstrap standard errors are computed in parallel.
#' @param yObs numerical vector: observed random sample from the mixture.
#' @param nboot integer: number of bootstrap replications for computing the standard errors.
#' If nboot = 0, no standard errors are computed.
#' @param rho real in (0,1): parameter determining the quantile of the log-likelihood values to be used at each iteration.
#' @param maxiter non-negative integer: maximum number of iterations.
#' @param alpha real in (0,1): smoothing parameter.
#' @param nsim non-negative integer: number of replications used in the normal and lognormal updating.
#' @param nrepsInt non-negative integer: number of replications used in the Monte Carlo estimate of the normalizing constant.
#' @param xiInst non-negative real: shape parameter of the instrumental GPD.
#' @param betaInst non-negative real: scale parameter of the instrumental GPD.
#' @param eps non-negative real: tolerance for the stopping criterion of the noisy Cross-Entropy method.
#' @param r positive integer: length of window to be used in the stopping criterion.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return If nboot > 0, a list with the following elements:
#'
#' estPars: Cross-Entropy estimates.
#'
#' nit: number of iterations needed for convergence.
#' 
#' loglik: maximized log-likelihood.
#' 
#' bootPars: parameter estimates obtained for each bootstrap sample.
#' 
#' stddev: bootstrap standard errors.
#' 
#' If nboot = 0, only estPars, nit and loglik are returned.
#' @keywords dynamic mixture; Cross-Entropy; non-parametric bootstrap.
#' @export
#' @examples
#' res = CENoisyFitBoot(Metro2019,0,.05,20,.5,500,500,3,3,.01,5,'exp')

CENoisyFitBoot = function(yObs,nboot,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r=5,weight)
{
  if (weight == 'cau')
  {
  if (nboot >0)
  {
    temp = CENoisyFit(1,yObs,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
    res = temp
    v = res[[1]]
    nitObs = res[[2]]
    loglikObs = res[[3]]
    n_righe = dim(v)[1]
    results = v[n_righe,]
    estPars = c(results[1],exp(results[3]+results[4]^2/2),results[5],exp(results[7]+results[8]^2/2),exp(results[9]+results[10]^2/2),exp(results[11]+results[12]^2/2))
    n = length(yObs)
    Yboot = list()
    for (i in 1:nboot)
    {
      Yboot[[i]] = sample(yObs,n,replace=TRUE)
    }
    nboot.list <- sapply(1:nboot, list)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      n.cores <- 2L
    } else {
      n.cores <- parallel::detectCores()
    }
    clust <- parallel::makeCluster(n.cores)
    temp <- parallel::parLapply(clust,nboot.list,CENoisyFit,Yboot,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
    parallel::stopCluster(cl=clust)
    results = matrix(0,nboot,12)
    bootPars = matrix(0,nboot,6)
    nit = rep(0,nboot)
    loglik = rep(0,nboot)
    for (i in 1:nboot)
    {
      res = temp[[i]]
      v = res[[1]]
      nit[i] = res[[2]]
      loglik[i] = res[[3]]
      n_righe = dim(v)[1]
      results[i,] = v[n_righe,]
      bootPars[i,] = c(results[i,1],exp(results[i,3]+results[i,4]^2/2),results[i,5],exp(results[i,7]+results[i,8]^2/2),exp(results[i,9]+results[i,10]^2/2),exp(results[i,11]+results[i,12]^2/2))
    }
    varcov = cov(bootPars)
    stddev = sqrt(diag(varcov))
    res = list(estPars=estPars,nit=nitObs,loglik=loglikObs,bootPars=bootPars,stddev=stddev)
    return(res)
  }
  if (nboot == 0)
  {
    temp = CENoisyFit(1,yObs,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
    res = temp
    v = res[[1]]
    nit = res[[2]]
    loglik = res[[3]]
    n_righe = dim(v)[1]
    results = v[n_righe,]
    estPars = c(results[1],exp(results[3]+results[4]^2/2),results[5],
    exp(results[7]+results[8]^2/2),exp(results[9]+results[10]^2/2),
    exp(results[11]+results[12]^2/2))
    res = list(estPars=estPars,nit=nit,loglik=loglik)
    return(res)
  }
}

if (weight == 'exp')
  {
  if (nboot >0)
  {
  temp = CENoisyFit(1,yObs,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
  res = temp
  v = res[[1]]
  nitObs = res[[2]]
  loglikObs = res[[3]]
  n_righe = dim(v)[1]
  results = v[n_righe,]
  estPars = c(exp(results[1]+results[2]^2/2),results[3],exp(results[5]+results[6]^2/2),exp(results[7]+results[8]^2/2),exp(results[9]+results[10]^2/2))
  n = length(yObs)
  Yboot = list()
  for (i in 1:nboot)
  {
    Yboot[[i]] = sample(yObs,n,replace=TRUE)
  }
  nboot.list <- sapply(1:nboot, list)
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    n.cores <- 2L
  } else {
    n.cores <- parallel::detectCores()
  }
  clust <- parallel::makeCluster(n.cores)
  temp <- parallel::parLapply(clust,nboot.list,CENoisyFit,Yboot,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
  parallel::stopCluster(cl=clust)
  results = matrix(0,nboot,10)
  bootPars = matrix(0,nboot,5)
  nit = rep(0,nboot)
  loglik = rep(0,nboot)
  for (i in 1:nboot)
  {
    res = temp[[i]]
    v = res[[1]]
    nit[i] = res[[2]]
    loglik[i] = res[[3]]
    n_righe = dim(v)[1]
    results[i,] = v[n_righe,]
    bootPars[i,] = c(exp(results[i,1]+results[i,2]^2/2),results[i,3],exp(results[i,5]+results[i,6]^2/2),exp(results[i,7]+results[i,8]^2/2),exp(results[i,9]+results[i,10]^2/2))
  }
  varcov = cov(bootPars)
  stddev = sqrt(diag(varcov))
  res = list(estPars=estPars,nit=nitObs,loglik=loglikObs,bootPars=bootPars,stddev=stddev)
  return(res)
  }

if (nboot == 0)
{
  temp = CENoisyFit(1,yObs,rho,maxiter,alpha,nsim,nrepsInt,xiInst,betaInst,eps,r,weight)
  res = temp
  v = res[[1]]
  nit = res[[2]]
  loglik = res[[3]]
  n_righe = dim(v)[1]
  results = v[n_righe,]
  estPars = c(exp(results[1]+results[2]^2/2),results[3],exp(results[5]+results[6]^2/2),exp(results[7]+results[8]^2/2),exp(results[9]+results[10]^2/2))
  res = list(estPars=estPars,nit=nit,loglik=loglik)
  return(res)
  }
}
}

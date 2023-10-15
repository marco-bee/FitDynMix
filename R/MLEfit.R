#' Estimating a dynamic mixture via MLE
#'
#' This function fits a dynamic mixture via standard maximum likelihood.
#' Currently only implemented for the lognormal - generalized Pareto case,
#' with Cauchy or exponential weight.
#' @param yObs numerical vector: observed sample.
#' @param bootreps non-negative integer: number of bootstrap replications. If equal to 0, no standard errors are computed.
#' @param intTol non-negative scalar: threshold for stopping the computation of the integral in the normalization
#' constant: if the integral on the interval from n-1 to n is smaller than intTol, the approximation procedure stops.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return 
#'
#' MLEpars vector: maximum likelihood estimates and maximized
#' log-likelihood.
#'
#' MLEboot matrix: maximum likelihood estimates obtained in
#' each bootstrap replication.
#'
#' sdMLE vector: bootstrap standard deviation of the MLEs.
#' 
#' @details Starting values for mu and sigma are the lognormal MLEs computed
#' with the observations below the median. Initial values for xi and
#' tau are the GPD MLEs obtained with the observations above the median.
#' For the location and scale parameter of the Cauchy, we respectively use
#' the first quartile and abs(log(sd(x)/2)). For the parameter of the exponential, we use
#' abs(log(sd(x)/2)).
#' @keywords dynamic mixture; maximum likelihood.
#' @seealso [AMLEfit]
#' @export
#' @examples
#' \donttest{
#' mixFit <- MLEfit(Metro2019,0,,'cau')}
#' @references{
#'   \insertRef{bee22b}{FitDynMix}
#' }
#'
#'
#' @importFrom Rdpack reprompt

MLEfit <- function(yObs,bootreps,intTol=1e-4,weight)
{
  n = length(yObs)
  mediana = median(yObs)
  y1 = yObs[yObs<mediana]
  y2 = yObs[yObs>=mediana]
  mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
  sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
  xi0 = evir::gpd(y2,mediana)$par.ests['xi']
  beta0 = evir::gpd(y2,mediana)$par.ests['beta']
  if (bootreps > 0)
  {
    if (weight == 'cau')
    {	
      muc0 = quantile(yObs,.5)
      tau0 = abs(log(sd(yObs)/2))
      x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
      res <- optim(x0Lik,dynloglik, gr=NULL,yObs,intTol,'cau',method='L-BFGS-B',
                   lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,50),control=list(fnscale=-1))
      estMLE <- c(res$par,res$value) # muc, tau, mu, sigma, xi, beta
      nreps.list <- sapply(1:bootreps, list)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        n.cores <- 2L
      } else {
        n.cores <- parallel::detectCores()
      }
      clust <- parallel::makeCluster(n.cores)
      MLEboot = matrix(0,bootreps,6)
      temp <- parallel::parLapply(clust,nreps.list, MLEBoot,yObs,intTol,'cau')
      parallel::stopCluster(cl=clust)
      for (i in 1:bootreps)
      {
        MLEboot[i,] = as.vector(unlist(temp[[i]]))
      }
      varcov = cov(MLEboot)
      stddev = sqrt(diag(varcov))
      out <- list(pars = estMLE, MLEboot = MLEboot, stddev = stddev)
    }
    if (weight == 'exp')
    {
      lambda0 = abs(log(sd(yObs)/2))
      x0Lik = as.numeric(c(lambda0,mu0,sigma0,xi0,beta0))
      res <- optim(x0Lik,dynloglik, gr=NULL,yObs,intTol,'exp',method='L-BFGS-B',
                   lower=c(.01,-Inf,.01,.01,.01),upper=c(Inf,Inf,Inf,Inf,Inf,Inf),control=list(fnscale=-1))
      estMLE <- c(res$par,res$value) # lambda, mu, sigma, xi, beta
      nreps.list <- sapply(1:bootreps, list)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE") {
        n.cores <- 2L
      } else {
        n.cores <- parallel::detectCores()
      }
      clust <- parallel::makeCluster(n.cores)
      MLEboot = matrix(0,bootreps,5)
      #    temp <- lapply(nreps.list, MLEBoot,yObs,intTol)
      temp <- parallel::parLapply(clust,nreps.list, MLEBoot,yObs,intTol,'exp')
      parallel::stopCluster(cl=clust)
      for (i in 1:bootreps)
      {
        MLEboot[i,] = as.vector(unlist(temp[[i]]))
      }
      varcov = cov(MLEboot)
      stddev = sqrt(diag(varcov))
    }
    out <- list(pars = estMLE, MLEboot = MLEboot, stddev = stddev)
  }  
  if (bootreps == 0)
  {
    if (weight == 'cau')
    {	
      muc0 = quantile(yObs,.5)
      tau0 = abs(log(sd(yObs)/2))
      x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
      res <- optim(x0Lik,dynloglik, gr=NULL,yObs,intTol,'cau',method='L-BFGS-B',
                   lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,50),control=list(fnscale=-1))
      estMLE <- c(res$par,res$value) # muc, tau, mu, sigma, xi, beta
      out <- list(pars = estMLE)
    }
    if (weight == 'exp')
    {	
      lambda0 = abs(log(sd(yObs)/2))
    x0Lik = as.numeric(c(lambda0,mu0,sigma0,xi0,beta0))
    res <- optim(x0Lik,dynloglik, gr=NULL,yObs,intTol,'exp',method='L-BFGS-B',
                 lower=c(.01,-Inf,.01,.01,.01),upper=c(Inf,Inf,Inf,Inf,Inf,Inf),control=list(fnscale=-1))
    estMLE <- c(res$par,res$value) # lambda, mu, sigma, xi, beta
    out <- list(pars = estMLE)
    }
  }
  return(out)
}

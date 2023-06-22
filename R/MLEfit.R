#' Estimating a dynamic mixture via MLE
#'
#' This function fits a dynamic mixture via standard maximum likelihood.
#' Currently only implemented for the lognormal - generalized Pareto case.
#' Bootstrap standard errors are computed via parallel computing.
#' @param yObs numerical vector: observed sample.
#' @param bootreps non-negative integer: number of bootstrap replications. If equal to 0, no standard errors are computed.
#' @param intTol non-negative scalar: threshold for stopping the computation of the integral in the normalization
#' constant: if the integral on the interval from n-1 to n is smaller than intTol, the approximation procedure stops.
#' @return 
#'
#' MLEpars (7 x 1) vector: maximum likelihood estimates and
#' maximized log-likelihood.
#'
#' MLEboot (bootreps x 6) matrix: maximum likelihood estimates obtained in
#' each bootstrap replication.
#'
#' sdMLE (6 x 1) vector: bootstrap standard deviation of the MLEs.
#' 
#' @details Starting values for mu and sigma are the lognormal MLEs computed
#' with the observations below the median. Initial values for xi and
#' tau are the GPD MLEs obtained with the observations above the median.
#' For the location and scale parameter of the Cauchy, we respectively use the first quartile and log(sd(x)/2).
#' @keywords dynamic mixture; maximum likelihood.
#' @seealso [AMLEfit]
#' @export
#' @examples
#' \donttest{
#' mixFit <- MLEfit(Metro2019,0)}
#' @references{
#'   \insertRef{bee22b}{FitDynMix}
#' }
#'
#'
#' @importFrom Rdpack reprompt

MLEfit <- function(yObs,bootreps,intTol=1e-4)
{
  n = length(yObs)
  mediana = median(yObs)
  y1 = yObs[yObs<mediana]
  y2 = yObs[yObs>=mediana]
  mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
  sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
  xi0 = evir::gpd(y2,mediana)$par.ests['xi']
  beta0 = evir::gpd(y2,mediana)$par.ests['beta']
  muc0 = quantile(yObs,.5)
  tau0 = log(sd(yObs)/2)
  x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
  res <- optim(x0Lik,dynloglik, gr=NULL,yObs,intTol,method='L-BFGS-B',lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,50),control=list(fnscale=-1))
  estMLE <- c(res$par,res$value) # muc, tau, mu, sigma, xi, beta
  if (bootreps > 0)
  {
    nreps.list <- sapply(1:bootreps, list)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      n.cores <- 2L
    } else {
      n.cores <- parallel::detectCores()
    }
    clust <- parallel::makeCluster(n.cores)
    MLEboot = matrix(0,bootreps,6)
#    temp <- lapply(nreps.list, MLEBoot,yObs,intTol)
    temp <- parallel::parLapply(clust,nreps.list, MLEBoot,yObs,intTol)
    parallel::stopCluster(cl=clust)
    for (i in 1:bootreps)
    {
      MLEboot[i,] = as.vector(unlist(temp[[i]]))
    }
    varcov = cov(MLEboot)
    stddev = sqrt(diag(varcov))
    out <- list(pars = estMLE, MLEboot = MLEboot, stddev = stddev)
    return(out)
  }
  if (bootreps == 0)
  {
    out <- list(pars = estMLE)
  }
}

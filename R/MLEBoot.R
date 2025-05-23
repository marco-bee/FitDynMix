#' Bootstrap standard errors of MLEs
#'
#' This function creates bootstrap samples of input data and fits a dynamic mixture via
#' maximum likelihood.
#' @param x list of integers: indices of replications.
#' @param y numerical vector: observed data.
#' @param intTol threshold for stopping the computation of the integral in the normalization
#' constant: if the integral on the interval from n-1 to n is smaller than intTol, the approximation procedure stops.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return A list with the following elements:
#'
#' MLE: maximum likelihood estimates obtained from each bootstrap sample.
#'
#' errors: number of times the MLE algorithm breaks down.
#' @details MLEs are computed by means of the optim function. When it breaks
#' down, the sample is discarded and a new one is generated. The function keeps
#' track of the number of times this happens.
#' @export
#' @examples
#' \donttest{bootMLEs <- MLEBoot(1,Metro2019,1e-04,'exp')}

MLEBoot = function(x,y,intTol,weight)
{
  n = length(y)
  if (weight == 'cau')
  {
    i = 1
    j = 1
    while (i <= 1)
    {
      yboot = sample(y,n,replace=TRUE)
      mediana = median(yboot)
      y1 = yboot[yboot<mediana]
      y2 = yboot[yboot>=mediana]
      mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
      sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
      xi0 = evir::gpd(y2,mediana)$par.ests['xi']
      beta0 = evir::gpd(y2,mediana)$par.ests['beta']
      muc0 = quantile(yboot,.25)
      tau0 = abs(log(sd(yboot)/2))
      #      print(intTol)
      x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
      MLE <- tryCatch({
      res <- optim(x0Lik,dynloglik, gr=NULL,yboot,intTol,'cau',method='L-BFGS-B',lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,Inf),control=list(fnscale=-1))
      MLE <- res$par # muc, tau, mu, sigma, xi, beta
      },
      error = function(e) {'errore'}) -> condizione
      if(inherits(condizione,"character",which=TRUE)==0) i = i + 1
      if(inherits(condizione,"character",which=TRUE)==1) j = j + 1
    }
    return(list(MLE=MLE,errori=j))
  }
if (weight == 'exp')
  {
  i = 1
  j = 1
  while (i <= 1)
  {
    yboot = sample(y,n,replace=TRUE)
    mediana = median(yboot)
    y1 = yboot[yboot<mediana]
    y2 = yboot[yboot>=mediana]
    mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
    sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
    xi0 = evir::gpd(y2,mediana)$par.ests['xi']
    beta0 = evir::gpd(y2,mediana)$par.ests['beta']
    lambda0 = log(sd(yboot)/2)
    x0Lik = as.numeric(c(lambda0,mu0,sigma0,xi0,beta0))
    MLE <- tryCatch({
    res <- optim(x0Lik,dynloglik, gr=NULL,yboot,intTol,'exp',method='L-BFGS-B',
                 lower=c(.01,-Inf,.01,.01,.01),upper=c(Inf,Inf,Inf,Inf,Inf,Inf),control=list(fnscale=-1))
    MLE <- res$par # lambda, mu, sigma, xi, beta
  },
  error = function(e) {'errore'}) -> condizione
  if(inherits(condizione,"character",which=TRUE)==0) i = i + 1
  if(inherits(condizione,"character",which=TRUE)==1) j = j + 1
  }
  return(list(MLE=MLE,errori=j))
  }
}

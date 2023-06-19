#' Creating bootstrap samples from observed data and computing MLEs
#'
#' This function creates bootstrap samples of input data and fits a dynamic mixture via
#' standard maximum likelihood.
#' @param x list of integers: indices of replications.
#' @param y numerical vector: observed data.
#' @param intTol threshold for stopping the computation of the integral in the normalization
#' constant: if the integral on the interval from n-1 to n is smaller than intTol, the approximation procedure stops.
#' @return A list with the following elements:
#'
#' MLE: maximum likelihood estimates obtained from each bootstrap sample.
#'
#' errors: number of times the MLE algorithm breaks down.
#' @details MLEs are computed by means of the optim function. When it breaks
#' down, the sample is discarded and a new one is generated. The function keeps
#' track of the number of times this happens.
#' @keywords dynamic mixture; MLE; non-parametric bootstrap.
#' @export
#' @examples
#' bootMLEs <- MLEBoot(1,Metro2019,1e-02)

MLEBoot = function(x,y,intTol)
{
  n = length(y)
  # j = 0
    MLE <- tryCatch({
      yboot = sample(y,n,replace=TRUE)
      mediana = median(yboot)
      y1 = yboot[yboot<mediana]
      y2 = yboot[yboot>=mediana]
      mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
      sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
      xi0 = evir::gpd(y2,mediana)$par.ests['xi']
      beta0 = evir::gpd(y2,mediana)$par.ests['beta']
      muc0 = quantile(yboot,.25)
      tau0 = log(sd(yboot)/2)
#      print(intTol)
      x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
      res <- optim(x0Lik,dynloglik, gr=NULL,yboot,intTol,method='L-BFGS-B',lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,100),control=list(fnscale=-1))
      print(res)
      MLE <- res$par # muc, tau, mu, sigma, xi, beta
    },
    error = function(e) {'errore'})
 #   if(class(condizione) == "character") j = j + 1
  return(MLE)
}

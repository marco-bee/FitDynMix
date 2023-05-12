#' Log-likelihood of a Lognormal-GPD dynamic mixture
#'
#' This function evaluates the log-likelihood of a Lognormal-GPD dynamic mixture,
#' approximating the normalizing constant via Monte Carlo simulation.
#' @param x (6 by 1) numerical vector: values of the parameters CA1, CA2, meanlog, sdlog, xi, beta.
#' @param y vector: points where the function is evaluated.
#' @param nreps non-negative integer: number of replications to be used in the computation of the integral in the normalizing
#' constant.
#' @return Log-likelihood of the lognormal-GPD mixture evaluated at y.
#' @keywords dynamic mixture.
#' @export
#' @examples
#' llik <- dynloglikMC(c(1,2,0,1,.25,3.5),Metro2019,10000)


dynloglikMC = function (x, y, nreps)
{
  muc <- x[1]
  tau <- x[2]
  mu <- x[3]
  sigma <- x[4]
  xi <- x[5]
  beta <- x[6]
  f <- function(x) (evir::dgpd(x, xi, 0, beta) - dlnorm(x, mu, sigma)) * atan((x - muc)/tau)
  p <- pcauchy(y, muc, tau)
  temp <- (1 - p) * dlnorm(y, mu, sigma) + p * evir::dgpd(y, xi, 0, beta)
  Z = nConst_MC(x, nreps)
  llik <- sum(log(temp/Z))
  # if (Z<0 | min(temp)<0)
  # {
  #   print('Shit!')
  #   print(c(muc,tau,mu,sigma,xi,beta))
  #   print(c(min(temp),max(temp)))
  # }
  return(llik)
}

#' Log-likelihood of a Lognormal-GPD dynamic mixture
#'
#' This function evaluates the log-likelihood of a Lognormal-GPD dynamic mixture,
#' with Cauchy or exponential weight, approximating the normalizing constant via Monte Carlo simulation.
#' @param x (6 by 1) numerical vector: values of the parameters CA1, CA2, meanlog, sdlog, xi, beta.
#' @param y vector: points where the function is evaluated.
#' @param nreps non-negative integer: number of replications to be used in the computation of the integral in the normalizing
#' constant.
#' @param xiInst non-negative real: shape parameter of the instrumental GPD.
#' @param betaInst non-negative real: scale parameter of the instrumental GPD.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return Log-likelihood of the lognormal-GPD mixture evaluated at y.
#' @keywords dynamic mixture.
#' @export
#' @examples
#' llik <- dynloglikMC(c(1,2,0,1,.25,3.5),Metro2019,10000,3,3,'exp')

dynloglikMC = function (x, y, nreps, xiInst, betaInst, weight)
{
  if (weight == 'cau')
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
    Z = nConst_MC(x, nreps, xiInst, betaInst, 'cau')
    llik <- sum(log(temp/Z))
  }
  
  if (weight == 'exp')
  {
    lambda <- x[1]
    mu <- x[2]
    sigma <- x[3]
    xi <- x[4]
    beta <- x[5]
    f <- function(x) (dlnorm(x,mu,sigma) - evir::dgpd(x, xi, mu=0, beta)) *
      exp(-lambda*x)
    p <- pexp(y,lambda)  
    temp <- (1 - p) * dlnorm(y, mu, sigma) + p * evir::dgpd(y, xi, 0, beta)
    Z = nConst_MC(x, nreps, xiInst, betaInst, 'exp')
    llik <- sum(log(temp/Z))
  }
  return(llik)
}

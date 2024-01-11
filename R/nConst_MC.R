#' Monte carlo approximation of the normalizing constant of a Lognormal-GPD dynamic mixture
#'
#' This function evaluates via Monte Carlo simulation the normalizing constant of a Lognormal-GPD dynamic mixture,
#' with Cauchy or exponential weight.
#' @param x if weight is equal to 'cau', (6 by 1) numerical vector: values of \eqn{\mu_c}, \eqn{\tau}, \eqn{\mu}, \eqn{\sigma}, \eqn{\xi}, \eqn{\beta};
#' if weight is equal to 'exp', (5 by 1) numerical vector: values of \eqn{\lambda}, \eqn{\mu}, \eqn{\sigma}, \eqn{\xi}, \eqn{\beta}.
#' @param nreps  number of replications to be used in the computation of the integral in the normalizing
#' constant.
#' @param xiInst real: value of the shape parameter of the instrumental density. Default value equal to 3.
#' @param betaInst non-negative real: value of the scale parameter of the instrumental density. Default value equal to 3.
#' @param weight 'cau' or 'exp': name of weight distribution.
#' @return Normalizing constant of the density of the lognormal-GPD mixture.
#' @keywords dynamic mixture.
#' @export
#' @examples
#' nconst <- nConst_MC(c(1,2,0,0.5,0.25,3.5), 10000, 3, 3, 'cau')

nConst_MC = function (x, nreps, xiInst=3, betaInst=3, weight)
{
  if (weight =='cau')
  {
    muc <- x[1]
    tau <- x[2]
    mu <- x[3]
    sigma <- x[4]
    xi <- x[5]
    beta <- x[6]
    f <- function(x) (evir::dgpd(x, xi, 0, beta) - dlnorm(x, mu, sigma)) * atan((x - muc)/tau)
    xsim = evir::rgpd(nreps,xiInst,0,betaInst)
    num = f(xsim)
    den = evir::dgpd(xsim,xiInst,0,betaInst)
    I = mean(num/den)
    Z <- 1 + (1/pi) * I
  }
  if (weight =='exp')
  {
    lambda <- x[1]
    mu <- x[2]
    sigma <- x[3]
    xi <- x[4]
    beta <- x[5]
    f <- function(x) (dlnorm(x,mu,sigma) - evir::dgpd(x, xi, mu=0, beta)) *
      exp(-lambda*x)
    xsim = evir::rgpd(nreps,xiInst,0,betaInst)
    num = f(xsim)
    den = evir::dgpd(xsim,xiInst,0,betaInst)
    I = mean(num/den)
    Z <- 1 + I
  }
  return(Z)
}

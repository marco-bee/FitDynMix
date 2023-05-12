#' Density of a Lognormal-GPD dynamic mixture
#'
#' This function evaluates the normalizing constant of a Lognormal-GPD dynamic mixture via Monte Carlo simulation.
#' @param x (6 by 1) numerical vector: values of the parameters CA1, CA2, meanlog, sdlog, xi, beta.
#' @param nreps  number of replications to be used in the computation of the integral in the normalizing
#' constant.
#' @param xiInst real: value of the shape parameter of the instrumental density. Default value equal to 3.
#' @param betaInst non-negative real: value of the scale parameter of the instrumental density. Default value equal to 3.
#' @return Normalizing constant of the density of the lognormal-GPD mixture.
#' @keywords dynamic mixture.
#' @export
#' @examples
#' nconst <- nConst_MC(x, 10000)

nConst_MC = function (x, nreps, xiInst=3, betaInst=3)
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
  return(Z)
}

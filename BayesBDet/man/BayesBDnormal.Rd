\name{BayesBDnormal}
\alias{BayesBDnormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian boundary estimation for continuous intensity images
}
\description{
Estimate boundaries in a continuous intensity image.
}
\usage{
BayesBDnormal(obs, inimean, nrun, nburn, J, slice)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
	The noisy observation which is a list with the following required elements:
	\itemize{
	\item{intensity: }{observed intensity at each pixel.}
	\item{theta.obs, r.obs: }{the location of the pixel at which the intensity is observed, using polar coordinates with respect to a reference point.}
	\item{center: }{the reference point for polar coords (theta.obs, r.obs).}}
}
  \item{inimean}{
	a constant to specify the initial mean functions in the Bayesian estimation.
}
  \item{nrun}{
	the number of MCMC samples to keep for estimation.
}
  \item{nburn}{
	the number of initial MCMC samples to discard.
}
  \item{J}{
	truncation number of the Gaussian process kernel. The number of eigenfunctions is \eqn{2J + 1}.
}
}
\value{
	\item{estimate }{Posterior mean estimate of image boundary at theta values.}
	\item{theta }{A grid of 200 values on \eqn{[0,2\pi]} at which to retrun the estimated boundary.}
	\item{lower, upper}{The lower and upper bounds of a \eqn{95\%} uniform credible band for the image boundary.}
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}
\author{

}
\examples{
par(mfrow = c(1, 2))
gamma.fun = ellipse(a = 0.35, b = 0.25)
obs = parnormobs(m = 100, mu.in = 4, mu.out = 1,sd.in = 1.5, sd.out = 1, design = 'J', gamma.fun)
rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
image(rotate(obs$intensity), axes = FALSE, asp = 1, main = 'observation')
cppsamp = BayesBDnormal(obs, .4, 4000, 1000, 10, TRUE)
x = cppsamp$estimate*cos(cppsamp$theta)+obs$center[1]
y = cppsamp$estimate*sin(cppsamp$theta) +.5obs$center[2]
lines(x,y, lty = 2, lwd = 3)
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
require(plotrix)
my.radial = function(r, theta, ...){
radial.plot(c(r[order(theta)]), c(theta[order(theta)]),
rp.type = "p", show.grid.label = TRUE, radial.lim = c(0, 0.5),
...)
}
my.radial(gamma.fun(theta.plot), theta.plot, line.col = 1, lty = 1, lwd = 1, show.grid = FALSE)
x = cppsamp$upper*cos(cppsamp$theta)
y = cppsamp$upper*sin(cppsamp$theta)
polygon(x, y,fillOddEven = TRUE, col = 'gray', border = NA)
x = cppsamp$lower*cos(cppsamp$theta)
y = cppsamp$lower*sin(cppsamp$theta) 
polygon(x, y,fillOddEven = TRUE, col = 'white', border = NA)
x = cppsamp$estimate*cos(cppsamp$theta)
y = cppsamp$estimate*sin(cppsamp$theta)
lines(x,y, lty = 2, lwd = 3)
}


\name{BayesBDet-package}
\alias{BayesBDet-package}
\alias{BayesBDet}
\docType{package}
\title{
  Bayesian Inference for Image Boundaries
}
\description{
  This package provides tools for carrying out a Bayesian analysis of image boundaries.  Functions are provided
  for both binary (Bernoulli) and continuous (Gaussian) images.  Examples, along with an interactive shiny function
  illustrate how to perform simulations, analyze custom data, and plot estimates and credible intervals.    
}
\details{
Package: BayesBDet
License: GPL (>= 3) 
Imports: Rcpp (>= 0.12.5), RcppArmadillo, shiny, plotrix
LinkingTo: Rcpp, RcppArmadillo
}
\author{
Nicholas Syring

Maintainer: Nicholas Syring <nasyrin@gmail.com>
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}
\examples{
  \dontrun{
## Analyzing a continuous image using BayesBDnormal
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
\dontrun{
## Analyzing a binary image using BayesBDbinary
set.seed(12345)
par(mfrow = c(1, 2))
gamma.fun = ellipse(a = 0.35, b = 0.25)
obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2, design = 'J', gamma.fun)
rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
image(rotate(obs$intensity), axes = FALSE, asp = 1, main = 'observation')
cppsamp = BayesBDbinary(obs, .4, 4000, 1000, 10, TRUE)
x = cppsamp$estimate*cos(cppsamp$theta)+obs$center[1]
y = cppsamp$estimate*sin(cppsamp$theta)+obs$center[2]
lines(x,y, lty = 2, lwd = 3)
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
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
\dontrun{
## Launch a shiny session for easy use of the BayesBDet package
BayesBDshiny()
}
}

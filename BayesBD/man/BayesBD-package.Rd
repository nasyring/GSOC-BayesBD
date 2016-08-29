\name{BayesBD-package}
\alias{BayesBD-package}
\alias{BayesBD}
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
Package: BayesBD
License: GPL (>= 3) 
Imports: Rcpp (>= 0.12.5), RcppArmadillo, shiny, plotrix
LinkingTo: Rcpp, RcppArmadillo
}
\author{
Nicholas Syring, Meng Li

Maintainer: Nicholas Syring <nasyrin@gmail.com>
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}
\examples{
  \dontrun{
## Analyzing a continuous image using BayesBDnormal
set.seed(12345)
gamma.fun = ellipse(a = 0.35, b = 0.25)
norm.obs = parnormobs(m = 100, mu.in = 4, mu.out = 1, sd.in = 1.5, sd.out = 1, design = 'J', center = c(0.5,0.5), gamma.fun)
norm.samp = BayesBDnormal(norm.obs, .4, 4000, 1000, 10, TRUE, FALSE)
plotBD(norm.obs, norm.samp)
\dontrun{
## Analyzing a binary image using BayesBDbinary
set.seed(12345)
gamma.fun = ellipse(a = 0.35, b = 0.25)
bin.obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2, design = 'J', center = c(0.5,0.5), gamma.fun)
bin.samp = BayesBDbinary(bin.obs, .4, 4000, 1000, 10, TRUE, FALSE)
plotBD(bin.obs, bin.samp)
\dontrun{
## Launch a shiny session for easy use of the BayesBD package
BayesBDshiny()
}
}

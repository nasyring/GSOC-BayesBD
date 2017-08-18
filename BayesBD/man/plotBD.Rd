\name{plotBD}
\alias{plotBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of posterior boundary estimates and data
}
\description{
Produces plots of image data, posterior boundary estimates, \cr
and 95% uniform credible bands.
}
\usage{
plotBD(fitted.image, plot.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted.image}{An object containing the \cr output of either fitBinImage or fitContImage.}
  \item{plot.type}{takes values 1, 2 or 3:\cr 1 to plot the data only; \cr 
	2 to plot the estimate and 95\% uniform credible bands; \cr and 3 
	to plot the data with the estimated boundary overlaid.}
}
\value{
	There is no output to console; plotBD produces plots.
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}
\examples{
\dontrun{
set.seed(12345)
gamma.fun = ellipse(a = 0.35, b = 0.25)
bin.obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2,
 design = 'J', center = c(0.5,0.5), gamma.fun)
bin.fit = fitBinImage(image = bin.obs, nrun=1000, nburn=1000,
 J=10, ordering='I', slice = FALSE, outputAll=FALSE)
par(mfrow = c(1,3))
plotBD(bin.fit, 1)
plotBD(bin.fit, 2)
plotBD(bin.fit, 3)
}
}


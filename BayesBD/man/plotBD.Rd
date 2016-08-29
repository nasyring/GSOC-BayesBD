\name{plotBD}
\alias{plotBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of posterior boundary estimates and data
}
\description{
Produces plots of image data, posterior boundary estimates, and 95% uniform credible bands.
}
\usage{
plotBD(obs, samples)
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
  \item{samples}{
	output from either BayesBDbinary or BayesBDnormal.
}
}
\value{
	There is no output to console. plotBD produces two plots in one window, on the left is the image data, and on the right is the posterior mean boundary, 95% uniform credible bands, and true boundary if available.
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}
\examples{
\dontrun{
set.seed(12345)
gamma.fun = ellipse(a = 0.35, b = 0.25)
bin.obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2, design = 'J', center = c(0.5,0.5), gamma.fun)
bin.samp = BayesBDbinary(bin.obs, .4, 1000, 0, 10, TRUE, FALSE)
plotBD(bin.obs, bin.samp)
}
}



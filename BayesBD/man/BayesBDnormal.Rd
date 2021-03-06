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
BayesBDnormal(obs, inimean, nrun, nburn, J, ordering_mu, 
ordering_sigma, mask, slice, outputAll)
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
  \item{ordering_mu}{
	Indicates which Gaussian distribution has larger mean intensity: "I", the Gaussian distribution inside the boundary; "O", the Gaussian distribution outside the boundary; "N", no ordering information is available. 
}
  \item{ordering_sigma}{
	Indicates which Gaussian distribution has larger intensity variance: "I", the Gaussian distribution inside the boundary; "O", the Gaussian distribution outside the boundary; "N", no ordering information is available. 
}
  \item{mask}{
	Logical vector (same length as obs$intensity) to indicate region of interest.  Should this data point be included in the analysis?
}
\item{slice}{
	boolean where TRUE means that slice sampling will be used to sample Fourier basis function coefficients and FALSE means that Metropolis-Hastings will be used instead.
}
  \item{outputAll}{
	boolean controlling the amount of output produced, see value below.
}
}
\value{
If outputAll is FALSE,
	\item{estimate }{Posterior mean estimate of image boundary at theta values.}
	\item{theta }{A grid of 200 values on \eqn{[0,2\pi]} at which to retrun the estimated boundary.}
	\item{lower, upper}{The lower and upper bounds of a \eqn{95\%} uniform credible band for the image boundary.}
If outputAll is TRUE, same as above, and additionally,
	\item{musig.smp}{posterior samples of \eqn{\mu_1}, \eqn{\mu_2}, \eqn{\sigma_1}, and \eqn{\sigma_2}.}
	\item{coef.smp}{posterior samples of Fourier basis function coefficients.}
}
\seealso{
   \code{\link{fitContImage}}
}
\references{
Li, M. and Ghosal, S.(2015) "Bayesian Detection of Image Boundaries." arXiv 1508.05847.
}


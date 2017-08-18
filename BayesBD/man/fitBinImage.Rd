\name{fitBinImage}
\alias{fitBinImage}

\title{
Data pre-processing and binary image analysis
}
\description{
This function can be used to analyze a binary image in .png or .jpeg,
 or an image represented as a list object in the format of par2obs.
}
\usage{
fitBinImage(image, gamma.fun = NULL, center = NULL, inimean = NULL, nrun,
 nburn, J, ordering, mask = NULL, slice, outputAll)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
	This may be a string representing the path to a .png or .jpeg file, or a list object in the same format as par2obs output, with intensity, r.obs, theta.obs, and center the required list contents.
}
  \item{gamma.fun}{
	This is a function, like triangle2 or ellipse, denoting the true boundary.  It is optional and only used when the image input refers to a .png or .jpeg file.
}
  \item{center}{
	This is required if the image input refers to a .png or .jpeg file, otherwise it is unused.
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
  \item{ordering}{
	Indicates which Bernoulli distribution has larger success probability: "I", the Bernoulli distribution inside the boundary; "O", the Bernoulli distribution outside the boundary; "N", no ordering information is available. 
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
\item{output}{}
If outputAll is FALSE,
	\item{estimate }{Posterior mean estimate of image boundary at theta values.}
	\item{theta }{A grid of 200 values on \eqn{[0,2\pi]} at which to retrun the estimated boundary.}
	\item{lower, upper}{The lower and upper bounds of a \eqn{95\%} uniform credible band for the image boundary.}
If outputAll is TRUE, same as above, and additionally,
	\item{pi.smp}{posterior samples of \eqn{\pi_1} and \eqn{\pi_2}.}
	\item{coef.smp}{posterior samples of Fourier basis function coefficients.}
\item{image}{the input image passed to fitBinImage.}
\item{obs}{the processed image data passed to BayesBDbinary.}
}
\seealso{
   \code{\link{par2obs}}
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
bin.fit = fitBinImage(image = bin.obs, nrun=1000,
 nburn=1000, J=10, ordering='I', slice = FALSE, outputAll=FALSE)
par(mfrow = c(1,3))
plotBD(bin.fit, 1)
plotBD(bin.fit, 2)
plotBD(bin.fit, 3)
}
}


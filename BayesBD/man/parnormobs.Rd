\name{parnormobs}
\alias{parnormobs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Gaussian intensity images
}
\description{
The generated data \eqn{Y ~ N(\mu, \sigma)} is the image intensity, where the params \eqn{(\mu, \sigma)} are determined by whether the location \eqn{X}, given by polar coords \eqn{(r,\theta)}, is inside the boundary or not. 
}
\usage{
parnormobs(m, mu.in, mu.out, sd.in, sd.out, design, gamma.fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
	\eqn{m * m} observations will be generated over the unit square centered at (0,0).
}
  \item{mu.in}{
	The mean intensity for pixels inside the image boundary.
}
  \item{mu.out}{
	The mean intensity for pixels outside the image boundary.
}
  \item{sd.in}{
	The standard deviation of intensity for pixels inside the image boundary.
}
  \item{sd.out}{
	The standard deviation of intensity for pixels outside the image boundary.
}
  \item{design}{
	Taking values: 'D' for deterministic (equally-spaced grid) design or 'J' for jitteredly random design.
}
  \item{gamma.fun}{
	The function to generate boundaries, see ellipse or triangle2.
}
}
\examples{
set.seed(2015)
# use ellipse boundary
gamma.fun = ellipse(a = 0.35, b = 0.25)
obs = parnormobs(m = 100, mu.in = 1, sd.in = 1, mu.out = 0, sd.out = 1, design = 'J', gamma.fun)
# rotate a matrix
rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
# rotate & image it - square (asp = 1)
image(rotate(obs$intensity), axes = FALSE, asp = 1, main = 'observation')
}


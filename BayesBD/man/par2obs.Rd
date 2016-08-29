\name{par2obs}
\alias{par2obs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate binary intensity images 
}
\description{
The generated data \eqn{Y ~ Bernoulli(p)} is the image intensity, where the success probability \eqn{p} is determined by whether the location \eqn{X}, given by polar coords \eqn{(r,\theta)}, is inside the boundary or not. 
}
\usage{
par2obs(m, pi.in, pi.out, design, gamma.fun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
	\eqn{m * m} observations will be generated over the unit square centered at (0,0).
}
  \item{pi.in}{
	The success probability, \eqn{P(Y_i = 1)}, where \eqn{Y_i} is intensity of pixel \eqn{i} if the location is inside the boundary.
}
  \item{pi.out}{
	The success probability, \eqn{P(Y_i = 1)}, where \eqn{Y_i} is intensity of pixel \eqn{i} if the location is outside the boundary.
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
obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2, design = 'J', gamma.fun)
# rotate a matrix
rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
# rotate & image it - square (asp = 1)
image(rotate(obs$intensity), axes = FALSE, asp = 1, main = 'observation')
}


\name{ellipse}
\alias{ellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate ellipse boundaries
}
\description{
Generate general ellipse boundaries with semidiameters a and b, centered at (r0, theta0), with the a axis rotated by phi relative to the polar axis. 
}
\usage{
ellipse(a, b, r0 = 0, theta0 = 0, phi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, b}{
	semidiameter parameters where a > b.
}
  \item{r0, theta0}{
	the polar form of the center of the ellipse
}
  \item{phi}{
	the angle rotated by the major axis.
}
}

\examples{
\dontrun{
gamma.fun = ellipse(a = 0.35, b = 0.25)
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
x = gamma.fun(theta.plot)*cos(theta.plot)
y = gamma.fun(theta.plot)*sin(theta.plot)
plot(x,y,type = 'l', axes=TRUE, frame.plot=FALSE)
}
}

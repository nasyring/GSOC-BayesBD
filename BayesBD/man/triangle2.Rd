\name{triangle2}
\alias{triangle2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate triangle boundaries
}
\description{
Generate boundaries of equilateral triangles of height S.
}
\usage{
triangle2(S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
	height of the generated triangle.
}
}

\examples{
gamma.fun = triangle(0.5)
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
require(plotrix)
my.radial <- function(r, theta, ...){
  radial.plot(c(r[order(theta)]), c(theta[order(theta)]),
              rp.type = "p", show.grid.label = TRUE, radial.lim = c(0, 0.5),
              ...)
}
my.radial(gamma.fun(theta.plot), theta.plot, main = 'Ellipse boundary',
          line.col = 1, lty = 2, lwd = 2, show.grid = FALSE)
}



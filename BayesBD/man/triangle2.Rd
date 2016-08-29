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
\dontrun{
gamma.fun = triangle(0.5)
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
x = gamma.fun(theta.plot)*cos(theta.plot)
y = gamma.fun(theta.plot)*sin(theta.plot)
plot(x,y,type = 'l', axes=TRUE, frame.plot=FALSE)
}
}



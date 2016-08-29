\name{eigenfun}
\alias{eigenfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the Fourier series
}
\description{
eigenfun(n, x) generates the first n fourier series, which are \eqn{\{1/\sqrt{2*\pi}, 1/\sqrt{\pi}*cos(x), 1/\sqrt{\pi}*cos(x), ...\}} evaluated at x.
}
\usage{
eigenfun(n,x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
	The number of eigenfunctions used to approximate the boundary.  It should be a odd number.
}
\item{x}{
	The argument at which the eigenfunctions are evaluated. 
}
}


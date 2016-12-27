\name{besselIs}
\alias{besselIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bessel Function
}
\description{
Exponentially scaled Bessel function of first kind.
}
\usage{
besselIs(x, nu, expon.scaled = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	numeric, \eqn{\geq 0}.

}
  \item{nu}{
	numeric; The order (may be fractional) of the Bessel function.
}
  \item{expon.scaled}{
	This argument does not affect the result; besselIs always returns the exponentially scaled version of the Bessel function of the first kind.
}
}
\seealso{ \code{\link[base]{besselI}} }



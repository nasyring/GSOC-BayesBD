\name{unisliceL}
\alias{unisliceL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Slice Sampler
}
\description{
Performs slice sampling for lambda parameter. This is a helper function used by BayesBDbinary and BayesBDnormal.
}
\usage{
unisliceL(x0, gx0, i_J, tauini, anini, alphalambda, betalambda, lambdaini)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
	initial value of lambda.
}
  \item{gx0}{
	likelihood at x0.
}
  \item{i_J}{
	number of basis functions (2*i_J+1).
}
  \item{tauini}{
	tau parameter in covariance function.
}
  \item{anini}{
	2J+1 vector of eigenfunction coefficients.
}
  \item{alphalambda}{
	hyperparameter for sampling tau, usually 1.
}
  \item{betalambda}{
	hyperparameter for sampling tau, usually 1.
}
  \item{lambdaini}{
	current value of lambda.
}
}


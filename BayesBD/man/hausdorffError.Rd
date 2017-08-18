\name{hausdorffError}
\alias{hausdorffError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hausdorff Error of boundary estimate 
}
\description{
Utility function to calculate the Hausdorff error between true image boudnary and estimated image boundary in a simulated image.
}
\usage{
hausdorffError(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	the output of fitBinImage or fitContImage
}
}
\examples{
set.seed(12345)
gamma.fun = ellipse(a = 0.35, b = 0.25)
norm.obs = parnormobs(m = 100, mu.in = 4, mu.out = 1, sd.in = 1.5, sd.out = 1, design = 'J',
 center = c(0.5,0.5), gamma.fun)
\dontrun{
norm.samp = fitContImage(norm.obs, NULL,NULL,.4, 1000, 0, 10,"I","I",rep(1,10000), FALSE, FALSE)
par(mfrow = c(1,3))
plotBD(norm.samp, 1)
plotBD(norm.samp, 2)
plotBD(norm.samp, 3)
hausdorffError(norm.samp)
}
}


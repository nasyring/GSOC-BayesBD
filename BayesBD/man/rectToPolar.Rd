\name{rectToPolar}
\alias{rectToPolar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate binary intensity images 
}
\description{
Utility function to convert rectangular (Euclidean) coordinates to polar coordinates.
}
\usage{
rectToPolar(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	x-axis in rectangular coords
}
  \item{y}{
	y-axis in rectangular coords	
}
}
\examples{
x = runif(100,-1,1)
y = runif(100,-1,1)
polar_coords = rectToPolar(x,y)
}


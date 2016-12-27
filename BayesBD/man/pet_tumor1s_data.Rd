\name{pet_tumor1s_data}
\alias{pet_tumor1s_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Positron Emission Tomography data of a tumor.  
}
\description{
This is an example data set.  Use the script below for analysis.
}

\examples{
\dontrun{
# Load the data.  Check the path before running.  
# The file is in the data folder of your BayesBD package.
tumor1 = read.csv("pet_tumor1s_data.csv", header=FALSE) 
# Run code down to the plot function
z = as.matrix(tumor1)
x = z
y = z
nc = ncol(z)
nr = nrow(z)
for(i in 1:nr){
	for(j in 1:nc){
		x[i,j] = i/nr - 0.5/nr
		y[i,j] = j/nc - 0.5/nc
	}
}
intensity = as.vector(z)
x = as.vector(x)
y1 = as.vector(y)
y = -y1+1
center = c(0.7,0.5)
r.obs = sqrt((x-.7)^2+(y-.5)^2)
theta.obs = atan2((y-.5),(x-.7))
subset = ifelse(intensity!=1,1,0)
tumor1.data = list(intensity = intensity, 
theta.obs = theta.obs,r.obs=r.obs,center=center,x=x,y=y)
plotBD(tumor1.data)
tumor1.samp = BayesBDnormal(tumor1.data, .01, 4000, 1000, 10,"N","N",subset, FALSE, FALSE)
# Plot the results for the first tumor.
plotBD(tumor1.data, tumor1.samp)
}
}
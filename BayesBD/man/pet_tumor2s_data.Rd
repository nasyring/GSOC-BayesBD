\name{pet_tumor2s_data}
\alias{pet_tumor2s_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Positron Emission Tomography data of a tumor.  
}
\description{
This is an example data set.  Use the script below for analysis.
}
\examples{
\dontrun{
# Load the data.  Be sure to check the path.
tumor2 = read.csv("pet_tumor2s_data.csv", header=FALSE) 
# Run code down to the plot function
z = as.matrix(tumor2)
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
center = c(0.4,0.55)
r.obs = sqrt((x-.4)^2+(y-.55)^2)
theta.obs = atan2((y-.55),(x-.4))
theta.obs = ifelse(theta.obs<0,theta.obs+2*pi,theta.obs)
subset = ifelse(intensity!=1,1,0)
tumor2.data = list(intensity = intensity, 
theta.obs = theta.obs,r.obs=r.obs,center=center,x=x,y=y)
plotBD(tumor2.data)
tumor2.samp = BayesBDnormal(tumor2.data, .1, 4000, 1000,
 10,"N","N",subset, FALSE, FALSE)
# Plot the results for the second tumor
plotBD(tumor2.data, tumor2.samp)
# Now, we will analyze the data from the second tumor 
# within the estimated boundary.
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
# Run code down to the plot function.
r.est = function(theta){
	thetas = c(theta.plot,2*pi)
	r.thetas = c(tumor2.samp$upper,tumor2.samp$upper[1])
	s = sort(c(theta,thetas))
	w = which(s==theta)
	lt = s[w-1]
	ut = s[w+1]
	lr = r.thetas[w-1]
	ur = r.thetas[w]
	r_est = ((theta - lt)/(ut-lt))*ur+((ut - theta)/(ut-lt))*lr
	return(r_est[1])
}
app.r.est = function(theta) apply(matrix(theta,length(theta),1),1,r.est)
r_ests = app.r.est(theta.obs)
r_ests = matrix(r_ests,length(theta.obs),1)
subset = ifelse(r.obs<=r_ests,1,0)
tumor2.samp = BayesBDnormal(tumor2.data, .1, 4000, 1000, 10,
"N","N",subset, FALSE, FALSE)
# Plot the results of the second run
plotBD(tumor2.data, tumor2.samp)
}
}
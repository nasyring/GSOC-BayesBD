lebesgueError = function(fit){
	gamma.fun = fit$gamma.fun
	est.boundary = function(theta){
		theta.seq = fit$output$theta
		post.est = fit$output$estimate
		s=sort(c(theta.seq,2*pi,theta))
		index = which(s==theta)[1]
		if(theta==0 || theta==2*pi){ 
			value = post.est[1]
		} else {
			a = theta.seq[index-1]
			b = ifelse(index<201,theta.seq[index],2*pi)
			va = post.est[index-1]
			vb = ifelse(index<201,post.est[index],post.est[1])
		}
		return(((theta-a)/(b-a))*vb+(1-((theta-a)/(b-a)))*va)
	}
	app.est.boundary = function(u){
		u = matrix(u,length(u),1)
		apply(u,1,est.boundary)
	}
	integrand = function(theta){.5*(abs(app.est.boundary(theta)-gamma.fun(theta))*(app.est.boundary(theta)+gamma.fun(theta)))}
	return(integrate(integrand, lower = 0, upper = 2*pi)$value)
}
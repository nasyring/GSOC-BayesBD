hausdorffError = function(fit){
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
	hausdorff.dist = function(x.points, y.points){
		l =nrow(x.points)
		d1 = rep(0,l)
		d2 = rep(0,l)
		for(i in 1:l){
			d1[i] = min(sqrt((x.points[i,1]-y.points[,1])^2+(x.points[i,2]-y.points[,2])^2))
			d2[i] = min(sqrt((x.points[i,1]-y.points[,1])^2+(x.points[i,2]-y.points[,2])^2))
		}
		sup.x.inf.y = max(d1)
		sup.y.inf.x = max(d2)
		return(max(sup.x.inf.y, sup.y.inf.x))
	}
	intensity = fit$obs$intensity
	gamma.fun.radii = gamma.fun(fit$obs$theta.obs)
	ground.truth = ifelse(fit$obs$r.obs<gamma.fun.radii,1,0)
	predicted = app.est.boundary(fit$obs$theta.obs)	
	x.gamma = gamma.fun.radii*cos(fit$obs$theta.obs)
	y.gamma = gamma.fun.radii*sin(fit$obs$theta.obs)
	x.est = predicted*cos(fit$obs$theta.obs)
	y.est = predicted*sin(fit$obs$theta.obs)
	hausdorff.est = hausdorff.dist(cbind(x.gamma, y.gamma), cbind(x.est, y.est))
	return(hausdorff.est)
}
fitContImage = function(image, gamma.fun = NULL, center = NULL, inimean=NULL, nrun, nburn, J, ordering_mu, ordering_sigma, mask = NULL, slice, outputAll){

	is.string <- function(input) {
    		is.character(input) & length(input) == 1
	}
	if(is.string(image)){
		if(any(substr(image,start = nchar(image)-2,stop = nchar(image))=="jpg",substr(image,start = nchar(image)-3,stop = nchar(image))=="jpeg")){
			img = readJPEG(image)
		}else if(substr(image,start = nchar(image)-2,stop = nchar(image))=="png"){
			img = readPNG(image)
		}else {
			return('The image file must have a .png or .jpeg/.jpg file type.')
		}
		
		if(length(dim(img))>2){
			img = matrix(img[,,1],dim(img)[1],dim(img)[2])
		}else {
			img = img
		}

		n1 = nrow(img)
		n2 = ncol(img)

		img_flip = img
	
		for(i in 1:n1){
			for(j in 1:n2){
				img_flip[i,j] = img[n1-i+1,j]
			}
		}

		img = img_flip

		if(any(center[1]>n2,center[2]>n1,center[1]<0,center[2]<0)){
			return(paste('The center should be a pixel (x,y) between 0 and ', ncol(img), 'for x, and 0 and ', nrow(img), ' for y.'))
		}

		r.obs = img
		theta.obs = img
		intensity = img
		for(i in 1:n1){
			for(j in 1:n2){
				r.obs[i,j] = sqrt(((i - center[2])/n1)^2 + ((j - center[1])/n2)^2)
				theta.obs[i,j] = atan2((i - center[2])/n1, (j - center[1])/n2)
				theta.obs[theta.obs < 0] = theta.obs[theta.obs < 0] + 2*pi
			}
		}

		r.obs = as.vector(r.obs)
		theta.obs = as.vector(theta.obs)
		intensity = as.vector(intensity)
		intensity = 200*intensity	
	
		obs = list(r.obs = r.obs, theta.obs = theta.obs, intensity = intensity, center=center)

		if(is.null(mask)){
			mask = rep(1,length(obs$intensity))
		}else {
			center = obs$center
			new.r.obs = obs$r.obs[mask==1]
			new.theta.obs = obs$theta.obs[mask==1]
			new.intensity = obs$intensity[mask==1]
			obs = list(r.obs=new.r.obs, theta.obs = new.theta.obs, intensity=new.intensity, center=center)
		}

if(is.null(inimean)){
			ini.mean.estimator = function(r){
				obs.in = obs$intensity[which(obs$r.obs<r)]
				obs.out = obs$intensity[which(obs$r.obs>=r)]
				mu.in = mean(obs.in)
				mu.out = mean(obs.out)
				sd.in = sd(obs.in)
				sd.out = sd(obs.out)
				log_lhood = sum(pnorm(obs.in, mu.in, sd.in, log.p=TRUE))+sum(pnorm(obs.out, mu.out, sd.out, log.p=TRUE))
				if(any(ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='I' & sd.in>sd.out, 
					ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='I' & sd.in>sd.out,
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='N',
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='N',
					ordering_mu=='N' & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='N' & ordering_sigma=='I' & sd.in>sd.out,
					ordering_mu=='N' & ordering_sigma=='N')){
				
					return(log_lhood)

				}else{
					return(-Inf)
				}
			}

			eval_seq = seq(from=min(obs$r.obs)+0.02,to=min(.5,max(obs$r.obs)-0.02), by=0.01)
			len=length(eval_seq)
			evals = rep(0,len)
			for(i in 1:len){evals[i]=ini.mean.estimator(eval_seq[i])}
			inimean = eval_seq[which.max(evals)]
		}

		output = BayesBDnormal(obs, inimean, nrun, nburn, J,  ordering_mu, ordering_sigma, mask = rep(1,length(obs$intensity)), slice, outputAll)

		return(list(image = image, center = center, gamma.fun = gamma.fun, output = output, obs = obs))

	}else if(is.list(image)){

		if(is.null(mask)){
			mask = rep(1,length(image$intensity))
		}else {
			center = image$center
			new.r.obs = image$r.obs[mask==1]
			new.theta.obs = image$theta.obs[mask==1]
			new.intensity = image$intensity[mask==1]
			image = list(r.obs=new.r.obs, theta.obs = new.theta.obs, intensity=new.intensity, center=center)
		}

if(is.null(inimean)){
			ini.mean.estimator = function(r){
				obs.in = image$intensity[which(image$r.obs<r)]
				obs.out = image$intensity[which(image$r.obs>=r)]
				mu.in = mean(obs.in)
				mu.out = mean(obs.out)
				sd.in = sd(obs.in)
				sd.out = sd(obs.out)
				log_lhood = sum(pnorm(obs.in, mu.in, sd.in, log.p=TRUE))+sum(pnorm(obs.out, mu.out, sd.out, log.p=TRUE))
				if(any(ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='I' & sd.in>sd.out, 
					ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='I' & sd.in>sd.out,
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='I' & mu.in>mu.out & ordering_sigma=='N',
					ordering_mu=='O' & mu.in<mu.out & ordering_sigma=='N',
					ordering_mu=='N' & ordering_sigma=='O' & sd.in<sd.out,
					ordering_mu=='N' & ordering_sigma=='I' & sd.in>sd.out,
					ordering_mu=='N' & ordering_sigma=='N')){
				
					return(log_lhood)

				}else{
					return(-Inf)
				}
			}

			eval_seq = seq(from=min(image$r.obs)+0.02,to=min(.5,max(image$r.obs)-0.02), by=0.01)
			len=length(eval_seq)
			evals = rep(0,len)
			for(i in 1:len){evals[i]=ini.mean.estimator(eval_seq[i])}
			inimean = eval_seq[which.max(evals)]
		}

		output = BayesBDnormal(image, inimean, nrun, nburn, J,  ordering_mu, ordering_sigma, mask = rep(1,length(image$intensity)), slice, outputAll)
		
		return(list(image = image, center = center, gamma.fun = gamma.fun, output = output, obs = list(r.obs=as.vector(image$r.obs), theta.obs=as.vector(image$theta.obs), intensity=as.vector(image$intensity), center=center)))

	}else {
		return("Input image is not a compatible image file nor a compatible list object.")
	}
}

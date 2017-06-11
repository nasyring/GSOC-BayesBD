fitBinImage = function(image, gamma.fun = NULL, center = NULL, inimean, nrun, nburn, J, ordering, mask = NULL, slice, outputAll){

	is.string <- function(input) {
    		is.character(input) & length(input) == 1
	}
	if(is.string(image)){
		if(substr(image,start = nchar(image)-3,stop = nchar(image))=="jpeg"){
			img = readJPEG(image)
		}else if(substr(image,start = nchar(image)-2,stop = nchar(image))=="png"){
			img = readPNG(image)
		}else {
			return('The image file must have a .png or .jpeg file type.')
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

		if(is.null(mask)){mask = rep(1,length(intensity))}	
	
		obs = list(r.obs = r.obs, theta.obs = theta.obs, intensity = intensity)

		output = BayesBDbinary(obs, inimean, nrun, nburn, J, ordering, mask, slice, outputAll)

		return(list(image = image, center = center, gamma.fun = gamma.fun, output = output, obs = obs))

	}else if(is.list(image)){

		if(is.null(mask)){mask = rep(1,length(image$intensity))}

		output = BayesBDbinary(image, inimean, nrun, nburn, J, ordering, mask, slice, outputAll)
		return(list(image = image, center = center, gamma.fun = gamma.fun, output = output, obs = list(r.obs=as.vector(image$r.obs), theta.obs=as.vector(image$theta.obs), intensity=as.vector(image$intensity))))
	}else {
		return("Input image is not a compatible image file nor a compatible list object.")
	}
}
plotBD =
function (fitted.image, plot.type) 
{

	is.string <- function(input) {
    		is.character(input) & length(input) == 1
	}

	if(is.string(fitted.image$image)){
		if(substr(fitted.image$image,start = nchar(fitted.image$image)-2,stop = nchar(fitted.image$image))=="jpg"){
			img = readJPEG(fitted.image$image)
		}else if(substr(fitted.image$image,start = nchar(fitted.image$image)-2,stop = nchar(fitted.image$image))=="png"){
			img = readPNG(fitted.image$image)
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


		if(any(fitted.image$center[1]>n2,fitted.image$center[2]>n1,fitted.image$center[1]<0,fitted.image$center[2]<0)){
			return(paste('The center should be a pixel (x,y) between 0 and ', ncol(img), 'for x, and 0 and ', nrow(img), ' for y.'))
		}

		y1 = 1:n1
		x1 = 1:n2
		y=NULL
		x=NULL	
		for(i in 1:n2){
			y = c(y,y1)
			x = c(x,rep(i,n1))
		}	

		intensity = img
		intensity = as.vector(intensity)

		estimate.x = fitted.image$output$estimate*cos(fitted.image$output$theta)*n2+fitted.image$center[1]
		estimate.y = fitted.image$output$estimate*sin(fitted.image$output$theta)*n1+fitted.image$center[2]
		upper.x = fitted.image$output$upper*cos(fitted.image$output$theta)*n2+fitted.image$center[1]
		upper.y = fitted.image$output$upper*sin(fitted.image$output$theta)*n1+fitted.image$center[2]
		lower.x = fitted.image$output$lower*cos(fitted.image$output$theta)*n2+fitted.image$center[1]
		lower.y = fitted.image$output$lower*sin(fitted.image$output$theta)*n1+fitted.image$center[2]

		if(plot.type == 1){
			if(min(intensity)<0){ 
       				normalized = (intensity+abs(min(intensity)))/(max(intensity)-min(intensity))
  			}else{
       				normalized = (intensity-abs(min(intensity)))/(max(intensity)-min(intensity))
    			}
    			plot(x, y, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
		}else if(plot.type == 2){
    			plot(x, y, col = 'white', pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)   			
  			polygon(upper.x, upper.y, fillOddEven = TRUE, col = "gray", border = NA)
   			polygon(lower.x, lower.y, fillOddEven = TRUE, col = "white", border = NA)
   			lines(estimate.x, estimate.y, lty = 2, lwd = 3, col='blue')
   			if (!is.null(fitted.image$gamma.fun)) {
     				gamma.x = fitted.image$gamma.fun(fitted.image$output$theta) * cos(fitted.image$output$theta) + fitted.image$center[1]
    				gamma.y = fitted.image$gamma.fun(fitted.image$output$theta) * sin(fitted.image$output$theta) + fitted.image$center[2]
   				lines(gamma.x, gamma.y, lty = 1, lwd = 1)
			}
		}else if(plot.type == 3){
			if(min(intensity)<0){ 
       				normalized = (intensity+abs(min(intensity)))/(max(intensity)-min(intensity))
  			}else{
       				normalized = (intensity-abs(min(intensity)))/(max(intensity)-min(intensity))
    			}
    			plot(x, y, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
   			lines(estimate.x, estimate.y, lty = 2, lwd = 3, col='blue')
			if (!is.null(fitted.image$gamma.fun)) {
     				gamma.x = fitted.image$gamma.fun(fitted.image$output$theta) * cos(fitted.image$output$theta) + fitted.image$center[1]
    				gamma.y = fitted.image$gamma.fun(fitted.image$output$theta) * sin(fitted.image$output$theta) + fitted.image$center[2]
   				lines(gamma.x, gamma.y, lty = 1, lwd = 1)
			}
		}else {
			return("plot.type must be 1, 2, or 3.")
		}
	}else if(is.list(fitted.image$image)){

		x = fitted.image$image$r.obs*cos(fitted.image$image$theta.obs)+fitted.image$image$center[1]
		y = fitted.image$image$r.obs*sin(fitted.image$image$theta.obs)+fitted.image$image$center[2]
		estimate.x = fitted.image$output$estimate*cos(fitted.image$output$theta)+fitted.image$image$center[1]
		estimate.y = fitted.image$output$estimate*sin(fitted.image$output$theta)+fitted.image$image$center[2]
		upper.x = fitted.image$output$upper*cos(fitted.image$output$theta)+fitted.image$image$center[1]
		upper.y = fitted.image$output$upper*sin(fitted.image$output$theta)+fitted.image$image$center[2]
		lower.x = fitted.image$output$lower*cos(fitted.image$output$theta)+fitted.image$image$center[1]
		lower.y = fitted.image$output$lower*sin(fitted.image$output$theta)+fitted.image$image$center[2]

		if(plot.type == 1){

    			if(min(fitted.image$image$intensity)<0){ 
       				normalized = (fitted.image$image$intensity+abs(min(fitted.image$image$intensity)))/(max(fitted.image$image$intensity)-min(fitted.image$image$intensity))
    			}else{
      				normalized = (fitted.image$image$intensity-abs(min(fitted.image$image$intensity)))/(max(fitted.image$image$intensity)-min(fitted.image$image$intensity))
   			}
   			plot(x, y, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)

		}else if(plot.type == 2){
  			plot(x, y, col = 'white', axes = FALSE, xlab = '', ylab = '',asp = 1)   			
  			polygon(upper.x, upper.y, fillOddEven = TRUE, col = "gray", border = NA)
   			polygon(lower.x, lower.y, fillOddEven = TRUE, col = "white", border = NA)
   			lines(estimate.x, estimate.y, lty = 2, lwd = 3, col='blue')
   			if (!is.null(fitted.image$gamma.fun)) {
     				gamma.x = fitted.image$gamma.fun(fitted.image$output$theta) * cos(fitted.image$output$theta) + fitted.image$image$center[1]
    				gamma.y = fitted.image$gamma.fun(fitted.image$output$theta) * sin(fitted.image$output$theta) + fitted.image$image$center[2]
   				lines(gamma.x, gamma.y, lty = 1, lwd = 1)
			}
			
		}else if(plot.type == 3){
    			if(min(fitted.image$image$intensity)<0){ 
       				normalized = (fitted.image$image$intensity+abs(min(fitted.image$image$intensity)))/(max(fitted.image$image$intensity)-min(fitted.image$image$intensity))
    			}else{
      				normalized = (fitted.image$image$intensity-abs(min(fitted.image$image$intensity)))/(max(fitted.image$image$intensity)-min(fitted.image$image$intensity))
   			}
    			plot(x, y, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
   			lines(estimate.x, estimate.y, lty = 2, lwd = 3, col = 'blue')
			if (!is.null(fitted.image$gamma.fun)) {
     				gamma.x = fitted.image$gamma.fun(fitted.image$output$theta) * cos(fitted.image$output$theta) + fitted.image$image$center[1]
    				gamma.y = fitted.image$gamma.fun(fitted.image$output$theta) * sin(fitted.image$output$theta) + fitted.image$image$center[2]
   				lines(gamma.x, gamma.y, lty = 1, lwd = 1)
			}
		}else {
			return("plot.type must be 1, 2, or 3.")
		}

	}else {
		return("Input image is not a compatible image file nor a compatible list object.")
	}
}

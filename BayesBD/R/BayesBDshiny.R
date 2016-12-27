BayesBDshiny = function () 
{
    ui = pageWithSidebar(
	titlePanel("BayesBD"),
      sidebarPanel(selectInput(inputId = "shape", label = "Choose either an elliptical, triangular, or user-supplied boundary, or indicate that the ground truth is unknown.", 
        choices = c("ellipse", "triangle", "file", "unknown")), 
        fileInput(inputId = "shape_file", label = "Use a custom boundary. The file should be an R script of a function called taking\n\t\t as input an angle in [0, 2pi] and returning the radius of the boundary from a reference point."), 
        selectInput(inputId = "data_type", label = "Choose to simulate binary or Gaussian data or input data file below.", 
            choices = c("binary sim", "normal sim", "user binary data", 
                "user normal data")), fileInput(inputId = "data_file", 
            label = "Use image data from file for analysis. The file should be in csv format with no headers. The first column is for intensity, the second gives the radius from the origin, and the third gives the angle from the origin and positive x-axis, and the fourth column gives the subset of data to be analyzed, 1 to include the data point, and 0 to exclude it."), 
	numericInput(inputId = "centerx", label = "Input the X-coordinate and Y-coordinate of the reference point interior to the boundary function.", 
            value = 0.5, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "centery", label = "Y-coordinate of the reference point.", 
            value = 0.5, min = NA, max = NA, step = NA, width = NULL), 
        sliderInput(inputId = "n_burn", label = "Choose a number of posterior samples to burn", 
            value = 1000, min = 500, max = 1000), sliderInput(inputId = "n_run", 
            label = "Choose a number of posterior samples to keep", 
            value = 1000, min = 500, max = 2000), "Use the following inputs for Binary simulations", 
        sliderInput(inputId = "p_in", label = "Choose the Bernoulli success probability inside the image", 
            value = 0.5, min = 0, max = 1), sliderInput(inputId = "p_out", 
            label = "Choose the Bernoulli success probability outside the image", 
            value = 0.2, min = 0, max = 1),
        selectInput(inputId = "ordering", label = "Indicate which region of the image has higher average intensity.", 
            choices = c("Inside", "Outside", "Unknown")),
"Use the following inputs for Gaussian simulations", 
        numericInput(inputId = "mu_in", label = "Mean intensity inside image", 
            value = 1, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "sd_in", label = "Standard deviation inside image", 
            value = 1, min = 0, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "mu_out", label = "Mean intensity outside image", 
            value = 0, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "sd_out", label = "Standard deviation outside image", 
            value = 1, min = 0, max = NA, step = NA, width = NULL), 
        selectInput(inputId = "ordering_mu", label = "Indicate which region of the image has higher average intensity.", 
            choices = c("Inside", "Outside", "Unknown")),
        selectInput(inputId = "ordering_sd", label = "Indicate which region of the image has higher variation in intensity.", 
            choices = c("Inside", "Outside", "Unknown")),
        downloadButton('downloadData', 'Download'),
        actionButton(inputId = "go", label = "Update")), mainPanel( plotOutput("image")))
    

     server = function(input, output) {
        theta.plot = seq(from = 0, to = 2 * pi, length.out = 200)
        rotate <- function(x) t(apply(x, 2, rev))
        data = eventReactive(input$go, {
            center = c(input$centerx, input$centery)
            if (input$shape == "ellipse") {
                gamma.fun = ellipse(a = 0.35, b = 0.25)
            }
            else if (input$shape == "triangle") {
                gamma.fun = triangle2(0.5)
            }
            else if (input$shape == "file") {
                gamma.fun = source(input$shape_file$datapath)$value
            }
            else {
                gamma.fun = FALSE
            }
            if (input$data_type == "binary sim") {
                obs = par2obs(m = 100, pi.in = input$p_in, pi.out = input$p_out, 
                  design = "J", center, gamma.fun)
            }
            else if (input$data_type == "normal sim") {
                obs = parnormobs(m = 100, mu.in = input$mu_in, 
                  mu.out = input$mu_out, sd.in = input$sd_in, 
                  sd.out = input$sd_out, design = "J", center, 
                  gamma.fun)
            }
            else {
                obs.list = read.csv(input$data_file$datapath, header=FALSE)
                obs = list()
                obs$intensity = matrix(obs.list[, 1], length(obs.list[, 
                  1]), 1)
                obs$r.obs = matrix(obs.list[, 2], length(obs.list[, 
                  2]), 1)
                obs$theta.obs = matrix(obs.list[, 3], length(obs.list[, 
                  3]), 1)
                obs$center = center
		obs$sub = matrix(obs.list[, 4], length(obs.list[, 
                  4]), 1)
            }
            if (any(input$data_type == "binary sim", input$data_type == 
                "user binary data")) {
		d_obs = obs
		d_obs$sub = as.vector(d_obs$sub)
		d_obs$intensity = as.vector(d_obs$intensity)
		d_obs$r.obs = as.vector(d_obs$r.obs)
		d_obs$theta.obs = as.vector(d_obs$theta.obs)
		d_obs$intensity = d_obs$intensity[d_obs$sub==1]
		d_obs$r.obs = d_obs$r.obs[d_obs$sub==1]
		d_obs$theta.obs = d_obs$theta.obs[d_obs$sub==1]
		if(input$ordering=="Inside"){
			ordero = "I"
		}else if(input$ordering == "Outside"){
			ordero = "O"
		}else{
			ordero = "N"
		}
		if(input$data_type == "user binary data"){
			use_obs = d_obs
		}else{
			use_obs = obs
		}
                cppsamp = BayesBDbinary(use_obs, 0.1, input$n_run, 
                  input$n_burn, 10,ordero,rep(1,length(use_obs$theta.obs)), FALSE, FALSE)
            }
            else {
		d_obs = obs
		d_obs$sub = as.vector(d_obs$sub)
		d_obs$intensity = as.vector(d_obs$intensity)
		d_obs$r.obs = as.vector(d_obs$r.obs)
		d_obs$theta.obs = as.vector(d_obs$theta.obs)
		d_obs$intensity = d_obs$intensity[d_obs$sub==1]
		d_obs$r.obs = d_obs$r.obs[d_obs$sub==1]
		d_obs$theta.obs = d_obs$theta.obs[d_obs$sub==1]
		if(input$ordering_mu=="Inside"){
			order_mu = "I"
		}else if(input$ordering_mu == "Outside"){
			order_mu = "O"
		}else{
			order_mu = "N"
		}
		if(input$ordering_sd=="Inside"){
			order_sd = "I"
		}else if(input$ordering_sd == "Outside"){
			order_sd = "O"
		}else{
			order_sd = "N"
		}
		if(input$data_type == "user binary data"){
			use_obs = d_obs
		}else{
			use_obs = obs
		}     
            cppsamp = BayesBDnormal(use_obs, 0.1, input$n_run, 
                  input$n_burn, 10,order_mu,order_sd,rep(1,length(use_obs$theta.obs)), FALSE, FALSE)
            }
            x = cppsamp$estimate * cos(cppsamp$theta) + obs$center[1]
            y = cppsamp$estimate * sin(cppsamp$theta) + obs$center[2]
            lx = cppsamp$lower * cos(cppsamp$theta) + obs$center[1]
            ly = cppsamp$lower * sin(cppsamp$theta) + obs$center[2]
            ux = cppsamp$upper * cos(cppsamp$theta) + obs$center[1]
            uy = cppsamp$upper * sin(cppsamp$theta) + obs$center[2]
	    ru = sqrt((ux-center[1])^2 + (uy-center[2])^2)
	    r.est = function(theta){
		thetas = c(theta.plot,2*pi)
		r.thetas = c(ru,ru[1])
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
            r_ests = app.r.est(obs$theta.obs)
            r_ests = matrix(r_ests,length(r_ests),1)
            subset = ifelse(as.vector(obs$r.obs)<=r_ests,1,0)
            return(list(obs = obs, x = x, y = y, lx = lx, ly = ly, 
                ux = ux, uy = uy, gamma.fun = gamma.fun, shape = input$shape, subset = subset))
        })
        output$image = renderPlot({
            d = data()
            par(mfrow = c(1, 2), mar = c(3, 2, 2, 2))
  xobs = d$obs$r.obs*cos(d$obs$theta.obs)+d$obs$center[1]
    yobs = d$obs$r.obs*sin(d$obs$theta.obs)+d$obs$center[2]
    if(min(d$obs$intensity)<0){ 
       normalized = (d$obs$intensity+abs(min(d$obs$intensity)))/(max(d$obs$intensity)-min(d$obs$intensity))
    }else{
       normalized = (d$obs$intensity-abs(min(d$obs$intensity)))/(max(d$obs$intensity)-min(d$obs$intensity))
    }
    plot(xobs, yobs, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
            if (d$shape != "unknown") {
		image(matrix(0,100,100),asp = 1, axes=FALSE,col = 'white')    
		x = d$obs$gamma.fun(theta.plot)*cos(theta.plot)+d$obs$center[1]
		y = d$obs$gamma.fun(theta.plot)*sin(theta.plot)+d$obs$center[2]
		max_x = max(d$obs$r.obs*cos(d$obs$theta.obs)+d$obs$center[1])
		min_x = min(d$obs$r.obs*cos(d$obs$theta.obs)+d$obs$center[1])
		max_y = max(d$obs$r.obs*sin(d$obs$theta.obs)+d$obs$center[2])
		min_y = min(d$obs$r.obs*sin(d$obs$theta.obs)+d$obs$center[2])
		polygon(c(max_x+1, max_x+1, min_x-1, min_x-1), c(max_y+1, min_y-1, min_y-1, max_y+1), fillOddEven = TRUE, col = "white", border = NA)
		lines(x,y, lty = 1, lwd = 1)
            }
            else {
	image(matrix(0,100,100),asp = 1, axes=FALSE,col = 'white')    
		max_x = max(d$obs$r.obs*cos(d$obs$theta.obs)+d$obs$center[1])
		min_x = min(d$obs$r.obs*cos(d$obs$theta.obs)+d$obs$center[1])
		max_y = max(d$obs$r.obs*sin(d$obs$theta.obs)+d$obs$center[2])
		min_y = min(d$obs$r.obs*sin(d$obs$theta.obs)+d$obs$center[2])
		polygon(c(max_x+1, max_x+1, min_x-1, min_x-1), c(max_y+1, min_y-1, min_y-1, max_y+1), fillOddEven = TRUE, col = "white", border = NA)
	            lines(d$x, d$y, lty = 2, lwd = 3)
            }
            polygon(d$ux, d$uy, fillOddEven = TRUE, col = "gray", 
                border = NA)
            polygon(d$lx, d$ly, fillOddEven = TRUE, col = "white", 
                border = NA)
	    lines(d$x, d$y, lty = 2, lwd = 3)
            if (d$shape != "unknown") {
		lines(x, y, lty = 1, lwd = 1)
                lines(d$x, d$y, lty = 2, lwd = 3)
            }
        })
	output$downloadData <- downloadHandler(
    filename = function() { 'image_data.csv' },
    content = function(file) {
      write.csv(data()$subset, file)
    }
  )
    }
    return(shinyApp(ui = ui, server = server))
}



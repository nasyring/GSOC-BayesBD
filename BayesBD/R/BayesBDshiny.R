BayesBDshiny = function () 
{
    ui = pageWithSidebar(
	titlePanel("BayesBD"),
      sidebarPanel(selectInput(inputId = "shape", label = "Choose either an elliptical, triangular, or user-supplied boundary, or indicate that the ground truth is unknown.", 
        choices = c("ellipse", "triangle", "file", "unknown")),
	numericInput(inputId = "mean", label = "Input the default mean boundary radius to begin fitting the model.", 
            value = 0.2, min = 0, max = 1, step = NA, width = NULL),
	conditionalPanel( 
	condition = "input.shape == 'file'",
        fileInput(inputId = "shape_file", label = "Use a custom boundary. The file should be an R script of a function called gamma.fun taking as input an angle in [0, 2pi] and returning the radius of the boundary from a reference point.")), 
        selectInput(inputId = "data_type", label = "Choose to simulate binary or Gaussian data or input image file below.", 
            choices = c("binary sim", "normal sim", "user binary image", 
                "user continuous image")),
	conditionalPanel( 
	condition = "input.data_type == 'user binary image' || input.data_type == 'user continuous image'", 
	textInput(inputId = "data_file",value="", 
            label = "Use image from file. The file should be in .png or .jpeg format. Type the full path here."),
	actionButton(inputId = "go_plot", label = "Display Image")), 
	numericInput(inputId = "centerx", label = "Input the X-coordinate and Y-coordinate of the reference point interior to the boundary function.", 
            value = 0.5, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "centery", label = "Y-coordinate of the reference point.", 
            value = 0.5, min = NA, max = NA, step = NA, width = NULL), 
	conditionalPanel( 
	condition = "input.data_type == 'user binary image' || input.data_type == 'user continuous image'", 
        selectInput(inputId = "pre_fit", label = "Choose if you would like to fit the boundary twice to filter the background.", 
            choices = c("No", "Yes"))), 	
        sliderInput(inputId = "n_burn", label = "Choose a number of posterior samples to burn", 
            value = 1000, min = 500, max = 1000), sliderInput(inputId = "n_run", 
            label = "Choose a number of posterior samples to keep", 
            value = 2000, min = 1000, max = 4000), 
	conditionalPanel( 
	condition = "input.data_type == 'binary sim'",  
        sliderInput(inputId = "p_in", label = "Choose the Bernoulli success probability inside the image", 
            value = 0.5, min = 0, max = 1), sliderInput(inputId = "p_out", 
            label = "Choose the Bernoulli success probability outside the image", 
            value = 0.2, min = 0, max = 1)),
	conditionalPanel( 
	condition = "input.data_type == 'binary sim' || input.data_type == 'user binary image'",
        selectInput(inputId = "ordering", label = "Indicate which region of the image has higher average intensity.", 
            choices = c("Inside", "Outside", "Unknown"))),
	conditionalPanel( 
	condition = "input.data_type == 'normal sim'",
        numericInput(inputId = "mu_in", label = "Mean intensity inside image", 
            value = 1, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "sd_in", label = "Standard deviation inside image", 
            value = 1, min = 0, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "mu_out", label = "Mean intensity outside image", 
            value = 0, min = NA, max = NA, step = NA, width = NULL), 
        numericInput(inputId = "sd_out", label = "Standard deviation outside image", 
            value = 1, min = 0, max = NA, step = NA, width = NULL)), 
	conditionalPanel( 
	condition = "input.data_type == 'normal sim' || input.data_type == 'user continuous image'",
        selectInput(inputId = "ordering_mu", label = "Indicate which region of the image has higher average intensity.", 
            choices = c("Inside", "Outside", "Unknown")),
        selectInput(inputId = "ordering_sd", label = "Indicate which region of the image has higher variation in intensity.", 
            choices = c("Inside", "Outside", "Unknown"))),
        downloadButton('downloadData', 'Download'),
        actionButton(inputId = "go", label = "Update")), 
	mainPanel( 
	 	verbatimTextOutput("info"),
		par(mfrow=c(1,2)),
		plotOutput("image1", click = "plot_click"),
		plotOutput("image")
	)
  		
			
	)
    

     server = function(input, output) {
        theta.plot = seq(from = 0, to = 2 * pi, length.out = 200)
        pre_plot = eventReactive(input$go_plot, {
		image = input$data_file
		cppsamp = fitContImage(image, NULL, c(0,0), .1, 1, 
                  0, 10,'I','I',NULL, FALSE, FALSE)
		return(cppsamp)
	})
	output$image1 <- renderPlot({
    		plotBD(pre_plot(),1)
 	 })
	output$info <- renderText({
    		if(input$data_type == 'user binary image' || input$data_type == 'user continuous image'){paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)}
  	})
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
                gamma.fun = NULL
            }
            if (input$data_type == "binary sim") {
                image = par2obs(m = 100, pi.in = input$p_in, pi.out = input$p_out, 
                  design = "J", center, gamma.fun)
            }
            else if (input$data_type == "normal sim") {
                image = parnormobs(m = 100, mu.in = input$mu_in, 
                  mu.out = input$mu_out, sd.in = input$sd_in, 
                  sd.out = input$sd_out, design = "J", center, 
                  gamma.fun)
            }
            else {
                image = input$data_file
            }
            if (any(input$data_type == "binary sim", input$data_type == 
                "user binary image")) {
		if(input$ordering == 'Inside'){ordero = 'I'
		}else if(input$ordering == 'Outside'){ordero = 'O'
		}else {ordero = 'N'}
                cppsamp1 = fitBinImage(image, gamma.fun, center = center, input$mean, input$n_run, 
                  input$n_burn, 10,ordero,NULL, FALSE, FALSE)
            }
            else {
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
            cppsamp1 = fitContImage(image, gamma.fun, center, input$mean, input$n_run, 
                  input$n_burn, 10,order_mu,order_sd,NULL, FALSE, FALSE)
            }

		theta.plot = seq(from = 0, to = 2*pi, length.out = 200)

		r.est = function(theta){
			thetas = c(theta.plot,2*pi)
			r.thetas = c(cppsamp1$output$upper,cppsamp1$output$upper[1])
			s = sort(c(theta,thetas))
			w = max(which(s==theta))
			lt = s[w-1]
			ut = s[w+1]
			lr = r.thetas[w-1]
			ur = r.thetas[w]
			r_est = ((theta - lt)/(ut-lt))*ur+((ut - theta)/(ut-lt))*lr
			return(r_est[1])
		}

		app.r.est = function(theta) apply(matrix(theta,length(theta),1),1,r.est)

		r_ests = app.r.est(cppsamp1$obs$theta.obs)
		r_ests = matrix(r_ests,length(cppsamp1$obs$theta.obs),1)
		subset = ifelse(cppsamp1$obs$r.obs<=r_ests,1,0)

	    if(input$pre_fit == 'Yes'){



			if (any(input$data_type == "binary sim", input$data_type == "user binary image")) {
				cppsamp2 = fitBinImage(image, gamma.fun, center = center, input$mean, input$n_run, 
                  				input$n_burn, 10,'N',subset, FALSE, FALSE)							
			}else {
				cppsamp2 = fitContImage(image, gamma.fun, center, input$mean, input$n_run, 
                 			 input$n_burn, 10,'N','N',subset, FALSE, FALSE)
			}
		r_ests2 = app.r.est(cppsamp2$obs$theta.obs)
		r_ests2 = matrix(r_ests2,length(cppsamp2$obs$theta.obs),1)
		subset2 = ifelse(cppsamp2$obs$r.obs<=r_ests2,1,0)
		return(list(cppsamp1=cppsamp1,cppsamp2=cppsamp2,subset=subset,subset2=subset2))
	    }else {
		return(list(cppsamp1=cppsamp1,subset=subset))	
		}


        })
        output$image = renderPlot({
            d = data()
           if(input$pre_fit=='Yes'){

			par(mfrow = c(1,2))
			plotBD(d$cppsamp1,3)
			plotBD(d$cppsamp2,3)

		}else {

			plotBD(d$cppsamp1,3)
		
		}
        })
output$downloadData <- downloadHandler(
    filename = function() { 'image_data.txt' },
    content = function(file) {
	d=data()
    if(input$pre_fit == 'Yes'){
	out = cbind(as.vector(d$cppsamp1$obs$r.obs), as.vector(d$cppsamp1$obs$theta.obs), as.vector(d$cppsamp1$obs$intensity), as.vector(d$cppsamp1$subset),as.vector( d$cppsamp1$subset2))
    }else {
	out = cbind(as.vector(d$cppsamp1$obs$r.obs), as.vector(d$cppsamp1$obs$theta.obs), as.vector(d$cppsamp1$obs$intensity), as.vector(d$cppsamp1$subset))
    }
      write.table(out, file)
    }
  )
    }
    return(shinyApp(ui = ui, server = server))
}



BayesBDshiny <-
function(){
require(plotrix)
my.radial = function(r, theta, ...){
radial.plot(c(r[order(theta)]), c(theta[order(theta)]),
rp.type = "p", show.grid.label = TRUE, radial.lim = c(0, 0.5),
...)
}
ui = fluidPage(
	selectInput(inputId = "shape", label = "Choose either an elliptical or triangular boundary or indicate that you have uploaded a file below.", choices = c("ellipse", "triangle", "file")),
	fileInput(inputId = "shape_file", label = "Use a custom boundary. The file should be an R script of a function called taking as input an angle in [0, 2pi] and returning the radius of the boundary from a reference point."),
	selectInput(inputId = "data_type", label = "Choose to simulate binary or Gaussian data or input data file below.", choices = c("binary sim", "normal sim", "user binary data", "user normal data")),
	fileInput(inputId = "data_file", label = "Use image data from file. The file should be in table format with no headers. The first column is for intensity, the second gives the radius from the origin, and the third gives the angle from the origin and positive x-axis."),
	sliderInput(inputId = "n_burn",
	  label = "Choose a number of posterior samples to burn",
	  value = 1000, min = 500, max = 1000),
	sliderInput(inputId = "n_run",
	  label = "Choose a number of posterior samples to keep",
	  value = 1000, min = 500, max = 1000),
	"Use the following inputs for Binary simulations",
	sliderInput(inputId = "p_in",
	  label = "Choose the Bernoulli success probability inside the image",
	  value = 0.5, min = 0, max = 1),
	sliderInput(inputId = "p_out",
	  label = "Choose the Bernoulli success probability outside the image",
	  value = 0.2, min = 0, max = 1),
	"Use the following inputs for Gaussian simulations",
	numericInput(inputId = "mu_in", label = "Mean intensity inside image", value = 1, min = NA, max = NA, step = NA,
  	width = NULL),
	numericInput(inputId = "sd_in", label = "Standard deviation inside image", value = 1, min = 0, max = NA, step = NA,
  	width = NULL),
	numericInput(inputId = "mu_out", label = "Mean intensity outside image", value = 0, min = NA, max = NA, step = NA,
  	width = NULL),
	numericInput(inputId = "sd_out", label = "Standard deviation outside image", value = 1, min = 0, max = NA, step = NA,
  	width = NULL),
	actionButton(inputId = "go", label = "Update"),
	plotOutput("image")
)

server = function(input, output){
	theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
	rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees
	data = eventReactive(input$go, {	

	 if(input$shape=="ellipse"){
		gamma.fun = ellipse(a = 0.35, b = 0.25)
	 }else if(input$shape=="triangle"){
		gamma.fun = triangle2(0.5)
	 }else{
		gamma.fun = source(input$shape_file$datapath)$value 
	 }	
	 if(input$data_type=='binary sim'){
	 	obs = par2obs(m = 100, pi.in = input$p_in, pi.out = input$p_out, design = 'J', gamma.fun)
	 }else if(input$data_type == 'normal sim'){
		obs = parnormobs(m = 100, mu.in = input$mu_in, mu.out = input$mu_out,sd.in = input$sd_in, sd.out = input$sd_out, design = 'J', gamma.fun)
       }else{
		obs.list = read.table(input$data_file$datapath) 
		obs = list()
		obs$intensity = matrix(obs.list[,1],sqrt(length(obs.list[,1])),sqrt(length(obs.list[,1])))
	 	obs$r.obs = matrix(obs.list[,2],sqrt(length(obs.list[,2])),sqrt(length(obs.list[,2])))
	 	obs$theta.obs = matrix(obs.list[,3],sqrt(length(obs.list[,3])),sqrt(length(obs.list[,3])))
	 } 
	  if(any(input$data_type == 'binary sim',input$data_type == 'user binary data')){
	  	cppsamp = BayesBDbinary(obs, .4, input$n_run, input$n_burn, 10, FALSE)
	  }else {
	  	cppsamp = BayesBDnormal(obs, .4, input$n_run, input$n_burn, 10, FALSE)
	  }
	  x = cppsamp$estimate*cos(cppsamp$theta)+obs$center[1]
	  y = cppsamp$estimate*sin(cppsamp$theta) +obs$center[2]
	  lx = cppsamp$lower*cos(cppsamp$theta)
        ly = cppsamp$lower*sin(cppsamp$theta)
	  ux = cppsamp$estimate*cos(cppsamp$theta)
	  uy = cppsamp$estimate*sin(cppsamp$theta)
	  return(list(obs = obs, x=x, y=y, lx = lx, ly = ly, ux = ux, uy = uy, gamma.fun = gamma.fun))
	})
	

	output$image = renderPlot({
	  d = data()
	  par(mfrow = c(1, 2))
	  image(rotate(d$obs$intensity), axes = FALSE, asp = 1, main = 'observation')
	  lines(d$x,d$y, lty = 2, lwd = 3)
	  my.radial(d$gamma.fun(theta.plot), theta.plot, line.col = 1, lty = 1, lwd = 1, show.grid = FALSE)
	  polygon(d$ux, d$uy,fillOddEven = TRUE, col = 'gray', border = NA)
	  polygon(d$lx, d$ly,fillOddEven = TRUE, col = 'white', border = NA)
	  lines(d$x-d$obs$center[1],d$y-d$obs$center[2], lty = 2, lwd = 3)
      })
}


return(shinyApp(ui = ui, server = server))}

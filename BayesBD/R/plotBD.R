plotBD =
function (obs, samples = list()) 
{
if(is.null(samples$estimate)){
xobs = obs$r.obs*cos(obs$theta.obs)+obs$center[1]
    yobs = obs$r.obs*sin(obs$theta.obs)+obs$center[2]
    if(min(obs$intensity)<0){ 
       normalized = (obs$intensity+abs(min(obs$intensity)))/(max(obs$intensity)-min(obs$intensity))
    }else{
       normalized = (obs$intensity-abs(min(obs$intensity)))/(max(obs$intensity)-min(obs$intensity))
    }
    plot(xobs, yobs, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
}else{
    par(mfrow = c(1, 2))
    xobs = obs$r.obs*cos(obs$theta.obs)+obs$center[1]
    yobs = obs$r.obs*sin(obs$theta.obs)+obs$center[2]
    if(min(obs$intensity)<0){ 
       normalized = (obs$intensity+abs(min(obs$intensity)))/(max(obs$intensity)-min(obs$intensity))
    }else{
       normalized = (obs$intensity-abs(min(obs$intensity)))/(max(obs$intensity)-min(obs$intensity))
    }
    plot(xobs, yobs, col = gray(normalized), pch = 15, cex = 0.375, axes = FALSE, xlab = '', ylab = '',asp = 1)
    x = samples$estimate * cos(samples$theta) + obs$center[1]
    y = samples$estimate * sin(samples$theta) + obs$center[2]
    theta.plot = samples$theta
    image(matrix(0,10,10), col = "white", axes = FALSE, asp = 1)
    lines(x, y, lty = 2, lwd = 3)
    x = samples$upper * cos(samples$theta) + obs$center[1]
    y = samples$upper * sin(samples$theta) + obs$center[2]
    polygon(x, y, fillOddEven = TRUE, col = "gray", border = NA)
    x = samples$lower * cos(samples$theta) + obs$center[1]
    y = samples$lower * sin(samples$theta) + obs$center[2]
    polygon(x, y, fillOddEven = TRUE, col = "white", border = NA)
    x = samples$estimate * cos(samples$theta) + obs$center[1]
    y = samples$estimate * sin(samples$theta) + obs$center[2]
    lines(x, y, lty = 2, lwd = 3)
    if (!is.null(obs$gamma.fun)) {
        x = obs$gamma.fun(theta.plot) * cos(theta.plot) + obs$center[1]
        y = obs$gamma.fun(theta.plot) * sin(theta.plot) + obs$center[2]
    }
    lines(x, y, lty = 1, lwd = 1)
}
}


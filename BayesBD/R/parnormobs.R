parnormobs <-
function (m, mu.in, mu.out, sd.in, sd.out, design, center, gamma.fun) {
    obs <- matrix(NA, m, m)
    if (design == "D") {
        x.axis = (col(obs) - 1)/m + 1/(2 * m)
        y.axis = (m - row(obs))/m + 1/(2 * m)
    }
    if (design == "J") {
        x.axis = (col(obs) - 1)/m + 1/(2 * m) + runif(m^2, min = -1/(2 * 
            m), max = 1/(2 * m))
        y.axis = (m - row(obs))/m + 1/(2 * m) + runif(m^2, min = -1/(2 * 
            m), max = 1/(2 * m))
    }
    if (design == "U") {
	x.axis = matrix(runif(m^2, 0, 1),m,m)
	y.axis = matrix(runif(m^2, 0, 1),m,m)
    }
    r.obs = sqrt((x.axis - center[1])^2 + (y.axis - center[2])^2)
    theta.obs <- atan2(y.axis - center[2], x.axis - center[1])
    theta.obs[theta.obs < 0] = theta.obs[theta.obs < 0] + 2 * 
        pi
    obsLabel = (r.obs < gamma.fun(theta.obs))
    n.In = sum(obsLabel)
    n.Out = sum(!obsLabel)
    obs[obsLabel] = rnorm(n.In, mean = mu.in, sd = sd.in)
    obs[!obsLabel] = rnorm(n.Out, mean = mu.out, sd = sd.out)
    return(list(intensity = obs, theta.obs = theta.obs, r.obs = r.obs, 
        center = center, x = x.axis, y = y.axis, gamma.fun = gamma.fun))
}

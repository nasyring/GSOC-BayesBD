par2obs <-
function (m, pi.in, pi.out, design, gamma.fun) 
{
    center = c(0.5, 0.5)
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
    r.obs = sqrt((x.axis - center[1])^2 + (y.axis - center[2])^2)
    theta.obs <- atan2(y.axis - 1/2, x.axis - 1/2)
    theta.obs[theta.obs < 0] = theta.obs[theta.obs < 0] + 2 * 
        pi
    obsLabel = (r.obs < gamma.fun(theta.obs))
    n.In = sum(obsLabel)
    n.Out = sum(!obsLabel)
    obs[obsLabel] = rbinom(n.In, size = 1, prob = pi.in)
    obs[!obsLabel] = rbinom(n.Out, size = 1, prob = pi.out)
    return(list(intensity = obs, theta.obs = theta.obs, r.obs = r.obs, 
        center = center, x = x.axis, y = y.axis))
}

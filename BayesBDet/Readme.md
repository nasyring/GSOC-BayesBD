Bayesian bounday detection 
================================

This package detects boundaries in images. 

### Version histories: 
v0.2 
- c++ implementation with ** times speed gain
- add Gaussian noised images 
- improve visualization 
- ... 
     
v0.1 
- implemented boundary detection methods for binary images 

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('BayesBDet', 'nasyring')
```

### Example

```S
set.seed(2015)
% ellipse boundary
gamma.fun = ellipse(a = 0.35, b = 0.25)
obs = par2obs(m = 100, pi.in = 0.5, pi.out = 0.2, design = 'J', gamma.fun)

# Bayes Estimation:
# it takes around 7min if runs 10000 iterations: saved in 'data.Rdata'
# BayesEst = BayesBD.binary(obs, n.run = 10000, n.burn = 1000)
data(data)

# visualize the estimates
theta.plot = seq(from = 0, to = 2*pi, length.out = 200)
gamma.hat.theta = BayesEst$gamma.hat(theta.plot)

## plotting utilities
require(plotrix)
my.radial <- function(r, theta, ...){
  radial.plot(c(r[order(theta)]), c(theta[order(theta)]),
              rp.type = "p", show.grid.label = TRUE, radial.lim = c(0, 0.5),
              ...)
}
# rotate a matrix
rotate <- function(x) t(apply(x, 2, rev))  # rotate closewise by 90 degrees

par(mfrow = c(1, 2))
# rotate & image it - square (asp = 1)
image(rotate(obs$intensity), axes = FALSE, asp = 1, main = 'observation')
my.radial(gamma.fun(theta.plot), theta.plot, line.col = 1, lty = 2, lwd = 2,
          main = 'Estimated boundary vs. True', show.grid = FALSE)
my.radial(gamma.hat.theta, theta.plot, add = TRUE,
         line.col = 'red', lty = 2, lwd = 2, show.grid = FALSE)
```

### Reference
Li, M. and Ghosal, S. (2016). Bayesian Detection of Image Boundaries. _arXiv preprint_.  (http://arxiv.org/abs/1508.05847)


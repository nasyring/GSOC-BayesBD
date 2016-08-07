triangle2 <-
function(S){
triangle.scalar <- function(theta, h) {
r=0
        if (any(theta >= 0 & theta <  pi/2, theta>= 11*pi/6 & theta < 2*pi)) {
      h = S*sin(pi/3)
            k = S*cos(pi/3)
hstar = k*tan(pi/6)
hh = h - hstar
r = hh / ((sin(theta) + (h/k)*cos(theta)))
        }
        if (theta>= pi/2 & theta < 7*pi/6) {
            theta = pi - theta
      h = S*sin(pi/3)
k = S*cos(pi/3)
hstar = k*tan(pi/6)
hh = h - hstar
r = hh / ((sin(theta) + (h/k)*cos(theta)))
        }
        if (theta>= 7*pi/6 & theta < 3*pi/2) {
           a = 3*pi/2-theta
           h = S*sin(pi/3)
k = S*cos(pi/3)
hstar = k*tan(pi/6)
     r = (hstar)/cos(a) 
        }
        if (theta>=3*pi/2 & theta< 11*pi/6) {
           a = theta-3*pi/2
           h = S*sin(pi/3)
k = S*cos(pi/3)
hstar = k*tan(pi/6)
     r = (hstar)/cos(a) 
        }
        return(r)
    }
    ret = function(theta) c(sapply(theta, function(theta) triangle.scalar(theta, 
        h)))
    return(ret)
}

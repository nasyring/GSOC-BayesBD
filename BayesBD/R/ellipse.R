ellipse <-
function (a, b, r0 = 0, theta0 = 0, phi = 0) 
{
    function(theta) {
        P = r0 * ((b^2 - a^2) * cos(theta + theta0 - 2 * phi) + 
            (a^2 + b^2) * cos(theta - theta0))
        R = (b^2 - a^2) * cos(2 * theta - 2 * phi) + a^2 + b^2
        Q = sqrt(2) * a * b * sqrt(R - 2 * r0^2 * sin(theta - 
            theta0)^2)
        r = (P + Q)/R
        return(r)
    }
}

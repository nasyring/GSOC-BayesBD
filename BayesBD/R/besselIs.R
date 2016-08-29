besselIs <-
function (x, nu, expon.scaled = FALSE) 
{
    besselI(x, nu, 1 + as.logical(expon.scaled))
}

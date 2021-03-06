\name{loglik_cDCC}
\alias{loglik_cDCC}
\title{
Loss function used in the cDCC robust estimation.
}
\description{
Loss function used in the cDCC robust estimation.
}
\usage{
loglik_cDCC(par, Qb, s, sigma)
}
\arguments{
  \item{par}{
Two-dimensional cDCC vector parameters 
}
  \item{Qb}{
Qbar matrix (obtained from the Robust_cDCC function)
}
  \item{s}{
Devolatilized returns
}
  \item{sigma}{
Sigma parameter (which is computed inside the Robust_cDCC function). In a three-dimensional case this value is equal to 0.8309765.
}
}
\details{
This function is used in the robust estimation. We can use it to evaluate the value of the robust cDCC loss function using several values of the vector parameters.
}
\value{
Returns the value of the loss function.
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and value-at-risk. Journal of Statistical Computation and Simulation 88.10 (2018): 1976-2000.
}
\author{
Carlos Trucíos
}
\examples{
\donttest{
# Estimating the parameters of the cDCC model in a robust way.
cDCC = Robust_cDCC(toyexampledata[,1:3])
param = cDCC[[1]]
Qbar = cDCC[[2]]
vol1 = fitted_Vol(param[1:3],toyexampledata[,1])
vol2 = fitted_Vol(param[4:6],toyexampledata[,2])
vol3 = fitted_Vol(param[7:9],toyexampledata[,3])
e = matrix(c(toyexampledata[,1]/vol1[1:nrow(toyexampledata)],
             toyexampledata[,2]/vol2[1:nrow(toyexampledata)],
             toyexampledata[,3]/vol3[1:nrow(toyexampledata)]), ncol=3)
      
loglik_cDCC(param[10:11],Qbar,e, 0.8309765)      
}
}
\keyword{DCC}
\keyword{Jump}
\keyword{Loss function}
\keyword{Outliers}
\keyword{Robustness}
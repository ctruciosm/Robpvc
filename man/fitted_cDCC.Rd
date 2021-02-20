\name{fitted_cDCC}
\alias{fitted_cDCC}
\title{
Conditional covariance matrix using the robust cDCC Estimator
}
\description{
Conditional covariance matrix using the robust cDCC estimator of Boud et al. (2013) with the modification introduced by Trucíos et at. (2018).
}
\usage{
fitted_cDCC(r, Qbar, params)
}
\arguments{
  \item{r}{
Matrix of time series returns.
}
  \item{Qbar}{
Qbar matrix obtained from Robust_cDCC function
}
  \item{params}{
Estimated parameters obtained from Robust_cDCC function
}
}
\details{
More details can be found in Boudt et al. (2013) and Trucíos et at. (2018).
}
\value{
The function returns the estimated conditional covariance matrix from t = 1, ..., T+1, where T is the length of the sample size.
}
\references{
  Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and Value-at-Risk. Journal of Statistical Computation and Simulation 88.10 (2018): 1976-2000.
}
\author{
Carlos Trucíos
}
\examples{
\donttest{
# Estimating the parameters of the cDCC model in a robust way.
cDCC = Robust_cDCC(toyexampledata[,1:3])
# Estimated conditional covariance matrix
H = fitted_cDCC(toyexampledata[,1:3], cDCC[[2]], cDCC[[1]])
# One-step-ahead conditional covariance matrix
H[[nrow(toyexampledata[,1:3]) + 1]]
}
}

\keyword{DCC}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness}
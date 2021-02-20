\name{rpvc_VaR}
\alias{rpvc_VaR}
\title{Value-at-Risk via Robust Principal Volatility Components}
\usage{
rpvc_VaR(rtn, m = 10, c = 0.99, k = 1, weights = NULL, alpha = c(0.01, 0.05))
}
\arguments{
\item{rtn}{Matrix of time series data}

\item{m}{lag parameter, by default m = 10}

\item{c}{threshold value, by default c = 0.99 (observation with largest 1 percent Mahalanobis distance are considered as extreme observations).}

\item{k}{number of volatility components}

\item{weights}{vector of portfolio weights, if NULL the equal-weight portffolio is used}

\item{alpha}{risk leveñ, iby default alpha = c(0.01, 0.05)}
}
\value{
One-step-ahead Value-at-Risk alpha risk level
}
\description{
Forecast the one-step-ahead Value-at-Risk (alpha risk level) using the robust principal volatility component procedure as in Trucíos et al. (2019).
}
\examples{
\donttest{
rpvc_VaR(toyexampledata, m = 10, c = 0.99, k = 1, alpha = c(0.01, 0.025, 0.05))
}
}

\references{
  Trucíos, Carlos, Luiz K. Hotta, and Pedro Valls-Pereira. On the robustness of the principal volatility components. Journal of Empirical Finance 52.1 (2019): 1201-2019.
}


\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness} 
\keyword{Value-at-Risk}
\keyword{Volatility Components}

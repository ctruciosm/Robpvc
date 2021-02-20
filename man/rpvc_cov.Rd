\name{rpvc_cov}
\alias{rpvc_cov}
\title{Conditional Covariance Matrix using Robust Principal Volatility Components}
\usage{
rpvc_cov(rtn, m = 10, c = 0.99, k = 1)
}
\arguments{
\item{rtn}{Matrix of time series data}

\item{m}{lag parameter, by default m = 10}

\item{c}{threshold value, by default c = 0.99 (observation with largest 1 percent Mahalanobis distance are considered as extreme observations).}

\item{k}{number of volatility components}
}
\value{
One-step-ahead conditional covariance matrix
}
\description{
Estimates the one-step-ahead conditional covariance matrix using the robust principal volatlity component approach of Trucíos et al. (2019) along with the robust GARCH/cDCC procedure of Boudt et al. (2013) with the modifications introduced by Trucíos et al. (2017, 2018).
}
\examples{
\donttest{
rpvc_cov(toyexampledata, m = 10, c = 0.99, k = 1)
}
}

\references{
 Boudt, Kris, Jon Danielsson, and Sébastien Laurent. Robust forecasting of dynamic conditional correlation GARCH models. International Journal of Forecasting 29.2 (2013): 244-257.
 
 Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation 87.16 (2017): 3152-3174.
  
  Trucíos, Carlos, Luiz K. Hotta, and Esther Ruiz. Robust bootstrap densities for dynamic conditional correlations: implications for portfolio selection and Value-at-Risk. Journal of Statistical Computation and Simulation 88.10 (2018): 1976-2000.
  
  Trucíos, Carlos, Luiz K. Hotta, and Pedro Valls-Pereira. On the robustness of the principal volatility components. Journal of Empirical Finance 52.1 (2019): 1201-2019.
}

\keyword{Conditional Covariance Matrix}
\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness} 
\keyword{Volatility Components}

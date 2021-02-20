\name{Robpvc-package}
\alias{Robpvc-package}
\docType{package}
\title{
  Robust Princial Volatility Components
}
\description{
Principal volatility components and one-step-ahead conditional covariance matrix in a robust-to-outliers way as in Trucíos et al. (2019). The package also implements the robust cDCC (Dynamic Conditional Correlation) estimator of Boudt et al. (2013) with the modification introduced by Trucíos et al. (2018).
}
\details{
This package provides a robust procedure to obtain principal volatility components. Aditionally, the conditional covariance matrix obtained via robust principal volatility components along with the robust GARCH/cDCC procedure of Boudt et al. (2013) with the modification introduced by Trucíos et at. (2017, 2018) is also implemented.
This procedure has shown good finite sample properties in both Monte Carlo experiments and empirical data. See;  Trucíos et at. (2019) for more information.
}
\author{
Carlos Trucíos <carlos.trucios@facc.ufrj.br>
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

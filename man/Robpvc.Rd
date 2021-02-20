\name{Robpvc}
\alias{Robpvc}
\title{Robust Principal Volatility Components}
\usage{
Robpvc(rtn, m = 10, c = 0.99)
}
\arguments{
\item{rtn}{Matrix of time series data}

\item{m}{lag parameter, by default m = 10}

\item{c}{threshold value, by default c = 0.99 (observation with largest 1 percent Mahalanobis distance are considered as extreme observations).}
}
\value{
None
}
\description{
Estimates the principal volatility components in a robust way using the procedure of Trucíos et al. (2019)
}
\examples{
\donttest{
Robpvc(toyexampledata)
}
}
\references{
  Trucíos, Carlos, Luiz K. Hotta, and Pedro Valls-Pereira. On the robustness of the principal volatility components. Journal of Empirical Finance 52.1 (2019): 1201-2019.
}

\keyword{Jump}
\keyword{Outliers}
\keyword{Robustness} 
\keyword{Volatility Components} 

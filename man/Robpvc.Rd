% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpvc.R
\name{Robpvc}
\alias{Robpvc}
\title{Robust Principal Volatility Components}
\usage{
Robpvc(rtn, m = 10, c = 0.99)
}
\arguments{
\item{rtn}{Matrix of time series data}

\item{m}{lag parameter, by default m = 10}

\item{c}{threshold value, by default c = 0.99 which means that the observation with largest 1 percent Mahalanobis distance are considered as extreme observations.}
}
\value{
None
}
\description{
Estimates the principal volatility components in a robust way
}
\examples{
\donttest{
Robpvc(toyexampledata)
}

}
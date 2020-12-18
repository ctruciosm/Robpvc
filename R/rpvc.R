#' Robust Principal Volatility Components
#'
#' Estimates the principal volatility components in a robust way
#'
#' @param rtn Matrix of time series data
#' @param m lag parameter, by default m = 10
#' @param c threshold value, by default c = 0.99 which means that the observation with largest 1 percent Mahalanobis distance are considered as extreme observations.
#'
#' @return None
#'
#' @examples
#' \donttest{
#' Robpvc(toyexampledata)
#' }
#'
#' @export
#' @import Rcpp RcppArmadillo
#' @importFrom DetMCD DetMCD
#' @importFrom stats quantile
Robpvc = function (rtn, m = 10, c = 0.99) {
  if(!is.matrix(rtn)) rtn = as.matrix(rtn)
  nT = nrow(rtn)
  k = ncol(rtn)
  Cov_AUX = DetMCD(rtn)
  pre_x = rtn-matrix(rep(Cov_AUX$center,nT),ncol=k,byrow=T)
  RM = list()
  d2 = d2_aux = rep(0,nT)
  RM[[1]] = Cov_AUX$cov
  d2_aux = rowSums((pre_x %*% solve(RM[[1]])) * pre_x)
  Lim_aux = quantile(d2_aux,c)
  d2 = robRmRcpp(RM[[1]], pre_x, Lim_aux, d2_aux)
  Lim = quantile(d2,c)
  INDEX = ifelse(d2>Lim,1,0)
  Weight1 = ifelse(d2<=Lim,1,1/d2)
  Weight = c(rep(0,k),Weight1[c((k+1):nT)])
  S_W = sum(Weight)
  S_W1 = sum(Weight1)
  A = diag(tcrossprod(pre_x))
  ind = A
  Sum1 = matrix(0, k, k)
  for (kk in 1:m) {
    Sum2 = matrix(0, k, k)
    for (tt in 1:nT){
      sign = (A<= ind[tt])
      signnew = Weight1*c(rep(0,kk),sign[1:(nT-kk)])
      temp1 = rtn*matrix(sqrt(signnew),ncol=k,nrow=nT)
      temp2 = (crossprod(temp1)-sum(signnew)*RM[[1]])/S_W
      Sum2 = Sum2 +  temp2%*%temp2/nT
      #Sum2 = Sum2 +  Weight1[tt]*temp2%*%temp2/S_W1
    }
    Sum1 = Sum1 + Sum2
  }
  m2 = eigen(Sum1)
  Valu = m2$values
  Prop = Valu/sum(Valu)
  Vec = m2$vectors
  return(list(residuals = pre_x, values = m2$values, vectors = m2$vectors,COV = Cov_AUX$cov,index = INDEX, D2 = d2))
}


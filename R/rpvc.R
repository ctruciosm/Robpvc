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
  S_W1 = sum(Weight1)
  A = diag(tcrossprod(pre_x))
  ind = A
  Sum1 = matrix(0, k, k)
  for (kk in 1:m) {
    Weight = c(rep(0,kk),Weight1[c((kk+1):nT)])
    S_W = sum(Weight)
    Sum2 = matrix(0, k, k)
    for (tt in 1:nT){
      sign = (A<= ind[tt])
      signnew = Weight*c(rep(0,kk),sign[1:(nT-kk)]) 
      temp1 = rtn*matrix(sqrt(signnew),ncol=k,nrow=nT)
      temp2 = (crossprod(temp1)-sum(signnew)*RM[[1]])/S_W
      #Sum2 = Sum2 +  temp2%*%temp2/nT
      Sum2 = Sum2 +  Weight1[tt]*temp2%*%temp2/S_W1
    }
    Sum1 = Sum1 + Sum2
  }
  m2 = eigen(Sum1)
  return(list(residuals = pre_x, values = m2$values, vectors = m2$vectors,COV = Cov_AUX$cov,index = INDEX, D2 = d2))
}


#' @noRd
#' @importFrom stats cor qchisq pchisq
Q_bar = function(s, k = 30, gamma = 0.95){
  n = nrow(s)
  p = ncol(s)
  AUXL = 20
  Lt = rep(0,n)
  Su = matrix(0,p,p)
  for(i in 1:n){
    if (i<=k/2){
      Ct = cor(s[1:(k+1),],method="spearman")
    } else {
      if (i>n-k/2){
        Ct = cor(s[(n-k):n,],method="spearman")
      } else {
        Ct = cor(s[(i-k/2):(i+k/2),],method="spearman")
      }
    }
    SCt = 2*sin(1/6*pi*Ct)
    AUXL = try(t(s[i,])%*%solve(SCt)%*%s[i,],silent = TRUE)
    Lt[i] = ifelse( AUXL <= qchisq(gamma,p),1,0)       
    Su = Su + s[i,]%*%t(s[i,])*Lt[i]
  }
  CN = (p/(p*pchisq(qchisq(gamma,p),p+2) + (1-gamma)*qchisq(gamma,p)))
  RC = CN*Su/sum(Lt)  
  D = solve(sqrt(diag(diag(RC))))
  R = D%*%RC%*%D
  return(R)
}


#' @export
#' @import Rcpp RcppArmadillo
#' @importFrom stats constrOptim pchisq qchisq integrate dchisq
#' @importFrom RobGARCHBoot ROBUSTGARCH fitted_Vol
Robust_cDCC = function(r){
  n = nrow(r)
  p = ncol(r)
  coef = matrix(0, ncol = p, nrow = 3)
  vol = matrix(0,ncol = p, nrow = n+1)
  e = matrix(0, ncol = p, nrow = n)
  for (i in 1:p){
    coef[,i] = ROBUSTGARCH(r[,i])
    vol[,i] = fitted_Vol(coef[,i], r[,i])
    e[,i] = r[,i]/vol[1:n,i]
  }
  
  if (p == 2) sigma_ = 0.8257925
  if (p == 3) sigma_ = 0.8309765
  if (p == 4) sigma_ = 0.8384375
  if (p == 5) sigma_ = 0.8466635
  if (p > 5){
    integrand = function(x) { (p+4)*x/(2+x)*dchisq(x,p) }
    sigma_ = p/integrate(integrand,0,Inf)$value
  }
  Qbarra = Q_bar(e)
  parini = gridcDCC(Qbarra,e, sigma_)
  
  ra = matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb = c(0.00001, 0.00001,-0.9999)
  coef_cDCC = constrOptim(theta = parini, f = loglik_cDCC, grad = NULL, ui = ra, ci = rb, Qb = Qbarra,s = e, sigma = sigma_)$par
  
  coeff = c(as.vector(coef),coef_cDCC)
  return(list(coeff,Qbarra))
}


#' @export
#' @noRd
#' @import Rcpp RcppArmadillo
#' @importFrom RobGARCHBoot fitted_Vol
fitted_cDCC = function(r, Qbar, params){
  Dim = dim(r)
  n = Dim[1]
  p = Dim[2]
  H = list()
  vol = matrix(0,ncol = p, nrow = n+1)
  e = matrix(0, ncol = p, nrow = n)
  for (i in 1:p){
    coef = params[(3*(i-1)+1):(3*i)]
    vol[,i] = fitted_Vol(coef, r[,i])
    e[,i] = r[,i]/vol[1:n,i]
  }
  
  dccpar = params[(3*p+1):(3*p+2)]
  R = cor_cDCC(dccpar,Qbar,e)
  
  for(j in 1:(n+1)){
    H[[j]] = diag(vol[j,])%*%R[,,j]%*%diag(vol[j,])
  }
  return(H)
}


#' @export
#' @import Rcpp RcppArmadillo
#' @importFrom DetMCD DetMCD
#' @importFrom stats quantile
#' @importFrom utils tail
rpvc_cov = function(rtn, m = 10, c = 0.99, k = 1){
  Y = scale(as.matrix(rtn), scale=FALSE)
  n = nrow(Y)
  m1y = Robpvc(Y, m, c)
  M1 = m1y$vectors[,1:k]
  M2 = m1y$vectors[,-c(1:k)]
  X = Y%*%M1 
  if (k == 1){
    coeff = ROBUSTGARCH(X)
    sigma2 = tail(fitted_Vol(coeff, X),1)^2
    Hhat = sigma2*M1%*%t(M1) + m1y$COV%*%M2%*%t(M2) + M2%*%t(M2)%*%m1y$COV%*%M1%*%t(M1)  
  } else{
    coeff_AUX = Robust_cDCC(X)
    coeff = coeff_AUX[[1]]
    Qbarra = coeff_AUX[[2]]
    H = fitted_cDCC(X, Qbar = Qbarra, params = coeff)
    Hhat = M1%*%as.matrix(H[[n+1]])%*%t(M1) + m1y$COV%*%M2%*%t(M2) + M2%*%t(M2)%*%m1y$COV%*%M1%*%t(M1)  
  }
  return(Hhat)
}



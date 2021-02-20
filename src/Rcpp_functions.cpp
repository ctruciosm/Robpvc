#include<RcppArmadillo.h>
#include<Rcpp.h>
#include<Rmath.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
//' @noRd
//' @import Rcpp RcppArmadillo
//' @useDynLib Robpvc
// [[Rcpp::export]]
SEXP robRmRcpp(arma::mat sigma1, arma::mat y, double Lim, NumericVector aux){
  int N = y.n_rows, K = sigma1.n_cols;
  arma::mat sigmaux(K,K);
  NumericVector d2(N);

  d2[0] = arma::mat(y.row(0)*sigma1.i()*y.row(0).t())(0,0);
  sigmaux = sigma1;

  for(int i=1; i<N; i++){
    if(aux[i-1]<Lim){
      sigmaux = 0.01*y.row(i-1).t()*y.row(i-1) + 0.99*sigmaux;
    }else {
      sigmaux = 0.01*sigma1 + 0.99*sigmaux;
    }
    d2[i] = arma::mat(y.row(i)*sigmaux.i()*y.row(i).t())(0,0);
    }

  return Rcpp::wrap(d2);
}



// [[Rcpp::depends(RcppArmadillo)]]
//' @noRd
//' @import Rcpp RcppArmadillo
//' @useDynLib Robpvc
// [[Rcpp::export]]
SEXP gridcDCC(arma::mat Qb,arma::mat s, double sigma){
  NumericVector coeff(2),vi(2);
  double alfa1, beta1;
  double alfa1min = 0.01,  alfa1max = 0.3, beta1min = 0.65,  beta1max = 0.98,  nalfa1 = 5,  nbeta1 = 5;
  double ml = 100000000, nml;
  double lmalfa1 = (alfa1max-alfa1min)/nalfa1;
  double lmbeta1 = (beta1max-beta1min)/nbeta1;
  Rcpp::Function loglik_cDCC("loglik_cDCC");
  
  for(int nj=0; nj<nalfa1; nj++){
    for(int nk=0; nk<nbeta1; nk++){
      alfa1 = alfa1min+nj*lmalfa1;
      beta1 = beta1min+nk*lmbeta1; 
      
      if(alfa1+beta1<0.999){
        coeff[0] = alfa1;
        coeff[1] = beta1;
        nml=Rcpp::as<double>(loglik_cDCC(coeff,Qb,s, sigma));
        if (nml<ml){
          vi[0] = coeff[0];
          vi[1] = coeff[1];
          ml=nml;
        }
      }
    }
  }
  return(vi);
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @export
//' @noRd
//' @useDynLib Robpvc
//' @import Rcpp RcppArmadillo
// [[Rcpp::export]]
SEXP loglik_cDCC(arma::vec par,arma::mat Qb,arma::mat s, double sigma){
  int T = s.n_rows;
  int K = s.n_cols;
  Rcpp::Function pchisq("pchisq");
  Rcpp::Function qchisq("qchisq");
  double CN = K/(K*as<double>(pchisq(qchisq(0.9973,K),K+2)) + as<double>(qchisq(0.9973,K))*(1-0.9973));
  arma::vec lR1(T), lR2(T), lR3(T), d(T);
  arma::mat  R(K,K), Qt(K,K),Pt(K,K), iPt(K,K);
  R.zeros();  Qt.zeros();  Pt.zeros();
  R = Qb;
  Qt = Qb;
  lR1[0] = log(det(R));
  d[0] = arma::conv_to<double>::from(s.row(0)*R.i()*s.row(0).t());
  lR2[0] = (K+4)*sigma*log(1+d[0]/2);
  lR3[0] =lR2[0]+lR1[0];
  Pt = sqrt(diagmat(Qt));
  
  for(int t = 1; t<T; t++){
    if (1<(as<double>(qchisq(0.9973,K))/d[t-1])){ 
      Qt = (1-par[0]-par[1])*Qb+ par[0]*CN*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    } else {
      Qt = (1-par[0]-par[1])*Qb + par[0]*CN*K/d[t-1]*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    }
    Pt = sqrt(diagmat(Qt));
    iPt = Pt.i();
    R =  iPt*Qt*iPt;
    lR1[t] = log(det(R));
    d[t] = arma::conv_to<double>::from(s.row(t)*R.i()*s.row(t).t());
    lR2[t] = (K+4)*sigma*log(1+d[t]/2); 
    lR3[t] = lR2[t]+lR1[t];
  }
  return Rcpp::wrap(mean(lR3));
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @noRd
//' @useDynLib Robpvc
//' @import Rcpp RcppArmadillo
// [[Rcpp::export]]
SEXP cor_cDCC(arma::vec par,arma::mat Qb,arma::mat s){
  int T = s.n_rows;
  int K = s.n_cols;
  Rcpp::Function pchisq("pchisq");
  Rcpp::Function qchisq("qchisq");
  arma::vec d(T+1);
  arma::mat R(K,K),Qt(K,K),Pt(K,K), iPt(K,K);
  arma::cube Rpred(K,K,T+1);
  double CN = K/(K*as<double>(pchisq(qchisq(0.9973,K),K+2)) + as<double>(qchisq(0.9973,K))*(1-0.9973));
  R = Qb;
  Qt = Qb;
  d[0] = arma::conv_to<double>::from(s.row(0)*inv_sympd(R)*s.row(0).t());
  Pt = sqrt(diagmat(Qt));
  Rpred.slice(0) = R;
  for(int t = 1; t<T+1; t++){
    if (1<(as<double>(qchisq(0.9973,K))/d[t-1])){ 
      Qt = (1-par[0]-par[1])*Qb+ par[0]*CN*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    } else {
      Qt = (1-par[0]-par[1])*Qb + par[0]*CN*K/d[t-1]*Pt*s.row(t-1).t()*s.row(t-1)*Pt +par[1]*Qt;
    }
    Pt = sqrt(diagmat(Qt));
    iPt = Pt.i();
    R =  iPt*Qt*iPt;
    Rpred.slice(t) = R;
    if(t<T){
      d[t] = arma::conv_to<double>::from(s.row(t)*R.i()*s.row(t).t());
    }
  } 
  return Rcpp::wrap(Rpred);
}

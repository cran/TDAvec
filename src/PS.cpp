#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computePS(NumericMatrix D, int homDim, NumericVector scaleSeq,int p){
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      ++n_rows; 
    }
  }
  
  int L = scaleSeq.size();
  if (n_rows == 0) return NumericVector(L);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
  
  NumericVector pp = pow(y-x,p);
  NumericVector w = pp/sum(pp);
  
  NumericVector phi(L);
  for (int i=0;i<L;++i){
    NumericVector Lambda = pmax(pmin(scaleSeq[i] - x, y - scaleSeq[i]),0);
    phi[i] = sum(w*Lambda);
  }
  return phi; 
}


#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix BallDistanceVector(Rcpp::List Y)
  /* input a list of vector
   * return distance matrix for Ball Impurity
   * use L2 norm for vector distance
   */
{
  int n,i,j;
  n=Y.size(); //number of subjects
  Rcpp::NumericMatrix D(n,n);
  for(i=0;i<n;i++)
  {
    for(j=i+1;j<n;j++)
    {
      arma::vec tempi=Rcpp::as<arma::vec>(Y[i]);
      arma::vec tempj=Rcpp::as<arma::vec>(Y[j]);
      arma::vec diff=tempi-tempj;
      D(i,j)=arma::norm(diff,2);
      D(j,i)=D(i,j);
    }
  }
  return(D);
}



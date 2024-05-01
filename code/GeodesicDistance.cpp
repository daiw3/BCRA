#include <RcppArmadillo.h>
#include <stdio.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix vector_to_matrix(Rcpp::NumericVector values, Rcpp::NumericVector dimensions ) {
  
  // Create an empty matrix with dimensions "dimensions"
  Rcpp::NumericMatrix mat(dimensions[0], dimensions[1]);
  
  // Fill matrix "mat" by values of vector "values"
  std::copy(values.begin(), values.end(), mat.begin());  
  
  return mat;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix GeodesicDistance(Rcpp::List Y, 
                                     double epsil)
  /* input a list of matrices
   * return distance matrix for Ball Impurity
   * use geodesic distances
   */
{
  int n,i,j;
  n=Y.size(); //number of subjects
  
  Rcpp::NumericMatrix D(n,n);
  for(i=0;i<n;i++)
  {
    for(j=i+1;j<n;j++)
    {
      arma::mat tempi =  as<arma::mat>(Y[i]);
      arma::mat tempj = as<arma::mat>(Y[j]);
      
      // convert into correlation matrices
      arma::mat FC1 = arma::tanh(tempi);
      arma::mat FC2 = arma::tanh(tempj);
      
      // compute Q_1^{-1/2} via eigen value decmposition
      arma::mat U, V;
      arma::vec S;
      arma::svd(U, S, V, FC1, "dc");
      
      // lift very small eigen values
      int ii;
      int p = S.size(); //number of eigen_values
      for(ii=0;ii<p;ii++){
        if(S[ii] < epsil){
          S[ii] = epsil;
        }
      }
      arma::vec S_sqrt = arma::sqrt(S);
      // FC1^{-1/2} = u[s^{-1/2}]u'
      Rcpp::NumericMatrix S_mat(p, p);   
      int ptr;
      for(ptr=0;ptr<p;ptr++){
        S_mat(ptr, ptr) = 1.0/S_sqrt[ptr];       
      }
      arma::mat Stmp =  as<arma::mat>(S_mat);
      arma::mat FC1_mod = U * Stmp * U.t();
      arma::mat M = FC1_mod * FC2 * FC1_mod;
      
      arma::mat U2, V2;
      arma::vec S2;
      arma::svd(U2, S2, V2, M, "dc");
      
      double stmp = arma::sum(arma::square(arma::log(S2)));
      
      D(i,j)= sqrt(stmp);
      D(j,i)=D(i,j);
    }
  }
  return(D);
}

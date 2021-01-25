// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

//' Gradient for OT (Internal C function)
//'
//' @param Mat the cost matrix
//' @param psi the vector of potentials
//' @export
//[[Rcpp::export]]
arma::vec Grad(NumericMatrix Mat, arma::vec psi){

  unsigned int i,p=Mat.ncol(), q=Mat.nrow();
   mat X = mat(Mat.begin(), Mat.nrow(), p, true); //true is there to ensure that it is copied
 // Rcout << "The value is " << X << std::endl;
  X.each_col() -= psi ;
  //Rcout << "The value is " << X << std::endl;
  vec F(q);
F.zeros();
  //Rcout << "The value is " << p << std::endl;
    for(i=0;i<p;++i)
      F[(X.col(i)).index_min()]+=1;
    for (i=0;i<q;++i)
    {
      F[i]/=p;
      F[i]-=(double) 1/q;
    }

  return F;
}

//' Objective for OT (Internal C function)
//'
//' @param Mat the cost matrix
//' @param psi the vector of potentials
//' @export
//[[Rcpp::export]]
double Objective(NumericMatrix Mat, arma::vec psi ){
  unsigned int i,p=Mat.ncol(), q=Mat.nrow();
  mat X = mat(Mat.begin(), Mat.nrow(), p, true);
  X.each_col() -= psi ;
  vec  F(q);
  F.zeros();

// Rcout << "The value is " <<  F << std::endl;
  int tempor=0;
  for(i=0;i<p;++i)
  {

    tempor=(X.col(i)).index_min();
  //  Rcout << "The value is " << tempor << std::endl;
    F[tempor]+= (X.col(i))[tempor];
   //Rcout << "The value is " <<  F[tempor] << std::endl;
  }

  //Rcout << "The value is " << F << std::endl;
  F /= p;

 // Rcout << "The value is " << F << std::endl;
    F += psi/q;

  return -sum(F);
}

//' SGD Algorithm for OT (Internal C function)
//'
//' @param C the descent constant
//' @param Mat the cost matrix
//' @export
//[[Rcpp::export]]
NumericVector SGD_OT(double C, const NumericMatrix  Mat){
  int N=Mat.ncol();
  int q=Mat.nrow() ;
  int i;
  NumericVector tildev(q), v(q), tempor(q),tempor2(q);
  for( i = 0; i < N; i++)
  {
   tempor=Mat(_,i);
   tempor=tempor-tildev;
    tempor2.fill((double) -1/q);
    tempor2[which_min(tempor)]+= 1;
    tildev = tildev -  C/sqrt(i+1)*tempor2 ;
    v = tildev/(i+1) + (double) i/(i+1)*v;

  }
  return v;
}


//'  Algorithm for plotting results of OT (Internal C function)
//'
//' @param Mat the cost matrix
//' @param psi the vector of potentials
//' @export
//[[Rcpp::export]]
NumericVector Plot(NumericMatrix Mat, arma::vec psi){
  unsigned int i,p=Mat.ncol();
  mat X = mat(Mat.begin(), Mat.nrow(), p, true); //true is there to ensure that it is copied
  // Rcout << "The value is " << X << std::endl;
  X.each_col() -= psi ;
  //Rcout << "The value is " << X << std::endl;
  NumericVector F(p);
  for(i=0;i<p;++i)
    F[i]=(X.col(i)).index_min();
  return F;
}

 // [[Rcpp::export]]
 arma::mat rGausResid( int n, int dim) {
   arma::mat Y = arma::randn(n, dim);
   Y= Y- repmat(mean(Y,0),n,1);
   arma::mat Sigma = Y.t()*Y/(n-1); // n-1 to copy the R function cov
   Y = Y * solve(chol(Sigma), eye(dim,dim));
   return  Y ;
 }










// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//



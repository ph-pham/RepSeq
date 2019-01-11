#include <Rcpp.h>

using namespace Rcpp;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
//' compute muScore
//'
//' This function returns multivariate score for each row of the input mumerical data set using the O'Brien (1984)'s method.
//' @param df a data frame of counts with samples in columns and clonotypes in rows
//' @return a vector of multivariate scores computed for each row 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector muScoreOB (const Rcpp::DataFrame & df){
  NumericMatrix x = internal::convert_using_rfunction(df, "as.matrix");
  unsigned int outrows = x.nrow(), i = 0, j = 0;
  double d;
  Rcpp::NumericVector out(outrows);
  Rcpp::NumericMatrix tmp(outrows, outrows);
  for (i = 0; i < outrows - 1; i++){
    Rcpp::NumericVector v1 = x.row(i);
    for (j = i + 1; j < outrows; j++){
      d = sum( (v1 > x.row(j)) - (x.row(j) > v1) );
      tmp(i, j) = sgn(d);
      tmp(j, i) = -sgn(d);
    }
  }
  out = rowSums(tmp);
  return wrap(out);
}


//' compute pre-muScore matrix
//' This function returns a matrix of scores comparing pairwise rows for each sample of the input data set. 
//'
//' @param df a data frame of counts with samples in columns and clonotypes in rows. 
//' @return a matrix  
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix muScoreM (const Rcpp::DataFrame & df){
  NumericMatrix x = internal::convert_using_rfunction(df, "as.matrix");
  unsigned int outrows = x.nrow(), i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(outrows, outrows);
  for (i = 0; i < outrows - 1; i++){
    Rcpp::NumericVector v1 = x.row(i);
    for (j = i + 1; j < outrows; j++){
      d = sum( (v1 > x.row(j)) - (x.row(j) > v1) );
      out(i, j) = sgn(d);
      out(j, i) = -sgn(d);
    }
  }
  return wrap(out);
}
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix surface(NumericMatrix grid) {
  //Rcout << "The value is " << grid.nrow() << std::endl;

  for (int i = 0; i < thresholds.nrow(); i++) {
    
  };
}

/*
// [[Rcpp::export]]
NumericMatrix surface(NumericMatrix grid,
                      NumericMatrix thresholds,
                      NumericVector x1,
                      NumericVector x2,
                      NumericVector g1,
                      NumericVector g2,
                      NumericVector y,
                      NumericVector weight) {
  
  //for (int i = 0; i < thresholds.);
}
*/

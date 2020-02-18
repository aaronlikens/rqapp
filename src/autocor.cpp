// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

double acf(arma::colvec& x, int& lag, double& meanX);


//TODO:  explore whether putting things together in the same function
//      will improve performance even more.  Right now, this is not necessary
//      as teh function is already faster than R by an order of magnitude.  
//      I achieve this speed by skipping time-consuming error checking,s etc...

// [[Rcpp::export]]
arma::colvec autocor(arma::colvec& x, int& lag_max){
    arma::colvec acor(lag_max);
    double meanX = mean(x);
    for (int i = 0; i < lag_max; ++i){
        acor(i) = acf(x,i,meanX);
    }
    
    return acor/max(acor);
}

// compute autocorrelation for a given lag.
double acf(arma::colvec& x, int& lag, double & meanX){
    double autocv; 
    
    // Loop to compute autovariance
    autocv = 0.0;
    for (unsigned int i=0; i<(x.size() - lag); i++){
        
        autocv += ((x[i] - meanX) * (x[i+lag] - meanX)) / (x.size());
    }
    
    return autocv;
}
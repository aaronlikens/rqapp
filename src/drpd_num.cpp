// #define ARMA_NO_DEBUG //decomment to speed up code once tested
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' dprd_numeric
//' 
//' Cross-Recurrence diagonal profile of two continuous time-series
//' 
//' Quick method to explore the cross-recurrence diagonal profile of two
//' time series. It returns the recurrence observed for different delays, 
//' the maximal recurrence obsered and teh dealy at which it occured.
//' 
//' @param t1 a time series
//' @param t2 a time series
//' @ws a constant indicating the range of delays (positive and negative) to explore
//' @radius
//' @export
// [[Rcpp::export]]
List drpd_numeric(arma::vec t1, arma::vec t2, int ws, double radius){
    
    unsigned int drpd_length = 2*ws + 1;
    unsigned int counter = 0;
    arma::vec drpd((drpd_length));
    drpd.zeros();

    // lag in one direction
    for (int i = (-ws-1); i < -1; ++i){
        int ix = abs(i);
        arma::vec y = t2(arma::span(ix-1, t2.n_elem-1));
        arma::vec x = t1(arma::span(0, y.n_elem-1));
        
        arma::vec dif = arma::abs(x-y);
        arma::vec in_radius = dif(arma::find(dif <= radius));
        double recpoint = in_radius.n_elem;
        drpd(counter) = recpoint/y.n_elem;
        counter += 1;
    }

    // zero lag
    arma::vec dif = abs(t1-t2);
    arma::vec in_radius = dif(arma::find(dif <= radius));
    double recpoint = in_radius.n_elem;
    drpd(counter) = recpoint/t1.n_elem;
    counter +=1; 

    // lag in the other direction
    for (int i = 1; i <= ws; ++i){

        arma::vec x = t1(arma::span(i, t1.n_elem-1));
        arma::vec y = t2(arma::span(0, x.n_elem-1));
        
        arma::vec dif = arma::abs(x-y);
        arma::vec in_radius = dif(arma::find(dif <= radius));
        double recpoint = in_radius.n_elem;
        drpd(counter) = recpoint/y.n_elem;
        counter += 1;
    }
    
    double maxrec = arma::max(drpd);
    arma::uvec maxlag = find(drpd == maxrec);
    
    return List::create(Named("profile") = drpd,
                        Named("maxrec") = maxrec,
                        Named("maxlag") = maxlag + 1
                            );
}


/*** R
library(crqa)
set.seed(1239)
x = rnorm(1000)
set.seed(1234)
y = rnorm(1000)

cpp = drpd_numeric(x,y,10,.01)
crq = drpdfromts(x,y,10,'continuous',.01)
*/


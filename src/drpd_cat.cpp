#include <Rcpp.h>
using namespace Rcpp;

// simple function for creating a vector of integers for use in indexing
IntegerVector seq_int(int from, int to);

// simple function to count the true values in a logical vector
double count_if(LogicalVector x);

const double tol = 0.0001;

// [[Rcpp::export]]
List drpd_cat(StringVector t1, StringVector t2, int ws) {
    unsigned int drpd_length = 2*ws + 1;
    unsigned int counter = 0;
    NumericVector drpd((drpd_length));
    
    // lag in one direction and then lag the other
    for (int i = (-ws-1); i < -1; ++i){
        int ix = abs(i);
        StringVector y = t2[seq_int(ix-1, t2.size()-1)];
        StringVector x = t1[ seq_int(0, y.size()-1)];
        drpd(counter) = (count_if(y==x))/y.size();
        counter += 1;
    }

    // compute recurrence at zero lag
    drpd(counter) = count_if(t1==t2)/t1.size();
    counter += 1;
    
    // lag in other direction
    for (int i = 1; i <= ws; ++i){
        StringVector x = t1[seq_int(i, t1.size()-1)];
        StringVector y = t2[seq_int(0, x.size()-1)];
        drpd(counter) = (count_if(y==x))/y.size();
        counter += 1;
    }
    
    // find maximum recurrence 
    double maxrec = max(drpd);
    
    // find lag of maximum recurrence
    unsigned int maxlag = which_max(drpd);
    
    // collect output into a list and return
    return List::create(Named("profile") = drpd,
                        Named("maxrec") = maxrec,
                        Named("maxlag") = maxlag);
}

// [[Rcpp::export]]
IntegerVector seq_int(int from, int to){
    
    IntegerVector sequence(abs(from-to)+1);
    sequence[0] = from;
    for (int i = 1; i < sequence.size(); ++i){
        sequence[i] = sequence[i-1]+1;
    }
    
    return sequence;
    
}

// [[Rcpp::export]]
double count_if(LogicalVector x) {
    int counter = 0;
    for(int i = 0; i < x.size(); i++) {
        if(x[i] == TRUE) {
            counter++;
        }
    }
    return counter;
}



/*** R
library(crqa)
set.seed(1234)
x = sample(letters,100000, replace = TRUE)
set.seed(1234)
y = sample(letters,100000, replace = TRUE)
drpd_cat(x,y,10)
drpdfromts(x,y,ws = 10, datatype = 'categorical')

*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//x     : time series
//tao   : time delay
//mmax  : maximum embedding dimension
//rtol  : near tolerance
//atol  : far tolerance, to accomdate for the fact that all neighbors may be far away
//        to begin with.

//reference:M. B. Kennel, R. Brown, and H. D. I. Abarbanel, Determining
//embedding dimension for phase-space reconstruction using a geometrical 
//construction, Phys. Rev. A 45, 3403 (1992). 
//author:"Merve Kizilkaya"

// TODO: This is pretty slow right now with several dim. That said, it is still several
//      times faster than base R. I think nonlinearTseries uses a better neighbor search algorithm
//      than this brute force approach.  Perhaps we shoud use one of those.
// TODO: another option to increase speed would be to examine only a random sample of
//       data points.  This is also done in Marwan's code in order to incre ase efficiency

// population standard deviation
double sdp(arma::vec x){
     
    return std::sqrt(arma::accu(arma::pow(x-arma::mean(x),2))/x.n_elem);
}


// do phase space reconstruction
arma::mat psr(arma::vec x, int m, int tao, int npoints){
    // int N = x.n_elem;
    // int M = x.n_elem - m*tao;
    int M = npoints;
    arma::mat Y(M,m);
    Y.zeros();
    
    for ( int i = 0; i < m; i++){
        int start = i*tao;
        // Rcout << "start is " << start << "\n";
        int stop = M-1 + (i*tao);
        // Rcout << "stop is " << stop << "\n";
    
        Y.col(i) = x.subvec(start,stop);
    }
    
    
    return Y;
    
}


// [[Rcpp::export]]
arma::vec fnn(arma::vec x, int tao, int mmax, double rtol, double atol){
    
    double Ra = sdp(x);
    arma::vec fn(mmax+1);
    fn.zeros();
    double D = 0.0;
    double R = 0.0;
    
    for (int m = 1; m < mmax + 1; ++m){
        // perform phase space reconstruction
        int M = x.n_elem - m*tao;
        arma::mat Y = psr(x,m,tao,M);
        for (int n = 0; n < M; ++n){
            
            arma::colvec ones(M);
            ones.ones();
            // create a matrix, y0(M,m), where each row is equal to 
            // row n of of the psr, Y
            arma::mat y0 = ones*Y.row(n);
            
            // substract y0 from Y to obtain the distance of row n from 
            // all other rows
            arma::vec distance = arma::sqrt(arma::sum(arma::pow(Y-y0,2),1));
            
            //TODO: MAYBE USE A FIND FUNCTION INSTEAD OF SORT
            //      NOTE: I TRIED SOMETHING LIKE THIS IN THE DEV VERSION BUT
            //      IT WAS SLOWER THAN THIS ONE.
            // arma::vec neardis = sort(distance);
            arma::uvec nearpos = sort_index(distance);
            
            // only examine the rth nearest neighbor where r = 1 (i.e. Kennel et al.,1992)
            unsigned int ix1 = n + (m)*tao;
            unsigned int ix2 = nearpos(1) + (m)*tao;
            D = std::abs(x(ix1) - x(ix2));
            R = std::sqrt(std::pow(D,2) + std::pow(distance(nearpos(1)),2));
            
            // comapre to tolerances to determine if the a false neighbor
            if (D/distance(nearpos(1)) > rtol || R/Ra > atol){
                
                fn(m) += 1;
            }
        }
    }
    
    // discard first element in fn because it is just a zero
    arma::vec result = fn.subvec(1,mmax);
    
    return (result/result(0))*100;
}



/***R
x = tseriesChaos::lorenz.ts
test = fnn(x,11,6,10,2)
print(test)
*/







































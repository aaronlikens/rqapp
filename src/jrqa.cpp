#define ARMA_NO_DEBUG 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// function for creating embedding matrix (phase-space reconstruction)
arma::mat jpsr(arma::vec& x, int dim, int tau);

// function to compute distance matrix
void jdist_mat(arma::mat& x1, arma::mat& x2, arma::mat& dist);

// function for creating the recurrence plot
void j_rp(arma::mat& x, double radius);

// function for computing the variance metrics of recurrence
List jline_stats(arma::mat& rp, int mindiagline, int minvertline, int t_win);

// extract diagonal line lengths
arma::vec jdiagonal_lines(arma::mat& rp, int mindiagline);

// extract vertical line lengths
arma::vec jvertical_lines(arma::mat& rp, int minvertline);


const double tol = 0.0001;


// [[Rcpp::export]]
List jrqa(arma::mat& ts, unsigned int embed = 1, unsigned int delay = 1,
          int normalize = 1, int rescale = 1, int mindiagline = 2,
          int minvertline = 2, int t_win = 1, double radius = 0.0001, 
          int whiteline = 0, int recpt = 0){
    
    
    // dm.zeros();
    arma::mat jrp(ts.n_rows, ts.n_rows);
    jrp.zeros();
    for (unsigned int i = 0; i < ts.n_cols; ++i){
        // normalize 
        arma::vec temp = ts.col(i);
        
        if (normalize > 0){
            switch(normalize){
            // normalize to unit interval
            case 1: temp = (temp - min(temp))/max(temp);
                
                // normlaize to a z-score
            case 2: temp = (temp - mean(temp))/stddev(temp);
                
            }
        }
        
        // embed time series according to chosen dimension and time delay
        arma::mat jpsr_ts = jpsr(temp, embed, delay);
        arma::mat dm(ts.n_rows, ts.n_rows);
        jdist_mat(jpsr_ts, jpsr_ts, dm);
       
        // arma::mat dm_rescale(dm.n_rows, dm.n_cols);
        
        double rescaledist;
        
        // find indices of the distance matrix that fall within the prescriped radius
        if (rescale > 0){
            switch(rescale){
            
            // rescale distance matrix to mean distance
            case 1: rescaledist = mean(mean(dm));
                dm = (dm/rescaledist)/100;
                
                // rescale distance matrix to max distance
            case 2: rescaledist = max(max(dm));
                dm = (dm/rescaledist)/100;
            }
        } else{
            dm = dm.t();
        }
        
        // create the recurrence plot matrix
        j_rp(dm, radius);
        
        if (i == 0){
            jrp = dm;
        }else{
            jrp = jrp % dm;
        }
        // jrp = dm;
    }
    
    List output = jline_stats(jrp, mindiagline, minvertline, t_win);

    
    // comute recurrence matrix
    if (recpt > 0){
        return List::create(Rcpp::Named("rp") = jrp, 
                            Rcpp::Named("rqa") = output);
    }else{
        return List::create(Rcpp::Named("rqa") = output);
    }
    
    
}



arma::mat jpsr(arma::vec& x, int dim, int tau){
    unsigned int embed = dim;
    unsigned int delay = tau;
    
    // find size of jpsr matrix and allocate memory
    arma::mat space(x.n_elem - (embed-1)*delay, embed);
    space.zeros();
    
    // insert first column 
    space.col(0) = x(arma::span(0, x.n_elem - embed*delay + delay -1));

    // main loop that does the embedding
    unsigned int start = 0;
    unsigned int stop = x.n_elem - embed*delay + delay - 1;
    for (unsigned int i = 1; i < embed; ++i){
        start += delay;
        stop  += delay;
        space.col(i) = x(arma::span(start,stop));
    }
    return space;
    
}

// TODO: This function could be additionally optomized
//      Currently  it is aobut 2-3 time slower than rdist in the fields package


void jdist_mat(arma::mat& x1, arma::mat& x2, arma::mat& dist){
    // TODO: distance matrix and recurrence plot could be combined in one for loop
    //      CAN'T BE DONE:  ENTIRE DISTANCE MATRIX MUST BE KNOWN BEFORE RESCALING, 
    //      IMPLING THAT RECURRENCE PLOT MUST BE CALCULATED AFTER THE FAC
    // TODO: add in check for matrix and throw error
    // TODO: also check that matrices have same dimension
    
    //NumericVector tempx1(x1.cols());
    // arma::mat dist(x1.n_rows, x2.n_rows);
    
    for ( unsigned int i = 0; i < x1.n_rows; ++i ){
        
        for (unsigned int j = 0; j < x1.n_rows; ++j){

            // dist(i,j) = sqrt(arma::accu(pow(x2.row(j) - x1.row(i),2)));
            dist(i,j) = sqrt(arma::accu(
                (x2.row(j) - x1.row(i))%(x2.row(j) - x1.row(i))
                                             ));
        }
        
    }
    dist.t();
    
    // return(dist);
}


void j_rp(arma::mat& x, double radius){
    
    // arma::mat recurrence_plot = as<arma::mat>(x);
    // find indices outside of specified range and flag as 999
    arma::uvec index1 = find(x > radius);
    // x.elem( index1 ).fill(999);
    
    // find all values inside the radius and replace with 1
    // arma::uvec index2 = find(x <= radius);
    // x.elem( index2 ).fill(1);
    x.elem( find(x <= radius) ).fill(1);
    // replace all 999s with 0s
    x.elem( index1 ).fill(0);

    // return x;
}




List jline_stats(arma::mat & rp, int mindiagline, int minvertline, int t_win){

    double rr = 0.0, det = 0.0, maxline = 0.0, meanline = 0.0, entropy = 0.0, rentropy = 0.0;
    double lam = 0.0, tt = 0.0, vmax = 0.0;
    arma::uword nrline;
    nrline = 0;
    double rp_sum = arma::accu(rp);
    
    // TODO: DONE: I need to make theiler window adjustments here
    //      In addition, I will follow Marwan and blank account recurrence on the
    //      diagonal. i.e. zero out the line of incidence before computing stats
    // TODO: For cross-recurrence, I don't wan to blank out LOI
    
    if (t_win > 0){
        for (int i = -t_win; i < t_win; ++i){
            rp.diag(i).zeros();
        }
        // for (int i = 0; i < t_win; ++i){
        //     rp.diag(-i).zeros();
        // }
    }else{
        rp.diag().zeros();
    }
    
    // calculate percent recurrence
    // double nrecurrence = arma::accu(rp);
    // if (arma::accu(rp) > 0){
    if(rp_sum > tol){
        // rr = 100*nrecurrence/(rp.n_rows*rp.n_rows);
        // rr = 100*arma::accu(rp)/(rp.n_rows*rp.n_rows);
        rr = 100*rp_sum/(rp.n_rows*rp.n_rows);
    }
    
    // create initialize vector to store lengths
    // int ndiags = x.rows();
    arma::vec dlengths = jdiagonal_lines(rp, mindiagline);
    arma::vec vlengths = jvertical_lines(rp, minvertline);

    
    // TODo: make decisions about upper, lower, both triangles. 
    //      For now, I will do both
    arma::vec jdiagonal_lines = dlengths(find(dlengths > tol));
    arma::vec vert_lines = vlengths(find(vlengths > tol));
    
    // estiamte number of lines parameter
    nrline = jdiagonal_lines.n_elem;
    // arma::vec u = (unique(jdiagonal_lines));
    // double n_unique = u.n_elem;
    
    
    
    // average length of lines
    if(jdiagonal_lines.n_elem > 0){
        meanline = mean(jdiagonal_lines);
        
        // maximum line length
        maxline = jdiagonal_lines.max();
        
        // perent determinism.  divide length by 29
        det = (arma::accu(jdiagonal_lines)/arma::accu(rp))*100;
        
        //TODO: entropy and normalized entropy.  These functions are not yet completely 
        //      equivalen those found in the CRQA package by Moreno.  In fact, I'm not sure
        //      that the calculations in Morenoa are correct outside of categorical recurrence
        //      NOTE: Entropy value seem to match Marwan.
        
        // calculate entropy and relative entropy
        arma::vec counts = arma::conv_to<arma::vec>::from(hist(jdiagonal_lines,arma::linspace(1,rp.n_rows,rp.n_rows)));
        arma::vec prob = counts/sum(counts);
        arma::vec p = prob(find(prob > 0));
        entropy = -arma::accu(p % arma::log(p));
    
        
        // relative entropy following form Moreno and Dale
        double nbins = p.n_elem;
        double denom = -1*std::log(1/nbins);
        rentropy = entropy/denom;
    }
    
    if(vert_lines.n_elem > 0){
        lam = (sum(vert_lines)/arma::accu(rp))*100;
        tt = mean(vert_lines);
        vmax = max(vert_lines);
    }
    
    // TODO: ADD IN OTHER MEASURES FROM MARWAN E.G., T1, T2, RTE, CLUST, TRANS
    return List::create(Named("rr") = rr,
                        Named("det") = det,
                        Named("div") = 1/maxline,
                        Named("nrline") = nrline,
                        Named("ratio") = det/rr,
                        Named("maxline") = maxline,
                        Named("meanline") = meanline,
                        Named("lam") = lam,
                        Named("tt") = tt,
                        Named("vmax") = vmax,
                        Named("entropy") = entropy,
                        Named("rentropy") = rentropy);
}



arma::vec jdiagonal_lines(arma::mat& rp, int mindiagline){
    
    arma::vec lengths(1);
    lengths(0) = 0;
    
    // start interation with first diagonal
    // extract diagonals one by and and find line lengths
    int start = rp.n_rows;
    int stop = rp.n_rows;
    start *= -1;
    start += 1;
    
    for (int i = start; i < stop; ++i){
        
        //  check if vertical line contains any recurrent points
        // if not, skip to the next line and avoid those computations
        if (arma::accu(rp.diag(i)) > 0.0001){
            
            // zero pad to capture patterns that start and end with 1
            arma::vec pad1(10);
            arma::vec d = arma::join_cols(arma::join_cols(pad1, rp.diag(i)),pad1);
            
            // find all zeros
            arma::uvec inds = find(d < 1);
            arma::vec k = arma::conv_to<arma::vec>::from(inds);
            
            // take the first difference of the indices of zeros so that 
            // difference between them reveals the length of interceding 
            // stretches of 1s
            arma::vec diffs = arma::diff(k)-1;
            
            // get all lines longer than minvertline
            arma::vec diaglines = diffs.elem(find(diffs >= mindiagline));
            
            
            // grow vector of lines lengths on every loop
            lengths = arma::join_vert(lengths, diaglines);
            
        }
        
    }
    
    
    
    return lengths;
}



arma::vec jvertical_lines(arma::mat& rp, int minvertline){
    
    arma::vec lengths(1);
    lengths(0) = 0;
    // start interation with first column
    // extract columns one by and and compute
    int n_cols = rp.n_cols;
    for (int i = 0; i < n_cols; ++i){
        
        //  check if vertical line contains any recurrent points
        // if not, skip to the next line and avoid those computations
        if (arma::accu(rp.col(i)) > 0.0001){
            // zero pad to capture patterns that start and end with 1
            // TODO: TRY TO ELIMINATE SOME UNWANTED STEPS HERE
            arma::vec pad1(10);
            // arma::vec pad2(10);
            // arma::vec d2 = rp.col(i);
            arma::vec d3 = arma::join_cols(pad1, rp.col(i));
            // arma::vec d = arma::join_cols(d3, pad2);
            arma::vec d = arma::join_cols(d3,pad1);
            
            // find all zeros
            arma::uvec inds = find(d < 1);
            
            arma::vec k = arma::conv_to<arma::vec>::from(inds);
            
            // take the first difference of the indices of zeros so that 
            // difference between them reveals the length of interceding 
            // stretches of 1s
            arma::vec diffs = arma::diff(k)-1;
            
            // get all lines longer than minvertline
            arma::vec vertlines = diffs.elem(find(diffs >= minvertline));
            
            
            // grow vector of lines lengths on every loop
            lengths = arma::join_vert(lengths, vertlines);
            
        }
        
    }
    
    
    
    return lengths;
}




/*** R
x = rep(0:2, 300)
x = cbind(x, x, x)

test = jrqa(ts = x, embed = 1, delay = 1, normalize = 1, rescale = 1,
                      mindiagline = 2, minvertline = 2, t_win = 1, radius = .00001, 
                      recpt = 0)

print(test$rqa)

x = sample(0:2, 300, replace = TRUE)
x = cbind(x,x,x)
test = jrqa(ts = x, embed = 1, delay = 1, normalize = 1, rescale = 1,
               mindiagline = 2, minvertline = 2, t_win = 1, radius = .00001, 
               recpt = 0)

print(test$rqa)
*/














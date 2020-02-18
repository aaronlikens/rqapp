#define ARMA_NO_DEBUG //decomment to speed up code once tested
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// TODO: currently, the reuslts do not match those from Wallot; however
//       because our other results match Marwan as well as Coco and Dale,
//      I think there must be differences in how Wallot is calculating the 
//      recurrence metrics. Therefore, I will default to Marwan.

// function for creating embedding matrix (phase-space reconstruction)
arma::mat mdpsr(arma::mat& x, int dim, int tau);

// function to compute distance matrix
void md_dist_mat(arma::mat& x1, arma::mat& x2, arma::mat& dist);

// function for creating the recurrence plot
void md_rp(arma::mat& x, double radius);

// function for computing the variance metrics of recurrence
List md_line_stats(arma::mat& md_rp, int mindiagline, int minvertline, int t_win);

// extract diagonal line lengths
arma::vec md_diagonal_lines(arma::mat& md_rp, int mindiagline);

// extract vertical line lengths
arma::vec md_vertical_lines(arma::mat& md_rp, int minvertline);


const double tol = 0.0001;


// [[Rcpp::export]]
List mdrqa(arma::mat& ts, unsigned int embed = 1, unsigned int delay = 1,
           int normalize=1, int rescale = 1, int mindiagline = 2, 
           int minvertline = 2, int t_win = 1, double radius=0.0001, 
           int recpt = 0){
    

    
    // normalize 
    if (normalize > 0){
        switch(normalize){
        //TODO: SEE HOW WALLOT NORMALIZES I.E., BY COLUMN OR THE ENTIRE MATRIX
        //      SEEMS LIKE EACH COLUMN SHOULD BE NORMALIZED ON ITS OWN
        
        // normalize to unit interval
        case 1: 
            for (unsigned int i = 0; i < ts.n_cols; i++){
                ts.col(i) = (ts.col(i)- arma::min(ts.col(i)))/arma::max(ts.col(i));
            }
            //ts = (ts - min(ts))/max(ts);
        // normlaize to a z-score
        case 2: 
            for (unsigned int i = 0; i < ts.n_cols; i++){
                ts.col(i) = (ts.col(i)- arma::mean(ts.col(i)))/arma::stddev(ts.col(i));
            }
            //ts = (ts - mean(ts))/stddev(ts);
        }
    }
    
    // embed time series according to chosen dimension and time delay
    //TODO: MAY NEED TO DECLARE THESE MATRICES OUTSIDE OF THE IF STATMENT BECAUSE OF SCOPE
    arma::mat dm(ts.n_rows,ts.n_rows);
    dm.zeros();
    md_dist_mat(ts,ts,dm);
    // arma::mat mdpsr_ts;
    if(embed > 1){
        arma::mat mdpsr_ts = mdpsr(ts, embed, delay);
        md_dist_mat(mdpsr_ts, mdpsr_ts, dm);

    }else{
        md_dist_mat(ts,ts,dm);
    }
    
    
    arma::mat dm_rescale(dm.n_rows, dm.n_cols);
    
    double rescaledist;
        
    // find indices of the distance matrix that fall within the prescriped radius
    if (rescale > 0){
        switch(rescale){
        
        // rescale distance matrix to mean distance
        case 1: rescaledist = mean(mean(dm));
                dm_rescale = (dm/rescaledist)/100;
        
        // rescale distance matrix to max distance
        case 2: rescaledist = max(max(dm));
                dm_rescale = (dm/rescaledist)/100;
        
        // rescale distance matrix to min distance
        case 3: rescaledist = min(min(dm));
                dm_rescale = (dm/rescaledist)/100;
        }
    } else{
        dm_rescale = dm;
    }
    
    // create the recurrence plot matrix
    md_rp(dm, radius);
    List output = md_line_stats(dm, mindiagline, minvertline, t_win);

    
    // comute recurrence matrix
    if (recpt > 0){
        return List::create(Rcpp::Named("rp") = dm, 
                            Rcpp::Named("rqa") = output);
    }else{
        return List::create(Rcpp::Named("rqa") = output);
    }
    
}


// [[Rcpp::export]]
arma::mat mdpsr(arma::mat& x, int dim, int tau){
    unsigned int embed = dim;
    unsigned int delay = tau;

    // find size of mdpsr matrix and allocate memory
    arma::mat space(x.n_rows - (embed-1)*delay, embed*x.n_cols);
    space.zeros();

    // insert first column
    arma::mat temp = x.submat(0, 0, x.n_rows - embed*delay + delay -1, x.n_cols-1);

    space.submat(0, 0, space.n_rows-1, x.n_cols-1) =
        x.submat(0, 0, x.n_rows - embed*delay + delay -1, x.n_cols-1);

    // main loop that does the embedding
    unsigned int rowstart = 0;
    unsigned int rowstop = x.n_elem - embed*delay + delay - 1;
    unsigned int colstart = 0;
    unsigned int colstop = x.n_cols-1;
    for (unsigned int i = 1; i < embed; ++i){
        rowstart += delay;
        rowstop  += delay;
        colstart += x.n_cols;
        colstop += x.n_cols;
        space.submat(0, colstart, space.n_rows-1, colstop) =
            x.submat(rowstart, 0, rowstop, x.n_cols-1);
    }
    return space;

}

// TODO: This function could be additional optomized
//      Currently  it is aobut 2-3 time slower than rdist in the fields package

void md_dist_mat(arma::mat& x1, arma::mat& x2, arma::mat& dist){
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

void md_rp(arma::mat& x, double radius){
    
    // arma::mat recurrence_plot = as<arma::mat>(x);
    // find indices outside of specified range and flag as 999
    arma::uvec index1 = find(x > radius);
    x.elem( index1 ).fill(999);
    
    // find all values inside the radius and replace with 1
    arma::uvec index2 = find(x <= radius);
    x.elem( index2 ).fill(1);
    
    // replace all 999s with 0s
    x.elem( index1 ).fill(0);

    // return x;
}



List md_line_stats(arma::mat & md_rp, int mindiagline, int minvertline, int t_win){

    double rr = 0.0, det = 0.0, maxline = 0.0, meanline = 0.0, entropy = 0.0, rentropy = 0.0;
    double lam = 0.0, tt = 0.0, vmax = 0.0;
    arma::uword nrline;
    nrline = 0;
    double md_rp_sum = arma::accu(md_rp);
    
    // TODO: DONE: I need to make theiler window adjustments here
    //      In addition, I will follow Marwan and blank account recurrence on the
    //      diagonal. i.e. zero out the line of incidence before computing stats
    // TODO: For cross-recurrence, I don't want to blank out LOI
    
    if (t_win > 0){
        for (int i = -t_win; i < t_win; ++i){
            md_rp.diag(i).zeros();
        }
        for (int i = 0; i < t_win; ++i){
            md_rp.diag(-i).zeros();
        }
    }else{
        md_rp.diag().zeros();
    }
    
    // calculate percent recurrence
    // double nrecurrence = arma::accu(md_rp);
    // if (arma::accu(md_rp) > 0){
    if(md_rp_sum > tol){
        // rr = 100*nrecurrence/(md_rp.n_rows*md_rp.n_rows);
        // rr = 100*arma::accu(md_rp)/(md_rp.n_rows*md_rp.n_rows);
        rr = 100*md_rp_sum/(md_rp.n_rows*md_rp.n_rows);
    }
    
    // create initialize vector to store lengths
    // int ndiags = x.rows();
    arma::vec dlengths = md_diagonal_lines(md_rp, mindiagline);
    arma::vec vlengths = md_vertical_lines(md_rp, minvertline);

    
    // TODo: make decisions about upper, lower, both triangles. 
    //      For now, I will do both
    arma::vec md_diagonal_lines = dlengths(find(dlengths > tol));
    arma::vec vert_lines = vlengths(find(vlengths > tol));
    
    // estiamte number of lines parameter
    nrline = md_diagonal_lines.n_elem;
    // arma::vec u = (unique(md_diagonal_lines));
    // double n_unique = u.n_elem;
    
    
    
    // average length of lines
    if(md_diagonal_lines.n_elem > 0){
        meanline = mean(md_diagonal_lines);
        
        // maximum line length
        maxline = md_diagonal_lines.max();
        
        // perent determinism.  divide length by 29
        det = (sum(md_diagonal_lines)/arma::accu(md_rp))*100;
        
        //TODO: entropy and normalized entropy.  These functions are not yet completely 
        //      equivalen those found in the CRQA package by Moreno.  In fact, I'm not sure
        //      that the calculations in Morenoa are correct outside of categorical recurrence
        //      NOTE: Entropy value seem to match Marwan.
        
        // calculate entropy and relative entropy
        arma::vec counts = arma::conv_to<arma::vec>::from(hist(md_diagonal_lines,arma::linspace(1,md_rp.n_rows,md_rp.n_rows)));
        arma::vec prob = counts/sum(counts);
        arma::vec p = prob(find(prob > 0));
        entropy = -arma::accu(p % arma::log(p));
    
        
        // relative entropy following form Moreno and Dale
        double nbins = p.n_elem;
        double denom = -1*std::log(1/nbins);
        rentropy = entropy/denom;
    }
    
    if(vert_lines.n_elem > 0){
        lam = (sum(vert_lines)/arma::accu(md_rp))*100;
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


arma::vec md_diagonal_lines(arma::mat& md_rp, int mindiagline){
    
    arma::vec lengths(1);
    lengths(0) = 0;
    
    // start interation with first diagonal
    // extract diagonals one by and and find line lengths
    int start = md_rp.n_rows;
    int stop = md_rp.n_rows;
    start *= -1;
    start += 1;
    
    for (int i = start; i < stop; ++i){
        
        //  check if vertical line contains any recurrent points
        // if not, skip to the next line and avoid those computations
        if (arma::accu(md_rp.diag(i)) > 0.0001){
            
            // zero pad to capture patterns that start and end with 1
            arma::vec pad1(10);
            arma::vec d = arma::join_cols(arma::join_cols(pad1, md_rp.diag(i)),pad1);
            
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


arma::vec md_vertical_lines(arma::mat& md_rp, int minvertline){
    
    arma::vec lengths(1);
    lengths(0) = 0;
    // start interation with first column
    // extract columns one by and and compute
    int n_cols = md_rp.n_cols;
    for (int i = 0; i < n_cols; ++i){
        
        //  check if vertical line contains any recurrent points
        // if not, skip to the next line and avoid those computations
        if (arma::accu(md_rp.col(i)) > 0.0001){
            // zero pad to capture patterns that start and end with 1
            // TODO: TRY TO ELIMINATE SOME UNWANTED STEPS HERE
            arma::vec pad1(10);
            // arma::vec pad2(10);
            // arma::vec d2 = md_rp.col(i);
            arma::vec d3 = arma::join_cols(pad1, md_rp.col(i));
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
x = cbind(x,x,x)
test = mdrqa(ts = x, embed = 1, delay = 1, normalize = 1, rescale = 1,
            mindiagline = 2, minvertline = 2, t_win = 1, radius = .00001,
            recpt = 0)

print(test$rqa)

x = sample(0:2, 300, replace = TRUE)
x = cbind(x,x,x)
test = mdrqa(ts = x, embed = 1, delay = 1, normalize = 1, rescale = 1,
           mindiagline = 2,minvertline = 2, t_win = 0, radius = .00001,
           recpt = 0)

print(test$rqa)
*/














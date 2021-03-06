

# ts is a N x p matrix
win_rqa = function(x, embed = 1, delay = 1, normalize = 1,
                  rescale = 1, mindiagline = 2, minvertline = 2,
                  t_win = 0, radius = 0.0001,  width, by){
    
    if(is.null(dim(x))){
        stop("\nx must be a matrix. Use win_rqa for vectors.")
    }
    if(ncol(x)>2){
        stop("\ncrqa can only be used on a bivariate series.\n
             win_mdrqa and win_jrqa may be more appropiate.")
    }
    
    N = length(x)
    start = 1
    stop = width
    win_result = NULL
    i = 1
    while (stop < N){
        tempx = x[start:stop,1]
        tempx = x[start:stop,2]
        win_result[[i]] = data.frame(rqa(ts1 = tempx, ts2 = tempx, 
                                embed = embed, delay = delay, 
                                normalize = normalize, rescale = rescale, 
                                mindiagline = mindiagline, 
                                minvertline = minvertline, t_win = t_win,
                                radius = radius, recpt = 0))
        
        i = i + 1
        start = start + by
        stop = stop + by
        
        
    }
    
    out = Reduce(rbind,win_result)
    names(out) = c('rr','det','div','nrline','ratio','maxline',
                          'meanline','lam','tt','vmax','entropy',
                          'rentropy')

    return(out)
    
}
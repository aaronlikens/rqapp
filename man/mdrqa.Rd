% FILE rqapp/man/rqa.Rd
\name{mdrqa}
\alias{mdrqa}
\title{
    Multi-Dimensional Recurrence Quantification Analysis
}

\description{
    This function performs multi-dimensional recurrence quantification anlaysis
    on bivariate and multivariate time series. 
}

\usage{
mdrqa(ts, embed, delay, normalize, rescale, mindiagline, minvertline, 
	t_win, radius, whiteline, recpt)	
}

\arguments{
	\item{ts}{a multivariate time series}
	\item{embed}{embedding dimension}
	\item{delay}{optimal time delay}
	\item{normalize}{should time series be normalized? (0 = no, 1 = unit 
		interval, 2 = z-score)}
	\item{rescale}{should distance matrix be rescaled? (0 = no, 1 = max 
		norm, 2 = min norm)}
	\item{mindiagline}{smallest number of diagonal points to be considered 
		a line}
	\item{minvertline}{smallest number of vertical points to be considered 
		a line}
	\item{t_win}{theiler window}
	\item{radius}{minimum distance within which points are considered 
		recurrent}
	\item{recpt}{should recurrence plot be returned? Not recommended for long
		 series.}

}

\details{
	This function performs multi-dimensional recurrence quantification anlaysis
	(MDRQA) on multivariate time series. Contrary to the appraoch used in 
	univariate series and JRQA, MDRQA uses direct measurement of state space
	variables in order to construct the recurrence plot without going through 
	the process of phase space reconstruction. 

	After constructing the recurrence plot for the multivariate series, a
	number of measures are computed to characterize recurrent stucture in the 
	time series.  These measures and their interpretation are well documented in
	the literature.  We provide simple definitions for each recurrence metric 
	below.  In addition, we provide references to standard readings including a 
	very readable introduction to RQA (i.e., Wallot et al., 2016; Marwan et 
	al., 2007).
}

\value{
	This function returns conventional JRQA variables and (optionally) the 
		joint recurrence plot. 
	\item{rr}{recurrence rate: the overall percentage of recurrent points}
	\item{det}{determinism: the percentage of recurrent points that fall on a
	    line}
	\item{div}{divergenc: inverse of determinism i.e. 1/det}
	\item{nrline}{number of lines: total number of lines in the upper triangle}
	\item{ratio}{ratio: percent determinism/percent recurrence i.e det/rr}
	\item{maxline}{longest line: the number points in the longest diagonal line}
	\item{meanline}{average line: average length of diagonal lines}
	\item{lam}{laminarity: perecentage of points that fall on vertical lines}
	\item{tt}{trapping time: average length of vertical lines}
	\item{vmax}{longest vertical line: the number of points in the longest 
	    vertical line}
	\item{entropy}{Shannon entropy: based on distribution of line lengths}
	\item{rentropy}{relative entropy: Shannon entropy normalized by number of
	    lines}
}

\references{
	Wallot, S., Roepstorff, A., & Mønster, D. (2016). Multidimensional 
    	Recurrence Quantification Analysis (MdRQA) for the analysis of 
	    multidimensional time-series: A software implementation in MATLAB and
	    its	application to group-level data in joint action. \emph{Frontiers in
	    psychology}, \emph{7}, 1835.
	
	Marwan, N., Romano, M. C. Theil, M., & Kurths, J. (2007). Recurrence plots
	    for the analysis of complex systems. \emph{Physics Reports}, \emph{438},
	    237-329.
}

\examples{
    \dontrun{
    	x <- sample(0:2, 100, replace = TRUE)
    	y <- sample(0:2, 100, replace = TRUE)
    	z <- sample(0:2, 100, replace = TRUE)
    	X = cbind(x,y,z)
    	x.recpt = mdrqa(ts = X, embed = 1, delay = 1, rescale = 0, normalize =
    	1, mindiagline = 2, minvertline = 2, t_win = 0, radius = 0.001, 
    	recpt = 1)
    	recurrence_plot(x.recpt$rp)
    
    }

}

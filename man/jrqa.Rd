% FILE rqapp/man/rqa.Rd
\name{jrqa}
\alias{jrqa}
\title{
    Joint Recurrence Quantification Analysis
}

\description{
    This function performs joint recurrence quantification anlaysis on
    bivariate and multivariate time series. 
}

\usage{
jrqa(ts, embed, delay, normalize, rescale, mindiagline, minvertline, 
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
	\item{whiteline}{not implemented}
	\item{recpt}{should recurrence plot be returned? Not recommended for long
		 series.}

}

\details{
	This function performs joint recurrence quantification anlaysis (JRQA) on 
	multivariate time series that have each been embedded (potentionally) in
	higher dimension than the originating series. Similar to the approach for
	univariate series, JRQA involves several steps:  First, for each time
	series, identify the optimal time \code{delay} as either the first zero
	crossing of the autocorrelation	function or the first minimum of the average
	mutual information function. Second, each time series is unfolded into
	\code{embed} dimensions by creating time-delayed copies of the original
	series. One method for determining the number of dimensions is by the
	method of False Nearest Neighbors \link{fnn}. Third, a distance matrix is
	computed among the embedded points of the series. A recurrence plot is
	constructed by vpassing the distance matrix through a heavyside function:
	distances less than or equal to the chosen \code{radius} are marked as 1
	(recurrent); distances falling outside the radius are marked as 0 (not 
	recurrent). 

	After constructing the recurrence plot for each series, the Hadamard product
	(pointwise product) is taken across all the resultant recurrence plots. Then
	a number of measures are computed to characterize recurrent stucture in the 
	time series.  These measures and their interpretation are well documented in 
	the literature.  We provide simple definitions for each recurrence metric 
	below.  In addition, we provide references to standard readings including a 
	very readable introduction to RQA (i.e., Webber & Zbilut, 2005; Marwan et 
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
	Webber, C. L., & Zbilut, J. P. (2005). Recurrence quantification analysis
		 of nonlinear dynamical time series. In S. Riley and G. C. Van Orden 
		 (eds). \emph{Tutorials in contemporary nonlinear methods for the
		 behavioral sciences}.
	
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
    	x.recpt <- jrqa(X, 1, 1, 0, 1, 2, 2, 0, .0001, 0, 1)
    	recurrence_plot(x.recpt$rp)
    
    }

}

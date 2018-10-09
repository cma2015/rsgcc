\name{gcc.corfinal}
\alias{gcc.corfinal}

\title{
	get the final correlaiton and p-value of Gini method
}
\description{
	Compare two correlations produced by GCC method for a gene pair, and choose one as the final output of GCC method.
}


\usage{
	gcc.corfinal(gcccor)
}


\arguments{
  \item{gcccor}{
	a list output by cor.pair function for GCC method.
}
}



\details{
	If the p-value is "NA", the correlation with absolute maximum value is selected; otherwise, the correlation with lower p-value is chosen.
}
\value{
  \item{gcc.fcor }{the final correlation of GCC.}
  \item{gcc.fpavlue }{the final pvalue of correlation.}
%% ...
}


\author{
Chuang Ma, Xiangfeng Wang
}


\seealso{
      \code{\link{onegcc}}, \code{\link{cor.pair}}.
}
\examples{

\dontrun{
   data(rsgcc)
   x <- rnaseq[1:4,]

   #compute correlation between 1th and 4th genes
   #significance level of the computed correlation 
   #is calcuated with 200 permutation tests.
   corpair <- cor.pair(c(1,4), GEMatrix = x, rowORcol = "row", 
                       cormethod = "GCC", pernum = 200, 
                       sigmethod = "two.sided") 

   #get the final correlation and p-value of GCC method 
   gcc.corfinal(corpair)  
}

}


\keyword{correlation}


\name{uniqueTissues}
\alias{uniqueTissues}

\title{
	get tissue information
}
\description{
	This function reads the sample names of genes and get unique tissue information for further tissue-specific genes finding and clustering.
}
\usage{
uniqueTissues(x)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	a numeric matrix containing gene expression value. The column labels are samples names. For two samples from the same tissue T, their names should be assigned as T.1 and T.2, respectively.
}
}


\value{
	A data matrix in which the elements is 0 (the sample not from the tissue) or 1 (the sample from the tissue)
}


\author{
Chuang Ma, Xiangfeng Wang
}



\seealso{
\code{\link{getsgene}}, \code{\link{gcc.tsheatmap}}.
}
\examples{

\dontrun{
   data(rsgcc)
   x <- rnaseq
   uniqueTissues(x)
}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tissue specific }

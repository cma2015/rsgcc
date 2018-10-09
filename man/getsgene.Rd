\name{getsgene}
\alias{getsgene}


\title{
	identify tissue(or condtion)-specific genes
}

\description{
	This function identifies tissue(or condition)- specific genes by considering the difference between the mean expression value of one tissue and the max expression value of other tissue.
}
\usage{
getsgene(x, Log = FALSE, Base = 2, AddOne = FALSE, 
         tsThreshold = 0.95, MeanOrMax = "Mean", Fraction = TRUE)
}


\arguments{
  \item{x}{
	a numeric matrix containing gene expression value. The column labels are samples names. For two samples from the same tissue T, their names should be assigned as T.1 and T.2, respectively.
}
  \item{Log}{
	logical indicating whether the gene expression value would be log-transformed.
}
  \item{Base}{
	a numeric value specifying the base of logarithm.
}
  \item{AddOne}{
	logical indicating if add one for avoding the problem of log-zero.
}
  \item{tsThreshold}{
	a numeric value giving the threshold of tissue specificity score. The tissue specificity score is 1, if the gene is only expressed in one tissue. Otherwise, the tissue specificity socre will be smaller than 1.
}
\item{MeanOrMax}{
	character "Mean" or "Max" indicate the mean or maximal expression value will be calculated for the tissue of interest. 
}
  \item{Fraction}{
	logical indicating whether the gene expression value would be scaled across tissues. 
}
}
\details{
	The tissue specificity score is calculated with the formula 1-min(R(1), R(2), ..., R(i),..., R(n)), where R(i) = M(i)/E(i), E(i) is the mean or maximal expression value of tissue i, and M(i) is the maximal expression value of other tissues. If the tissue specificity score higher than tsThreshold, then the gene is considered as tissue specifically expressed. 
	
	If Fraction is TRUE, the expression values of a gene is scaled accorss the tissues with the formula e(i)/(e(1)+e(2)+...+e(n)). e(i) is the expression value of the consider gene in ith sample.
}
\value{
	A list with following components:
    \item{csGenes }{a data matrix containing expression vlaues of tissue specific genes.}
    \item{csScoreMat }{a data matrix with three columns containg the gene index information from x, tissue specificity score and the tissue information with the tissue specificity score.}
}



\author{
Chuang Ma, Xiangfeng Wang.
}



\references{
[1] Chuang Ma, Xiangfeng Wang. Machine learning-based differential network analysis of transcriptomic data: a case study of stress-responsive gene expression in Arabidopsis thaliana. 2013 (Submitted).
}



\examples{

\dontrun{
   data(rsgcc)
   tsRes <- getsgene(rnaseq, tsThreshold = 0.75, MeanOrMax = "Mean", Fraction = TRUE)
}

}


\keyword{ tissue specific }

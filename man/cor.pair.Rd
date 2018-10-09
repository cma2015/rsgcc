\name{cor.pair}
\alias{cor.pair}
\alias{gcc.cor.pair}
\alias{gcc.corpair}

\title{
	compute the correlation between two genes
}



\description{
	This function can compute the correlation of a pair of genes with Gini correlation and four other correlation methods. The signficance level (p-value) of computed correlation can be estimated with the permutation test method.
}


\usage{
cor.pair(idxvec, 
      GEMatrix,
      rowORcol = c("row", "col"),
      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
      pernum = 0,
      sigmethod = c("two.sided", "one.sided"))
}



\arguments{
  \item{idxvec}{
	a numer vector containing two elements for the index of genes or samples in GEMatrix (e.g., c(1,2) ).
  }
  \item{GEMatrix}{
	a data matrix containing numeric variables. Example: rows correspond to genes and columns to samples. This parameter is the same as the "GEMatrix" defined for cor.matrix.
  }
  \item{rowORcol}{
	a character string ("row" or "col") indicating gene expression data will be extracted by rows or columns for correlation calculation. "row": correlation between two genes. "col": correlaiton between two samples.
  }
  \item{cormethod}{
	a character string that specifies the correlation method to be used for correlation calculation. 
  }
  \item{pernum}{
	the number of permutation test used for calcluating statistical significance level (i.e., p-value) of correlations. 
  }
  \item{sigmethod}{
	a character string ("two-sided" or "one-sided") specifying the method used to compute p-value for permutation test.
  }
}



\value{
	A list with the following components:
  \item{gcc.rankx }{a Gini correlation produced by using the rank information of the first gene (i.e., the first element in idxvec).}
  \item{gcc.ranky }{a Gini correlation produced by using the rank information of the second gene (i.e., the second element in idxvec).}
  \item{gcc.rankx.pvalue }{p-value of gcc.rankx.}
  \item{gcc.ranky.pvalue }{p-value of gcc.ranky.}
  \item{cor }{the correlation produced by "PCC", "SCC", "KCC" or "BiWt".}
  \item{pvalue }{the p-value of cor.}

}


\author{
Chuang Ma, Xiangfeng Wang
}


\note{
(1) To perform BiWt, the R package "biwt" should be installed in advance. 

(2) When the cormethod is defined as "GCC", this function will output a list with four numeric elements: gcc.rankx, gcc.ranky, gcc.rankx.pvalue, gcc.ranky.pvalue. Otherwise, it will output a list with two elements (cor and p-value)
}



\seealso{
   \code{\link{onegcc}}, \code{\link{cor.matrix}}, \code{\link{gcc.corfinal}}.
}




\examples{
\dontrun{
    data(rsgcc)        #load the sample data in rsgcc package
    x <- rnaseq[1:4,]  #construct a GEMatrix with the RNA-Seq data of the first four genes

   #compute correlation between the 1st and 4th genes
   corpair <- cor.pair(c(1,4), GEMatrix = x, rowORcol = "row", 
                       cormethod = "GCC", pernum = 0, 
                       sigmethod = "two.sided")  
   #compute correlation between the 1st and 4th genes, 
   #the p-value of correlation will be estimated with 2,000 permuation test.
   corpair <- cor.pair(c(1,4), GEMatrix = x, rowORcol = "row", 
                       cormethod = "GCC", pernum = 2000, 
                       sigmethod = "two.sided")  

   #compute correlation between the 1st and 4th samples
   corpair <- cor.pair(c(1,4), GEMatrix = x, rowORcol = "col", 
                       cormethod = "GCC", pernum = 0, 
                       sigmethod = "two.sided")  
}
}

\keyword{ correlation }


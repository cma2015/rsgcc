\name{adjacencymatrix}
\alias{adjacencymatrix}


\title{
adjacency matrix calculation
}
\description{
This function generates the adjacency matrix for network re-construction from gene expression data with different methods including Gini correlation (GCC), Pearson correlation (PCC), Spearman correlation (SCC), Kendall correlation (KCC),Tukey's biweight correlation coefficient (BiWt), mutual information (MI), and maximal information-based nonparametric exploration (MINE) statistic methods. Euclidean distance (ED) between two genes can also be calculated. It was implemented these methods in C language and parallel mode, and thus is greatly faster than the cor.matrix function.
}
\usage{
adjacencymatrix(mat, genes.row = NULL, genes.col = NULL, 
     method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"), 
     k = 3, cpus = 1, 
     saveType = "matrix", 
     backingpath = NULL, 
     backingfile = "adj_mat", 
     descriptorfile = "adj_desc", 
     ... )
}


\arguments{
  \item{mat}{a data matrix containing gene expression dataset where rows defines for genes and columns for samples.}
  \item{genes.row}{if genes.row and genes.col are not NULL, a subset of genes will be selected for correlation calcuation and set as the rownames of adjacency matrix. Currently, doesn't work for BiWt and MINE}
  \item{genes.col}{if genes.row and genes.col are not NULL, a subset of genes will be selected for correlation calcuation and set as the colnames of adjacency matrix.Currently, doesn't work for BiWt and MINE}
  \item{method}{a method used to calculate the association between a pair of genes.}
  \item{k}{the number of nearest neighbors to be considered for estimating the mutual information. Must be less than the number of columns of mat, and only work for the mutual information(MI) method.}
  \item{cpus}{the number of cpus will be used for calcuation.}
  \item{saveType}{the type (matrix or bigmatrix) specified for the output.}
  \item{backingpath}{the path used to save big matrix. If it is NULL, current working directory will be used. Works only when the saveType is "bigmatrix".}
  \item{backingfile}{the file name of big matrix. Works only when the saveType is "bigmatrix".}
  \item{descriptorfile}{the description file of big matrix. Works only when the saveType is "bigmatrix".}
  \item{\dots}{Further parameters passed for MINE method. More information can be found in R package minerva.}
}


\value{
  \item{value}{a matrix (or big.matrix) recording the associations between the gene pairs.}
}


\references{
[1] Ma C and Xiang XF. Application of the Gini correlation coefficient to infer regulatory relationships in transcriptome analysis, Plant Physiology, 2012, 160(1):192-203.

[2] Sales G and Romualdi C. Parmigene-a parallel R package for mutual information estimation and gene network reconstruction. Bioinformatics, 2012, 27:1876-1877.

[3] Davide Albanese, Michele Filosi, Roberto Visintainer, et al. minerva and minepy: a C engine for the MINE stuite and its R, Python and MATLAB wrappers. Bioinformatics, 2013, 29(3): 407-408.

[4] David N. Reshef, Yakir A. Reshef, Hilary K. Finucane, et al. Detecting novel associations in large data sets. Science, 2011, 334(6062): 1518-1524.

[5] Johanna Hardin, Aya Mitani, Leanne Hicks and Brian VanKoten. A robust measure of correlation between two genes on a microarray. BMC Bioinformatics, 2007, 8:220.

}


\note{
1) The mutural information estimation is based on k-nearest neighbor distance (Sales G and Romualdi C, 2012). Thus the parameter k only works for the mutual information method.

2) Two correlations can be produced by the GCC method by reciprocally using the rank and value information of one gene (or variable). Here the correlation with the maximum absolute values is selected for generating the adjacency matrix. 

3) More information about the big.matrix can be found in bigmemory package. 

}



\examples{
\dontrun{
 mat = matrix(rnorm(180), nrow = 10)
 rownames(mat) <- c(1:10)
 colnames(mat) <- c(1:18)
 mat

 #using GCC method to compute the correlation of all gene paris
 adjacencymatrix( mat, method = "GCC", cpus = 2 )

 #using GCC method to compute the correlation of a subset of gene pairs
 adjacencymatrix( mat = mat, genes.row = c(1:5), genes.col = c(5:8), method = "GCC", cpus = 2 )

 

 #for MI method, k works here.
 adjacencymatrix( mat, method = "MI", k= 3)
}

}


\keyword{ correlation }


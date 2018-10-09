\name{cor.matrix}
\alias{cor.matrix}
\alias{gcc.cor.matrix}
\alias{gcc.cormatrix}

\title{
	correlation calculation for a set of genes
}

\description{
	This function provides five correlation methods (GCC, PCC, SCC, KCC and BiWt) to calculates the correlations between a set of genes.
}


\usage{cor.matrix(GEMatrix,
      cpus = 1,
      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
      style = c("all.pairs", "pairs.between", "adjacent.pairs", "one.pair"),
      var1.id = NA,
      var2.id = NA,
      pernum = 0,
      sigmethod = c("two.sided", "one.sided"),
      output = c("matrix", "paired"))
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GEMatrix}{
	a data matrix containing the gene expression data of a set of genes. Each row of the GEMatrix corresponds to a gene, and each column corresponds to the expression level in a sample. 
}
  \item{cpus}{
	the number of cpus used for correlation calcluation.
}
  \item{cormethod}{
	a character string that specifies a correlation method to be used for correlation calculation.
}
  \item{style}{
	a character string that indicates the all or partial genes to be used for correlation calculation.
}
  \item{var1.id}{
	a numeric vector specifying the row numbers of genes.
}
  \item{var2.id}{
	a numeric vector specifying the row numbers of genes. Suppose the var1.id and var2.id are respectively c(1,2) and c(3,6), then the the correlation of gene pairs (G1,G3) and (G2,G6) will be calcuated. For styles of "pairs.between" and "one.pair", this parameter MUST be pre-defined. For the other styles, this parameter can be automatically defined by the program itself. 
}
  \item{pernum}{
	the number of permutation test used for calcluating statistical significance level (i.e., p-value) of correlations. 
}
  \item{sigmethod}{
	a character string ("two-sided" or "one-sided") that specifies the method used to compute p-value for permutation test. 
}
  \item{output}{
	a character string ("matrix" or "paired") that represents the output format of correlations. Specifiying the "matrix" will output two matrix for correlations and p-values, respectively. Specifiying the "paired" will output only one matrix, in which each row provides the information of gene pair, the correlation and p-value.
}

}



\details{
	Given a data matrix (e.g., microarray and RNA-Seq gene expression matrix), calculating correlation with GCC and other correlation methods for partial(or all) individuals (e.g., genes). The statistical significance (i.e., p-value) of each correlation is derived from the permutation test. Parallel computing options are also provided for speeding up the correlation calculation.  	

}



\value{
	A list with the following components:
  \item{corMatrix }{correlation of gene pairs shown in matrix form. This data matrix is generated only when the output format "matrix" is specified. }
  \item{pvalueMatrix }{p-value of correlations shown in matrix form. This data matrix is generated only when the output format "matrix" is specified. }
  \item{corpvalueMatrix }{correlation and p-values listed in one form. This data matrix is generated only when the output format "paired" is indicated. }
}




\author{
Chuang Ma, Xiangfeng Wang
}


\note{

(1) The rsgcc provides the RNA-Seq profiled expression level of 100 genes as a sample data to implement cor.matrix, cor.pari and other functions in the package. After running the command: data(rsgcc), the expression data of these genes will be loaded to the GEMatrix "rnaseq". The user can also load the GEMatrix from the gene expression file, which should be in a textual format of a gene expression matrix. An example of the gene expression file(e.g., "/home/rsgcc/geneExpFile.txt") is shown as follow: 

sample1	sample2	sample3	sample4

gene1	45	65	77	75

gene2	75	78	83	39

gene3	2	11	10	6


Then the GEMatrix can be obtained by load this gene expression file with the command: 
x <- as.matrix(read.table("/home/rsgcc/geneExpFile.txt"))

(1) var1.id and var2.id should be defined with the numeric vector format for "pairs.between", or "one.pair" styles. 

(2) To perform BiWt, the R package "biwt" should be installed in advance. 

(3) To perform the parallel compution, the "snowfall" package in R should be installed in advance.
}


\seealso{
	\code{\link{cor.pair}}, \code{\link{onegcc}}, \code{\link{cor.test}}.
}



\examples{

\dontrun{
   data(rsgcc)        #load the sample data in rsgcc package
   x <- rnaseq[1:4,]  #construct a GEMatrix with the RNA-Seq data of the first four genes

   #run on one CPU for all the possible gene pairs in the GEMatrix "x".
   #do not cacluate the p-value of computated correlations.
   cor.matrix(x, cpus = 1, 
              cormethod = "GCC", style = "all.pairs", 
              pernum = 0, sigmethod = "two.sided", 
              output = "matrix")  

   #run on two CPUs, snowfall package should be properly installed.
   #cacluate the p-value of correlations with the 2000 permutation tests.
   #output the results in "paired" format.
   cor.matrix(x, cpus = 2, 
              cormethod = "GCC", style = "all.pairs", 
              pernum = 2000, sigmethod = "two.sided", 
              output = "paired")  

   #calcuate correlation on the pairs between the 1st, 2nd and 3rd genes in the GEMatrix "x".
   cor.matrix(x, cpus = 1, 
              cormethod = "GCC", style = "pairs.between", 
              var1.id = c(1:3), var2.id = c(1:3),
              pernum = 2000, sigmethod = "two.sided", 
              output = "matrix")

  #calcuate correlation on the adjacent genes ((G1,G2), (G2,G3), (G3,G4),...) in the GEMatrix "x".
   cor.matrix(x, cpus = 1, 
              cormethod = "GCC", style = "adjacent.pairs", 
              pernum = 2000, sigmethod = "two.sided", 
              output = "matrix")

}

}

\keyword{ correlation }


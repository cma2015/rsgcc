\name{gcc.hclust}
\alias{gcc.hclust}


\title{
	hierarchical cluster
}


\description{
	Hierarchical cluster analysis of microarrany and RNA-Seq gene expression data with Gini correlation and four other correlation methods.
}


\usage{
gcc.hclust(x, 
      cpus = 1, 
      method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
      distancemethod = c("Raw", "Abs", "Sqr"),
      clustermethod = c("complete", "average", "median", 
                        "centroid", "mcquitty", "single", "ward"))
}


\arguments{
  \item{x}{a data matrix containing numeric variables, which is the same as the GEMatrix defined in the cor.matrix function.}
  \item{cpus}{the number of cpus used for computation.}
  \item{method}{a character string indicating the method to be used to calculate the associations.}
  \item{distancemethod}{a character string specifying the distance method to be used. Currently, three distance methods are available, include: "Raw" (1-cor)", "Abs" (1-|cor|), and "Sqr" (1-|cor|^2).}
  \item{clustermethod}{the distance measure to be used. This must be one of "complete", "average", "median", "centroid", "mcquitty", "single", or "ward".}
}

\details{
	This function generate the cluster tree with different distance measures for clustering analysis of microarray and RNA-Seq gene expression data by integrating the hclust function of stats package in R (http://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html). Similar to the hclust, the values output by gccdist can be directly used to plot cluster trees with plot function.  
}


\value{
	A list with the following components:
  \item{hc}{an object describes the tree information produced by the clustering process. This object is also a list with five components: "merge" is a numeric matrix with n-1 rows and 2 columns. n is the number of used individuals (e.g., genes). Row i describes the merging of clusters at step i of the clustering. "order" is a vector giving the order of individuals for tree cluster plotting. "height" is a vector with n-1 numeric values associated with the distance measure for the particular cluster method. "labels" are labels of the individuals being clustered. "method" is the distance measure used for cluster analysis. See details for the description in hclust function of stats package.}
  \item{dist}{a data matrix containing the distances between different genes.}
  \item{pairmatrix }{a data matrix including the correlation between different genes.}
}


\author{
	Chuang Ma, Xiangfeng Wang
}


\examples{

\dontrun{

   #obtain gene expression data of 10 genes.
   data(rsgcc)
   x <- rnaseq[1:10,]
      
   #hierarchical clustering analysis of these 10 genes with GCC method
   hc <- gcc.hclust(x, cpu=1, method = "GCC", 
                    distancemethod = "Raw", clustermethod = "complete")

   #plot cluster tree
   plot(hc$hc)
}

}


\keyword{ cluster }

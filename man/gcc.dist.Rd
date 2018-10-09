\name{gcc.dist}
\alias{gcc.dist}

\title{
	compute distance matrix for hierarchical clustering
}

\description{
	This function computes the distance between the rows of a data matrix with the specified distance method. 
}


\usage{
gcc.dist(x, 
      cpus = 1,
      method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
      distancemethod = c("Raw", "Abs", "Sqr"))
}

\arguments{
  \item{x}{a data matrix containing numeric variables, which is the same as the "GEMatrix" defined in the cor.matrix function.}
  \item{cpus}{the number of cpus used for computation.}
  \item{method}{a character string indicating the method to be used to calculate the associations. }
  \item{distancemethod}{a character string indicating the distance method to be used. Currently, three distance methods are available, include: "Raw" (1-cor)", "Abs" (1-|cor|), and "Sqr" (1-|cor|^2). }
  
}



\value{
	A list with the following components:
  \item{dist }{a data matrix containing the distances between different genes.}
  \item{pairmatrix }{a data matrix including the correlation between different genes.}
%% ...
}


\author{
Chuang Ma, Xiangfeng Wang
}


\seealso{
 \code{\link{cor.matrix}}, \code{\link{gcc.hclust}}, \code{\link{gcc.tsheatmap}}.
}



\keyword{ cluster }

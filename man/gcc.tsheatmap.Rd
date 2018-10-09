\name{gcc.tsheatmap}
\alias{gcc.tsheatmap}


\title{
	correlaiton and clustering analysis of tissue-specific genes 
}


\description{
	This function performs the correlaiton and clustering analysis of tissue-specific genes with expression data generated from microarray and RNA-Seq experiments.
}


\usage{
gcc.tsheatmap(x,
                           
      cpus = 1,
                        
      ## correlation method
      method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"),
                          
      distancemethod = c("Raw", "Abs", "Sqr"),
                           
      #cluster method
      clustermethod = c("complete", "average", "median", 
                        "centroid", "mcquitty", "single", "ward"),
                         
      #hcdata by output gcc.tsheatmap 
      rowhcdata = NULL,
      colhcdata = NULL,    
                          
                          
      keynote = "FPKM",
                        
      ## dendrogram control
      symm = FALSE,
                      
      ## data scaling
      scale = c("none","row", "column"),
      na.rm=TRUE,

      ## image plot
      revC = identical(Colv, "Rowv"),
      add.expr,

      ## mapping data to colors
      breaks,
      symbreaks=min(x < 0, na.rm=TRUE) || scale!="none",

      ## colors
      colrange = c("yellow", "red"), 
                          
      tissuecol= "heat.colors",

      ## block sepration
      colsep = 0.15,
      rowsep,
      sepcolor="white",
      sepwidth=c(0.05,0.05),
          
      ## level trace
      trace=c("none","column","row","both"),
      tracecol="cyan",
      hline=median(breaks),
      vline=median(breaks),
      linecol=tracecol,

      ## plot margins
      margins = c(5, 5),

      ## plot labels
      main = NULL,
      xlab = NULL,
      ylab = NULL,

      ## plot layout
      lmat = NULL,
      lhei = NULL,
      lwid = NULL,

      ## extras
      ...)
}


\arguments{
  \item{x}{
	a data matrix containing numeric variables. Example: rows may correspond to genes and columns to samples. 
}
  \item{cpus}{
	the number of cpus used for correlaiton calcluation. snowfall package in R needed to be installed in advance.
}
  \item{method}{
	a character string that specifies a correlation method to be used for association calculation.
}
  \item{distancemethod}{
	a character string specifying the distance method to be used. Currently, three distance methods are available, include: "Raw" (1-cor)", "Abs" (1-|cor|), and "Sqr" (1-|cor|^2).
}
  \item{clustermethod}{
	the distance measure to be used. This must be one of "complete", "average", "median", "centroid", "mcquitty", "single", or "ward".
}
  \item{rowhcdata}{
	the object of class hc generated from gcc.hclust for rows in x.
}

  \item{colhcdata}{
	the object of class hc generated from gcc.hclust for columns in x.
}


  \item{keynote}{
	a character string indicating the lable of color key.
}
  \item{symm}{
	logical indicating if x should be treated as a symmetrical matrix.
}
  \item{scale}{
	a character string specifying if the data values would be centered and scaled by rows or by columns, or none.
}
  \item{na.rm}{
	logical indicating whether the Nas should be eliminated.
}
  \item{revC}{
	logical indicating if the column order should be reversed for plotting.
}
  \item{add.expr}{
	expression that will be evaluated after the call to image.
}
  \item{breaks}{
	(optional)Either a integer number specifying the break points to be used, or a numeric vector indicating the splitting points for binning x into colors. 
}
  \item{symbreaks}{
	Boolean indicating whether breaks should be made symmetric about 0. This option works if the quanbreaks is FALSE.
}
  \item{colrange}{
	colors used for the image. It could be a function(i.e., heat.colors) or a vector of colors with at leat two elements (e.g., c("green", "black", "red")). 
}
  \item{tissuecol}{
	colors for tissues. tissuecol could be a function(i.e., heat.colors) or a vector of colors for different tissues.
}
  \item{colsep}{
	(optional) vectors of integers indicating which columns should be seperated from the preceding columans by a narrow space of color sepcor.
}
  \item{rowsep}{
	(optional) vectors of integers indicating which rows should be seperated from the preceding rows by a narrow space of color sepcor.
}
  \item{sepcolor}{
	(optional) color used to seperate rows or columns.
}
  \item{sepwidth}{
	(optional) A numeric vector containing two elements giving the width (colsep) or height (rowsep) for the seperation of columns or rows.
}
  \item{trace}{
	character string indicating a solid "trace" lined should be drawn across "rows", or "column", or "both" or "none".
}
  \item{tracecol}{
	color for trace
}
  \item{hline}{
	vector of values whithin cells where horizontal lines should be drawn with line col.
}
  \item{vline}{
	vector of values whithin cells where vertical lines should be drawn with line col.
}
  \item{linecol}{
	color for hline and vline.
}
  \item{margins}{
	a numeric vector containing 2 elements specifying the margins for column and row names, respectively. See (par(mar=*)).
}
  \item{main}{
	main title. defaults to none.
}
  \item{xlab}{
	x-axis label. defaults to none.
}
  \item{ylab}{
	y-axis label. defaults to none.
}
  \item{lmat}{
	position matrix for visual layout. 
}
  \item{lhei}{
	column height for visual layout.
}
  \item{lwid}{
	column width for visual layout. For instance, lwid = c(0.5, 0.05, 0.01, 0.5, 0.01, 0.05, 0.5)
}
  \item{\dots}{
	additonal arguments passed on to image.
}
}


\value{
	A list with the following components:
  \item{retval }{a list with components of "rowInd" (row index of heat map from x), "colInd" (column index of heat map from x), "call" (the match call), "carpet" (reordered and scaled 'x' values used generate the main 'carpet'), "rowDendrogram" (row dendrogram), "colDendrogram" (column dendrogram), "breaks" (break points for binning x), "col" (colors used), and "colorTable" (a data frame providing the lower and upper bound and color for each bin).}
  \item{hcr }{the values returned from gcc.hclust function for clustering individuals (e.g., genes) in row direction}
  \item{hcc }{the values returned from gcc.hclust function for clustering individuals (e.g., genes) in column direction}
}


\author{
Chuang Ma, Xiangfeng Wang
}



\seealso{
\code{\link{gcc.dist}}, \code{\link{cor.matrix}}, \code{\link{gcc.hclust}}, \code{\link{gcc.tsheatmap}}.
}


\examples{

\dontrun{
   data(rsgcc)
         
   #get expression matrix of tissue-specific genes
   tsRes <- getsgene(rnaseq, tsThreshold = 0.75, MeanOrMax = "Max", Fraction = TRUE)

   #heat map of tissue-specific genes
   thm <- gcc.tsheatmap(tsRes$tsgene, cpus = 1, method = "GCC", 
                    distancemethod = "Raw", clustermethod = "complete")
}
}




\keyword{cluster}
\keyword{tissue specific}

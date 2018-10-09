\name{gcc.heatmap}
\alias{gcc.heatmap}


\title{
	heat map 
}
\description{
	The heat map is a color imange representing the data in the a matrix. The dendrogram information are usually added to the left side and/or to the top for displaying the clustering information.
}
\usage{
gcc.heatmap(x, 
      cpus = 1, 
      method = c("GCC", "PCC", "SCC", "KCC", "BiWt", "MI", "MINE", "ED"), 
      distancemethod = c("Raw", "Abs", "Sqr"), 
      clustermethod = c("complete", "average", "median", 
                        "centroid", "mcquitty", "single", "ward"), 
  
      #hcdata output by gcc.tsheatmap function
      rowhcdata = NULL, 
      colhcdata = NULL, 

      keynote = "FPKM", 

      ## dendrogram control
      symm = FALSE, 
      Rowv = TRUE, 
      Colv = if (symm) "Rowv" else TRUE, 
      dendrogram = c("both", "row", "column", "none"), 

      ## data scaling
      scale = c("none", "row", "column"), 
      na.rm = TRUE, 
      revC = identical(Colv, "Rowv"), 
      add.expr, 

      #break points for binning values in x
      breaks = 16, 
      quanbreaks = TRUE, 
      symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",

      #colors 
      colrange = c("green", "black", "red"), 

      colsep, 
      rowsep,
      sepcolor = "white", 
      sepwidth = c(0.05, 0.05), 
      cellnote, 
      notecex = 1, 
      notecol = "cyan", 
      na.color = par("bg"), 
      trace = c("none", "column", "row", "both"), 
      tracecol = "cyan", 
      hline = median(breaks), 
      vline = median(breaks), 
      linecol = tracecol, 
      margins = c(5, 5), 
      ColSideColors, 
      RowSideColors, 
      cexRow = 0.2 + 1/log10(dim(x)[1]), 
      cexCol = 0.2 + 1/log10(dim(x)[2]), 
      labRow = NULL, 
      labCol = NULL, 

      #color key
      key = TRUE, 
      keysize = 0.65, 
      density.info = c("none", "histogram", "density"), 
      denscol = tracecol, 
      symkey = min(x < 0, na.rm = TRUE) || symbreaks, 
      densadj = 0.25, 

      #image information
      main = NULL, 
      xlab = NULL, 
      ylab = NULL, 
      lmat = NULL, 
      lhei = NULL, 
      lwid = NULL, 
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
  \item{Rowv}{
	logical determining if the row dendrogram should be reordered.
}
  \item{Colv}{
	logical determining if the columns dendrogram should be reordered.
}
  \item{dendrogram}{
	a character string indicating whether to draw the "none", "row", "column", "both" dendrograms.
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
	expression that will be evaluated after the call to image. Can be used to add components to the plot.
}
  \item{breaks}{
	(optional)Either a integer number specifying the break points to be used, or a numeric vector indicating the splitting points for binning x into colors. 
}
  \item{quanbreaks}{
	logical indicating if the splitting points for binning x in quantile scale. For instance, if quanbreaks is TRUE, the breaks would be quantile( unique(c(x)), probs = seq(0, 1, length = breaks), na.rm = TRUE). 
}
  \item{symbreaks}{
	Boolean indicating whether breaks should be made symmetric about 0. This option works when the quanbreaks is FALSE.
}
  \item{colrange}{
	colors used for the image. It could be a function(i.e., heat.colors) or a vector of colors with at leat two elements (e.g., c("green", "black", "red")). 
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
  \item{cellnote}{
	(optional)a matrix of character strings which will be placed within each color cell.
}
  \item{notecex}{
	(optional)numeric scaling factor for cellnot itmes.
}
  \item{notecol}{
	(optional)character string specifying the color of cellnote text. Default to "green".
}
  \item{na.color}{
	color to be used for missing value (NA). Defaults to the plot background color.
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
  \item{ColSideColors}{
	(optional)character string of colors for annotating the columns of heat map.
}
  \item{RowSideColors}{
	(optional)character string of colors for annotating the rows of heat map.
}
  \item{cexRow}{
	cex.axis for the row lables.
}
  \item{cexCol}{
	cex.axis for the column lables.
}
  \item{labRow}{
	character strings indicating the lables of rows. Default to rownames(x)
}
  \item{labCol}{
	character strings indicating the lables of columns. Default to colnames(x)
}
  \item{key}{
	logical indicating whether the color key would be draw.
}
  \item{keysize}{
	numeric value specifying the size of color key.
}
  \item{density.info}{
	character string indicating whether to superimpose a "histogram", a "density" plot, or not plot("none") on the color-key.
}
  \item{denscol}{
	character string giving the color for the density display specified by density.info, defaults to the same value as tracecol.
}
  \item{symkey}{
	Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values and to FALSE otherwise.
}
  \item{densadj}{
	Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key. Default to 0.25.
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
	position matrix for visual layout. See details from the help page of heatmap.2.
}
  \item{lhei}{
	column height for visual layout. See details from the help page of heatmap.2
}
  \item{lwid}{
	column width for visual layout. See details from the help page of heatmap.2
}
  \item{\dots}{
	additonal arguments passed on to image.
}
}


\details{
	This function plots the heat map of microarray and RNA-Seq gene expression data by modifying the scripts of heatmap.2 in R. The main modifications include: (1) designing several distance measures derived from Gini correlation and other correlation methods; (2) providing the option of quanbreaks for RNA-Seq data.
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

\note{
	This function clusters microarray and RNA-Seq gene expression data and plot heatmap by refining heatmap.2 function in gplots package. Therefore, most parameters and output values are defined similarly as those in heatmap.2.
}


\seealso{
\code{\link{gcc.dist}}, \code{\link{cor.matrix}}, \code{\link{gcc.hclust}}, \code{\link{gcc.tsheatmap}}.
}
\examples{

\dontrun{
   data(rsgcc)
   x <- rnaseq[1:50,]
   ghm <- gcc.heatmap(x, cpus = 1, method = "GCC", 
               distancemethod = "Raw", clustermethod = "complete", labRow = "")
}
  
}


\keyword{cluster}

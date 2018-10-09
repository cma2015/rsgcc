\name{rsgcc.gui}
\alias{rsgcc.gui}


\title{
	graphical user interface (GUI) of rsgcc package
}


\description{
	This function provides a graphical user interface (GUI) to perform the correlation and clustering analysis via a series of mouse actions without command-line based R programmming. The output of clustering informaiton in "CDT" format can be further visulized and analyzed by TreeView program.  
}


\usage{
      rsgcc.gui(margins = c(1, 1), labRow = "", labCol = "", 
                lwid = c(0.5, 0.05, 0.01, 0.5, 0.01, 0.05, 0.5), 
                keynote = "FPKM")
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{margins}{
	a numeric vector containing 2 elements specifying the margins for heat map. See (par(mar=*)).
}
  \item{labRow}{
	character strings indicating the lables of rows. Default to rownames(x).
}
  \item{labCol}{
	character strings indicating the lables of columns. Default to colnames(x).
}
  \item{lwid}{
	column width for visual layout.   
}
  \item{keynote}{
	a character string indicating the lable of color key.
}
}
\details{
	For heat map of ts-genes, rsgcc will run the gcc.tsheatmap function. the lwid could be c(0.5, 0.05, 0.01, 0.5, 0.01, 0.05, 0.5). The 2nd, 4th and 6th elements are column widths of color tissue bar, heat map and color key bar. The 3rd and 5th are the widths of separation for these three figures. The The first and last elements are the widths of "blank region" for displaying the labels of tissue and key bar.
	If the option "Find and cluster ts-genes" is not selected, rsgcc will call the gcc.heatmap for clustering analysis. In this case, the lwid should be a numeric vector with two elements(e.g.,c(0.65,4)). Here 0.65 is the width of column for row dendrogram. 4 is the width of heat map.

}


\author{
Chuang Ma, Xiangfeng Wang
}
\note{
	rsgcc.gui is built upon gWidgets package. Make sure "gWidgets" and "gWidgetsRGtk2" package is properly installed.

	The following is a guide of using the rsgcc GUI.

Step 1: Select a gene expression data file and load expression data with read.table or read.csv, which is decided by the program itself according to the file suffix (".txt", ".csv", or nothing). Each row of the table is one gene, and each column is the expression data of one sample. The column names are sample IDs indiciating the tissue informaiton (i.e., "T1.1", "T1.2", "T2"). After the data is alreadly loaded, you can click the option ("Display loaded data") to display the gene expression datain a new window. If the tissue-specific genes are interested for the clustering analysis, please select the option ("Find and cluster ts-genes") and specify the threshold of tissue-specificity score (Default 0.95. The threshold should be smaller than 1.0; the <ENTER> MUST be pressed to confirm the change ).

Step 2: Selecte a correlation method. Default: Gini correlation. 

Step 3: Specify a distance measure. Default: raw correlation (1-coef).

Step 4: Choose a cluster method. In current version, rsgcc includes seven cluster methods. More information about these cluster methods can be found at the help page of hclust function.

Step 5: Set a integer for the number of CPUs to be used. The snowfall package in R is needed for the paralleled computing to speed up the calculation of correlation coefficients. After you change the number, the <ENTER> MUST be pressed to confirm the change.

Step 6: Press the button "Start to run" to perform the correlation and clustering analysis of gene expression data. A heat map will be visualized at the right region of interface if the task is finished.

Step 7: Three bars can be slided to adjust colors in heat map. 

Step 8: After the appropriate colors have already been determined, you can "save correlation and cluster data". The correlations will be output to a file with three columns (gene1, gene2, correlation). The cluster information will be output into three files (the suffix are ".atr", ".gtr" and ".cdt") for visualizing and analyzing by Treeview program. A pdf file will also be generated for heat map. All these files in the same directory of gene expression data file. 
}


\seealso{
\code{\link{gcc.tsheatmap}}, \code{\link{gcc.heatmap}}.
}


\examples{

\dontrun{
   library("gWidgetsRGtk2")

   library(rsgcc)

   ## the GUI of rsgcc will show up after the GUI toolkit "gWidgetsRGtk2" is selected.
   rsgcc.gui()
}


}


\keyword{correlation}
\keyword{cluster}

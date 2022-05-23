# ___**rsgcc**: Gini methodology-based correlation and clustering analysis of microarray and RNA-Seq gene expression data___ </br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
![](https://tctechcrunch2011.files.wordpress.com/2014/06/apple_topic.png?w=220) </br>
**Brief introduction:**
This package provides functions for calculating associations between two genes with five correlation methods(e.g., the Gini correlation coefficient [GCC], the Pearson's product moment correlation coefficient [PCC], the Kendall tau rank correlation coefficient [KCC], the Spearman's rank correlation coefficient [SCC] and the Tukey's biweight correlation coefficient [BiWt], and three non-correlation methods (e.g., mutual information [MI] and the maximal information-based nonparametric exploration [MINE], and the euclidean distance [ED]). It can also been implemented to perform the correlation and clustering analysis of transcriptomic data profiled by microarray and RNA-Seq technologies. Additionally, this package can be further applied to construct gene co-expression networks (GCNs).</br>

## Installation ##
- **Note:** This is a variant of original R package rsgcc(https://cran.r-project.org/web/packages/rsgcc/index.html) which does not provide GUI interface because the installation of dependency package "gWidgetsRGtk2" is not convenient for server users.
```R
install.packages("devtools")
library(devtools)
install_github("cma2015/rsgcc")
```
- For users who would like to use GUI interface (**Windows system**)
```R
# install from cran
install.packages("rsgcc")
```
- For users who would like to use GUI interface (**ubuntu**)
```bash
sudo apt-get install libgtk2.0-dev
```
```R
install.packages("RGtk2")
install.packages("gWidgetsRGtk2")
install.packages("rsgcc")
```
- For users who would like to use GUI interface (**MAC OS X**)
**Step 1:** install R package [RGtk2](https://gist.github.com/sebkopf/9405675)
**Step 2:**
```R
install.packages("rsgcc")
```

## Quickly start
```R
library(rsgcc)
# calculate Gini correlation of two genes
data(rsgcc) #A gene expression matrix which contains 100 genes and 13 samples.
res <- cor.pair(idxvec = c(1,3), GEMatrix = rnaseq, rowORcol = "row", cormethod = "GCC")
res
```
**Note:** For more details about this package, please refer to https://github.com/cma2015/rsgcc/blob/master/rsgcc.pdf


## Citation
Ma, Chuang, and Xiangfeng Wang. "[**Application of the Gini correlation coefficient to infer regulatory relationships in transcriptome analysis.**](http://www.plantphysiol.org/content/early/2012/07/13/pp.112.201962.short)" Plant physiology (2012): pp-112.

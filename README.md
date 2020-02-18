## Background
Single cell RNA sequencing enables to obtain a transcriptome profiles for an individual cell. However, due to the technical limitations, the transcriptom profile through a single cell sequencing technology is vulnerable to a technical noise called a dropout event which results excessive zeros in the raw sequencing data. That is, the gene expression in the particular cells suffering dropout could severely deviate from the common gene expression patterns in the same cell type. The excessive zeros (or zero-inflated noise) limits the accurate and reliable downstream analysis of the single cell sequencing. PRIME (**PR**obabilistic **IM**putation to reduce dropout effects in a single cell **E**xpression) aims to effectively reduce dropout events in the single cell sequencing by removing the artificial zeros and keeping a biologically regulated zeros.

Overview of PRIME
![Alt text](./vignettes/Fig1.png?raw=true "Title")

<!-- PRIME: a probabilistic imputation method to reduce dropout effects in single cell RNA sequencing -->


## Installation guide 
** Note that since "devtools" has an issue with R version (>=3.6.0), please do the followin steps to install PRIME. Once "devtools" resolve the compatibility issue, we can provide the simple installation command. **

First, please do not use R version (>=3.6.0) to avoide compatibility issue. We tested R version 3.5.2.

Next, you need to install edgeR and multtest using the following commands:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("multtest")
BiocManager::install("edgeR")
```

Third, if you already installed Rtools, skip this step. Please install Rtools (ver 3.5) you can find Rtools at: https://cran.r-project.org/bin/windows/Rtools/

Next, please install caTools (ver 1.17.2) and install PRIME using the following commencs
```
install.packages("bitops")
pkgurl <- "https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.2.tar.gz"
install.packages(pkgurl , repos = NULL, type = "source")

install.packages("devtools")
library(devtools)
install_github("hyundoo/PRIME")
```



## Example
In this example, we utilize the test data from Buettner et al. (2015) and visualize each sample in a low-dimensional space through PCA. To run PRIME, first load R package and test data using the following commands:
```
library(PRIME)
data(testdata)
```
In order to impute dropouts through PRIME, simply run the following command:
```
PRIME_res <- PRIME(testdata)
```
PRIME accepts the count matrix, where each row corresponds to the genes and each column corresponds to the single cells. 

Note that it does no required to be normalized because PRIME() first normalizes the input data through cpm (counts per million) and transform it to log10 scale. PRIME() returns the same dimensional matrix.


To visualize the samples in a low dimensional space through PCA, use the following commands:
```
library(ggplot2)
pca_res <-prcomp(t(log10(1+PRIME_res)))
PC_all <- data.frame(pca_res$x, Sample=factor(label))
ggplot(PC_all, aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=1,alpha=0.7)
```

![Alt text](./vignettes/PCA_plot.png?raw=true "Title")

## Parameters
PRIME has several parameters and the default parameters result acceptable results for most cases. 
```
PRIME(sc_cnt, max_it = 5, err_max = 0.05, nPCs = 20, qth_min = 0.9, th_min = 0.85, min_nbr = 0.2, max_nbr = 1.25)
```
- sc_cnt: input count matrix, where each row correspondes to genes and each column corresponds to single cell samples
- max_it: The maximum number of iteration for imputing dropout events. PRIME generally runs 1 to 3 iterations.
- err_max: Stop criteria for each iteration. 
- nPCs: The number of principal components to estimate cell type similarities. Based on the emprical results, PRIME shows acceptable resutls for 15 to 20 PCs.
- qth_min: The minimum threshold to select neighboring nodes. 
- th_min: The weighting coefficient for percentile for similarity measurement (i.e., Pearson correlation and Euclidean distance). 
- min_nbr: Minimum percentage of neighboring nodes to identify local subnetwork (i.e., set of similar cells)
- max_nbr: Maximum percentage of neighboring nodes to identify local subnetwork (i.e., set of similar cells)


---
title: "Introduction to PRIME"
author: "Hyundoo Jeong and Zhandong Liu"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
## Background
Single cell RNA sequencing enables to obtain a transcriptome profiles for an individual cell. However, due to the technical limitations, the transcriptom profile through a single cell sequencing technology is vulnerable to a technical noise called a dropout event which results excessive zeros in the raw sequencing data. That is, the gene expression in the particular cells suffering dropout could severely deviate from the common gene expression patterns in the same cell type. The excessive zeros (or zero-inflated noise) limits the accurate and reliable downstream analysis of the single cell sequencing. PRIME (**PR**obabilistic **IM**putation to reduce dropout effects in a single cell **E**xpression) aims to effectively reduce dropout events in the single cell sequencing by removing the artificial zeros and keeping a biologically regulated zeros.


<!-- PRIME: a probabilistic imputation method to reduce dropout effects in single cell RNA sequencing -->


## Installation guide
To install PRIME, you need to install devtools and type the following command:
```
install.packages("devtools")
library(devtools)
install_github("hyundoo/PRIME")
```




## Example
In this example, we utilize the test data from Buettner et al. (2015). To run PRIME, the recommended input is simply the count matrix for the single cell RNA sequencing. Each row corresponds to the gene and each column corresponds to each cell. In the following example code, we visualize the single cells through PCA after imputing dropouts in the single cell RNA sequencing. 
```
library(PRIME)
data(testdata)

PRIME_res <- PRIME(testdata)

pca_res <-prcomp(t(log10(1+PRIME_res)))
plot(pca_res$x[,1], pca_res$x[,2], bg= c("red", "blue","green")[factor(label)],  type = "p", pch = 21, xlab="PC1", ylab="PC2")
legend("topleft", title="Cell types",
  	c("G1","S","G2M"), fill=c("red", "blue","green"), horiz=TRUE)
```

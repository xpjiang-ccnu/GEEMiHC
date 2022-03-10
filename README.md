# GEEMiHC

Type: Package

Title: Detecting sparse microbial association signals from longitudinal microbiome data based on generalized estimating equations

Version: 1.0

Author: Han Sun

Maintainer: Han Sun sunh529@mails.ccnu.edu.cn; Xingpeng Jiang xpjiang@mail.ccnu.edu.cn

Imports: phyloseq, cluster, compositions, permute, PGEE, vegan, ape, dirmult, aSPU, MiSPU, devtools 

Description: GEEMiHC is used for detecting sparse microbial association signals between microbiome and a host phenotype from longitudinal microbiome data.

License: GPL-2

Encoding: UTF-8

LazyData: true

URL: https://github.com/xpjiang-ccnu/GEEMiHC



## Introduction

This R package, **GEEMiHC**, can be used for detecting sparse microbial association signals adaptively from longitudinal microbiome data. It can be applied to datasets with diverse types of outcomes to study the association between diverse types of host phenotype and microbiome, such BMI (Gaussian distribution), disease status (Binomial distribution) or number of tumors (Poisson distribution). Considering cross-sectional data as a special case of longitudinal data, it can be also applied to cross-sectional data, in which case the results will be consistent with MiHC.



## Installation 

phyloseq:
```
BiocManager::install("phyloseq")
```

cluster:
```
install.packages("cluster")
```

compositions:
```
install.packages("compositions")
```

permute:
```
install.packages("permute")
```

PGEE:
```
install.packages("PGEE")
```

vegan:
```
install.packages("vegan")
```

ape:
```
install.packages("ape")
```

dirmult:
```
install.packages("dirmult")
```

aSPU:
```
install.packages("aSPU")
```

MiSPU:
```
install.packages("MiSPU")
```

devtools:
```
install.packages("devtools")
```


You may install `GEEMiHC` from GitHub using the following code: 

```
devtools::install_github("xpjiang-ccnu/GEEMiHC", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------



## Usage
```
GEEMiHC(y, id, covs, otu.tab, tree, model, Gamma=c(1,3,5,7,9), Lamda=matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T), comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
```



## Arguments
* _y_ - response variable (i.e., host phenotype of interest). Exponential family of distributions (e.g., Gaussian, Binomial, Poisson) outcomes.
* _id_ - A vector for identifying the sequence of subjects/clusters of longitudinal data.
* _covs_ - covariate (e.g., age, gender). Default is covs = NULL.
* _otu.tab_ - A matrix of the OTU table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
* _tree_ - A rooted phylogenetic tree. Default is tree = NULL, but we recommend the inclusion of phylogenetic tree information. If not, the weighted higher criticism tests cannot be considered.
* _model_ - "gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
* _Gamma_ - A subset consists of the candidate modulation schema for lower sparsity levels. Default is Gamma=c(1,3,5,7,9).
* _Lamda_ - The weight factor for candidate modulation schema in _Gamma_. Default is the equal weight for each candidate modulation schema, i.e., Lamda = matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T).
* _comp_ - An indicator if the OTU table contains absolute abundances or relative abundances. Default is comp=FALSE for absolute abundances.
* _CLR_ - An indicator if the OTU table needs to be converted using the centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
* _opt.ncl_ - A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
* _n.perm_ - A number of permutations. Default is n.perm=5000. 



## Values
_$rank.order_ - rank order for significant factors.

_$simes.pv.GEE.AR_  - The p-value for the Simes test of GEEMiHC with autoregressive structure.

_$simes.pv.GEE.EX_  - The p-value for the Simes test of GEEMiHC with exchange structure.

_$simes.pv.GEE.IN_  - The p-value for the Simes test of GEEMiHC with independence structure.

_$ind.pvs.GEEMiHC.AR_ - The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with autoregressive structure.

_$ind.pvs.GEEMiHC.EX_ - The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with exchange structure.

_$ind.pvs.GEEMiHC.IN_ - The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with independence structure.

_$ada.pvs_ - The p-values for global omnibus higher criticism tests of three GEEMiHC with different structure and aGEEMiHC.



## Example

Import requisite R packages: 

```
library(cluster)
library(permute)
library(phyloseq)
library(PGEE)  
library(GEEMiHC)  
```


Import example microbiome data:

```
data(CD_longitudinal)
otu.tab <- CD_longitudinal@otu_table
tree <- CD_longitudinal@phy_tree
y <- sample_data(CD_longitudinal)$label
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(CD_longitudinal)$age)
covs[,2] <- as.factor(sample_data(CD_longitudinal)$smoker)
id <- sample_data(CD_longitudinal)$id
```

Fit GEEMiHC:

```
set.seed(123)
out <- GEEMiHC(y, id, covs=covs, otu.tab=otu.tab, tree=tree, model="binomial", n.perm=5000)
out
```



## References
* Sun H, et al. Detecting sparse microbial association signals adaptively from longitudinal microbiome data based on generalized estimating equations. (under review)

* Koh H and Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. _Microbiome_ 2020;**8**(1):63.

* McMurdie PJ and Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;**8**(4):e61217

* Paradis E, et al. APE: Analyses of phylogenetics and evolution in R language. _Bioinformatics_ 2004;**20**(2):289-290.

* Reynolds A, et al. Clustering rules: A comparison of partitioning and hierarchical clustering algorithms. _J Math Model Algor_ 2006;**5**:475–504.

* Simes RJ. An improved Bonferroni procedure for multiple tests of significance. _Biometrika_ 1986;**73**(3):751-754.

* Vázquez-Baeza Y., et al. Guiding longitudinal sampling in IBD cohorts. _Gut_ 2018;**67**:1743-1745.

* Wang L. GEE analysis of clustered binary data with diverging number of covariates. _Ann. Statist._ 2011;**39**:389–417.

* Wang L, et al. Penalized generalized estimating equations for high-dimensional longitudinal data analysis. _Biometrics_ 2012;**68**(2):353-360.

* Wu C, et al. An adaptive association test for microbiome data. _Genome Med_ 2016;**8**(1):56.

---------------------------------------------------------------------------------------------------------------------------------------



# The function for generating simulated OTU count data

## SimulateOTU

### Description
We generate the OTUs count data simulated based on the Dirichlet-multinomial model according to real data.



### Usage
```
SimulateOTU(data, nSam, parameters, mu, size)
```



### Arguments
* _data_ - real data.

* _nSam_ - Sample size.

* _parameters_ - The estimated parameter based on a real microbiome data, including OTU proportions and overdispersion parameter.

* _mu_ - The mean of the negative binomial distribution.

* _size_ - The size of the negative binomial distribution.



### Values

_$OTU_ - OTU counts table simulated based on real data.



### Example

```
library(dirmult) 
data("throat.otu.tab", package = "MiSPU")
nOTU = 100
otu_sum <- apply(throat.otu.tab, 2, sum)
throat.otu.tab.100 <- throat.otu.tab[, order(otu_sum, decreasing = T)[1:nOTU]]
parameters <- dirmult(throat.otu.tab.100)
otu.tab <- SimulateOTU(throat.otu.tab.100, nSam = 50, parameters, mu = 1000, size = 25)
```



### References

* Chen J and Li H. Variable selection for sparse Dirichlet-multinomial regression with an application to microbiome data analysis. _Annals of Applied Statistics_ 2013;**7**(1).

* Sun H, et al. A powerful adaptive microbiome-based association test for microbial association signals with diverse sparsity levels. _Journal of Genetics and Genomics_ 2021;**48**(9):851-859.

* Wu C, et al. An adaptive association test for microbiome data. _Genome Med_ 2016;**8**(1):56.



## Statement

Our code mainly refers to R packages, _MiHC_, _MiSPU_ and _MiATDS_.

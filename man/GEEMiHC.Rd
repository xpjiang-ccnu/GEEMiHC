\name{GEEMiHC}
\alias{GEEMiHC}
\title{GEEMiHC}
\usage{
GEEMiHC(y, id, covs, otu.tab, tree, model, Gamma=c(1,3,5,7,9), Lamda=matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T), comp=FALSE, CLR=FALSE, opt.ncl=30, n.perm=5000)
}
\description{
This function tests longitudinal study of microbiome association test for sparse microbial association signals using aGEEMiHC. It can be applied to the simulation datasets with response variables of different traits and the real cross-sectional and longitudinal datasets to indicate its universality for various types of microbiome datasets.
}

\arguments{
  \item{y}{
Response variable (i.e., host phenotype of interest). Exponential family of distributions (e.g., Gaussian, Binomial, Poisson) outcomes.
}
  \item{id}{
A vector for identifying the sequence of subjects/clusters of longitudinal data.

}
  \item{covs}{
A data.frame (or matrix/vector) for covariate (e.g., age). Default is covs = NULL.
}
  \item{otu.tab}{
A matrix of the OTU table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
}
  \item{tree}{
A rooted phylogenetic tree. Default is tree = NULL, but we recommend the inclusion of phylogenetic tree information. If not, the weighted higher criticism tests cannot be considered.
}
  \item{model}{
"gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
}
  \item{Gamma}{
A subset consists of the candidate modulation schema for lower sparsity levels. Default is Gamma=c(1,3,5,7,9).
}

  \item{Lamda}{
An indicator to consider weighted high criticism tests or not. Default is W=TRUE to consider weighted higher criticism tests.
}
  \item{comp}{
The weight factor for candidate modulation schema in _Gamma_. Default is the equal weight for each candidate modulation schema, i.e., Lamda = matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T).
}
  \item{CLR}{
An indicator if the OTU table needs to be converted using the centered log-ratio (CLR) transformation. Default is CLR=FALSE for no CLR transformation.
}
  \item{opt.ncl}{
A upper limit to find the optimal number of clusters. Default is opt.ncl=30.
}
  \item{n.perm}{
A number of permutations. Default is n.perm=5000.
}
}

\value{
 rank.order: rank order for significant factors.

simes.pv.GEE.AR: The p-value for the Simes test of GEEMiHC with autoregressive structure.

simes.pv.GEE.EX: The p-value for the Simes test of GEEMiHC with exchange structure.

simes.pv.GEE.IN: The p-value for the Simes test of GEEMiHC with independence structure.

ind.pvs.GEEMiHC.AR: The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with autoregressive structure.

ind.pvs.GEEMiHC.EX: The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with exchange structure.

ind.pvs.GEEMiHC.IN: The p-values for the item-by-item unweighted and weighted higher criticism tests of GEEMiHC with independence structure.

ada.pvs: The p-values for global omnibus higher criticism tests of three GEEMiHC with different structure and aGEEMiHC.

}



\examples{
# Import requisite R packages
library(cluster)

library(phyloseq)
library(permute)
library(PGEE)
source("D:/R_pro/GEEMiHC/R/GEEMiHC.stat.R")


# Import example microbiome data
load(("D:/R_pro/GEEMiHC/data/CD_longitudinal.rda"))
otu.tab <- CD_longitudinal@otu_table
tree <- CD_longitudinal@phy_tree
y=CD_longitudinal@sam_data$label
covs <- data.frame(matrix(NA, length(y), 2))
covs[,1] <- as.numeric(sample_data(CD_longitudinal)$age)
covs[,2] <- as.factor(sample_data(CD_longitudinal)$smoker)
id <- sample_data(CD_longitudinal)$id

Lamda = matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = TRUE)

# Fit GEEMiHC
set.seed(123)
out <- GEEMiHC(y, id, covs = covs, otu.tab=otu.tab, tree = tree, model="binomial", Gamma = c(1,3,5,7,9), Lamda = Lamda, comp = FALSE, CLR = FALSE, opt.ncl = 30, n.perm = 5000)
out

}

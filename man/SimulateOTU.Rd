\name{SimulateOTU}
\alias{SimulateOTU}
\title{SimulateOTU}
\usage{
SimulateOTU(data, nSam, parameters, mu, size)
}
\description{
We generate the OTUs count table simulated based on the Dirichlet-multinomial model according to real data.
}

\arguments{
  \item{data}{Real data.
}
  \item{nSam}{Sample size.
}
  \item{parameters}{The estimated parameter based on a real microbiome data, including OTU proportions and overdispersion parameter.

}
  \item{mu}{The mean of the negative binomial distribution.
}
  \item{size}{The size of the negative binomial distribution.
}

}

\value{
OTU: OTU counts table simulated based on real data.

}



\examples{
#Simulated OTUs
library(dirmult)
data("throat.otu.tab", package = "MiSPU")
nOTU = 100
otu_sum <- apply(throat.otu.tab, 2, sum)
throat.otu.tab.100 <- throat.otu.tab[, order(otu_sum, decreasing = TRUE)[1:nOTU]]
parameters <- dirmult(throat.otu.tab.100)
otu.tab <- SimulateOTU(throat.otu.tab.100, nSam = 50, parameters=parameters, mu = 1000, size = 25)

}

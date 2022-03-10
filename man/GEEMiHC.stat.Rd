\name{GEEMiHC.stat}
\alias{GEEMiHC.stat}
\title{GEEMiHC.stat}
\usage{
GEEMiHC.stat(Zs, Gamma, Lamda, Ws=Ws)
}
\description{The calculation for the test statistic of GEEHC (i.e., uGEEHC and wGEEHC) and its modulated versions.

}

\arguments{
  \item{Zs}{The marginal standardized score statistic for each OTU.

}
  \item{Gamma}{A subset consists of the candidate modulation schema for lower sparsity levels. Default is Gamma=c(1,3,5,7,9).

}
  \item{Lamda}{The weight factor for candidate modulation schema in Gamma. Default is the equal weight for each candidate modulation schema, i.e., Lamda = matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T).

}
  \item{Ws}{The weight factor for each OTU.


}


\value{
maxHC: The test statistic for GEEHC and its modulated versions.
}


\examples{
GEEMiHC.stat(Zs=Zs, Gamma=c(1,3,5,7,9), Lamda=matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T), Ws=Ws)
}

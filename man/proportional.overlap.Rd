\name{proportional.overlap}
\alias{proportional.overlap}

\title{
Proportional similarity overlap index
}
\description{
This index is a measurement of the overlap of proportional similarity (PS) indices between 
pairs of taxa across samples. As the PS is a niche breadth index that measures the similarity
between a taxon's distribution and an environmental parameter, the PS overlap also considers
environmental information. A value of 1 indicates perfect overlap, and 0 indicates no overlap. 
}
\usage{
proportional.overlap(df, sampleInfo, envInfo, q = 1.65)
}

\arguments{
  \item{df}{
A matrix of taxa (rows) by samples (columns) as discrete counts per sample. Col 1 must be
a taxon identifier.
}
  \item{sampleInfo}{
A categorical variable identifying which samples correspond to which environment.
}
  \item{envInfo}{
A quantitative variable of some environmental parameter measured per sample.
}
  \item{q}{
A coefficient for the LOQ, set to 1.65. Decreasing or increasing q will affect which taxa
are flagged as being below the LOQ.
}
}
\details{
Taxon inputs must be as discrete counts across samples, with a taxon identifier in column 1. 
The output is a matrix of paired taxon overlap indices. Specifically, Proportional Overlap 
(PO) indices of PO[1,2] are rows and PO[2,1] are columns, for taxa 1 and 2, however these 
values are identical. Feinsinger's PS is calculated as per feinsingers.PS. 

Proportional overlap is calculated as follows:
								PO[i,j] = 1 - (|PS[i] - PS[j]|)/(PS[i] + PS[j])
whereby PS[i] is the PS of taxon i and PS[j] is the PS of taxon j. 

Taxa below the LOQ are flagged with an asterisk.

This function depends on reshape2 to produce a matrix of overlap comparisons.
}
\value{
An object of class "data.frame" that gives PO[1,2] as rows and PO[2,1] as columns.
}
\references{
Feinsinger et al. 1981. A simple measure of niche breadth. Ecology 62(1):27-32
}
\author{
Damien Finn
}
\note{
NA
}



\seealso{
feinsingers.PS
}
\examples{
data(df)
overlapdf <- df[48:53,]
sampleInfo <- c(rep("R1",10), rep("R2",10), rep("R3",10), rep("R4",10))
pH.grad <- c(2.1, 2.2, 2, 1.9, 2.1, 1.8, 1.9, 2, 2.1, 1.9, 3.5, 3.6, 3.5, 
             3.4, 3.6, 3.5, 3.5, 3.4, 3.7, 3.4, 6.6, 6.5, 6.4, 6.8, 7, 6.6, 
             6.8, 6.9, 7, 7.1, 8, 8.2, 7.9, 8.1, 7.8, 7.9, 8.3, 8.2, 8.1, 7.9) 
res <- proportional.overlap(overlapdf, sampleInfo, pH.grad)

}

\keyword{ niche }
\keyword{ feinsingers }

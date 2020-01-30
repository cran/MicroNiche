\name{levins.overlap}
\alias{levins.overlap}

\title{
Levin's niche overlap index
}
\description{
This index is a measurement of the overlap between pairs of taxa across samples. A value 
of 1 indicates perfect overlap, and 0 indicates no overlap.
}
\usage{
levins.overlap(df, q = 1.65)
}

\arguments{
  \item{df}{
A matrix of taxa (rows) by samples (columns) as discrete counts per sample. Col 1 must be
a taxon identifier.
}
  \item{q}{
A coefficient for the LOQ, set to 1.65. Decreasing or increasing q will affect which taxa
are flagged as being below the LOQ.
}
}
\details{
Taxon inputs must be as discrete counts across samples, with a taxon identifier in column 1.  
The output is a matrix of paired taxon overlap indices. Specifically, Levin's Overlap (LO) 
indices of LO[1,2] are rows and LO[2,1] are columns, for taxa 1 and 2. 

Levin's overlap is calculated as follows:
								LO[i,j] = sum(p[i]*p[j])/sum(p[i]^2)
whereby p[i] is the proportion of taxon i per sample and p[j] is the proportion of taxon j
per sample. 

Taxa below the LOQ are flagged with an asterisk.

This function depends on reshape2 to produce a matrix of overlap comparisons.
}
\value{
An object of class "data.frame" that gives LO[1,2] as rows and LO[2,1] as columns.
}
\references{
Ludwig and Reynolds. 1988 Statistical Ecology. Wiley and Sons, Milton Australia
}
\author{
Damien Finn
}
\note{
NA
}


\seealso{
levins.Bn
}
\examples{
data(df)
overlapdf <- df[48:53,]
res <- levins.overlap(overlapdf)

}

\keyword{ niche }
\keyword{ levins }

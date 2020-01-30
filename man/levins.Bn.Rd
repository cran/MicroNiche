\name{levins.Bn}
\alias{levins.Bn}

\title{
Levins' niche breadth index
}
\description{
This index is a measurement of niche breadth (Bn) of a taxon weighted by the number of differing
environments (R), whereby a value of 1 indicates generalist taxa that are equally abundant
across environments, and 0 indicates specialist taxa that favour specific environments.
}
\usage{
levins.Bn(df, R, sampleInfo, q = 1.65)
}

\arguments{
  \item{df}{
A matrix of taxa (rows) by samples (columns) as discrete counts per sample. Col 1 must be
a taxon identifier.
}
  \item{R}{
The number of different environments being compared.
}
  \item{sampleInfo}{
A categorical variable identifying which samples correspond to which environment.
}
  \item{q}{
A coefficient for the LOQ, set to 1.65. Decreasing or increasing q will affect which taxa
are flagged as being below the LOQ.
}
}
\details{
Taxon inputs must be as discrete counts across samples, with a taxon identifier in column 1.  
This function also performs null model testing against a distribution of Bn generated from 
R random proportions permuted 999 times. This allows the user to perform significance testing 
and to determine a probability that the taxa Bn differs from the null Bn. The output includes 
each taxon's Bn, a P value of significance testing against the null model distribution, and a 
Benjamin-Hochberg adjusted P value to account for false discovery rate. 

A plot of the null distribution with 0.05 and 0.95 quantiles highlighted as red dotted lines 
is also provided. 

Finally, this function also calculates the limit of quantification (LOQ) to identify taxa 
whose niche cannot be confidently measured. A plot of log abundance x taxon rank is provided
to visualise the distribution of the data, a fitted lognormal taxa rank model, and the 
LOQ determined based on the standard deviation of the lognormal model. 

Levin's niche breadth is calculated as follows:
								Bn[j] = (1/R)/sum(p[i]^2) 
whereby R is the number of differing environments and p[i] is the proportion of taxon j in
environment i. 

This function depends on ggplot2 to generate a plot of the null distribution and reshape2
to identify the taxa below the LOQ.
}
\value{
An object of class "data.frame" that gives the Bn, P value and adjusted P value for each taxon.
A column noting taxa that fall below the LOQ are identified as T.
}
\references{
Feinsinger et al. 1981. A simple measure of niche breadth. Ecology 62(1):27-32
}
\author{
Damien Finn
}
\note{
The null model testing becomes increasingly robust with increasing R. If R == 2 the results
should be interpreted with caution. 
}


\seealso{
levins.overlap
}
\examples{
data(df)
sampleInfo <- c(rep("R1",10), rep("R2",10), rep("R3",10), rep("R4",10))
res <- levins.Bn(df, 4, sampleInfo)
}

\keyword{ niche }
\keyword{ levins }

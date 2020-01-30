\name{feinsingers.PS}
\alias{feinsingers.PS}

\title{
Feinsinger's proportional similarity index
}
\description{
This index is a measurement of proportional similarity (PS) of a taxon in relation to an 
environmental parameter, such as pH or temperature. A value of 0 indicates an inverse 
relationship between taxon abundance and the environmental parameter, whilst 1 indicates a 
positive relationship between taxon abundance and the environmental parameter. A value of 
0.5 indicates no relationship.
}
\usage{
feinsingers.PS(df, R, sampleInfo, envInfo, q = 1.65)
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
This function also performs null model testing against a distribution of PS generated from 
R random proportions permuted 999 times. This allows the user to perform significance testing 
and to determine a probability that the taxa PS differs from the null PS. The output includes 
each taxon's PS, a P value of significance testing against the null model distribution, and 
a Benjamin-Hochberg adjusted P value to account for false discovery rate. 

A plot of the null distribution with 0.05 and 0.95 quantiles highlighted as red dotted lines 
is also provided. 

Finally, this function also calculates the limit of quantification (LOQ) to identify taxa 
whose niche cannot be confidently measured. A plot of log abundance x taxon rank is provided
to visualise the distribution of the data, a fitted lognormal taxa rank model, and the 
LOQ determined based on the standard deviation of the lognormal model.  

Feinsinger's PS is calculated as follows:
								PS[j] = 1 - 0.5*sum(|p[i] - r[i]|)
whereby p[i] is the proportion of taxon j in environment i and r[i] is the proportion of the
environmental parameter in i.

This function depends on ggplot2 to generate a plot of the null distribution and reshape2
to identify the taxa below the LOQ.
}
\value{
An object of class "data.frame" that gives the PS, P value and adjusted P value for each taxon.
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
proportional.overlap
}
\examples{
data(df)
sampleInfo <- c(rep("R1",10), rep("R2",10), rep("R3",10), rep("R4",10))
pH.grad <- c(2.1, 2.2, 2, 1.9, 2.1, 1.8, 1.9, 2, 2.1, 1.9, 3.5, 3.6, 3.5, 
             3.4, 3.6, 3.5, 3.5, 3.4, 3.7, 3.4, 6.6, 6.5, 6.4, 6.8, 7, 6.6, 
             6.8, 6.9, 7, 7.1, 8, 8.2, 7.9, 8.1, 7.8, 7.9, 8.3, 8.2, 8.1, 7.9) 
res <- feinsingers.PS(df, 4, sampleInfo, pH.grad)
}

\keyword{ niche }
\keyword{ feinsingers }

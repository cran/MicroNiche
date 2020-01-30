\name{df}
\alias{df}
\docType{data}
\title{
Data frame for package MicroNiche
}
\description{
This data frame consists of taxa (rows) by samples (columns). Taxa are separated into six
distributions (D1 - 6) with 10 'species' per distribution. Samples are separated into four
'environments' (R1 - 4) with 10 samples per environment.
}
\usage{data("df")}
\format{
  A data frame with 60 observations on the following 41 variables.
  \describe{
    \item{\code{Taxon}}{a factor with levels \code{D1S1} \code{D1S10} \code{D1S2} \code{D1S3} \code{D1S4} \code{D1S5} \code{D1S6} \code{D1S7} \code{D1S8} \code{D1S9} \code{D2S1} \code{D2S10} \code{D2S2} \code{D2S3} \code{D2S4} \code{D2S5} \code{D2S6} \code{D2S7} \code{D2S8} \code{D2S9} \code{D3S1} \code{D3S10} \code{D3S2} \code{D3S3} \code{D3S4} \code{D3S5} \code{D3S6} \code{D3S7} \code{D3S8} \code{D3S9} \code{D4S1} \code{D4S10} \code{D4S2} \code{D4S3} \code{D4S4} \code{D4S5} \code{D4S6} \code{D4S7} \code{D4S8} \code{D4S9} \code{D5S1} \code{D5S10} \code{D5S2} \code{D5S3} \code{D5S4} \code{D5S5} \code{D5S6} \code{D5S7} \code{D5S8} \code{D5S9} \code{D6S1} \code{D6S10} \code{D6S2} \code{D6S3} \code{D6S4} \code{D6S5} \code{D6S6} \code{D6S7} \code{D6S8} \code{D6S9}}
    \item{\code{R1S1}}{a numeric vector}
    \item{\code{R1S2}}{a numeric vector}
    \item{\code{R1S3}}{a numeric vector}
    \item{\code{R1S4}}{a numeric vector}
    \item{\code{R1S5}}{a numeric vector}
    \item{\code{R1S6}}{a numeric vector}
    \item{\code{R1S7}}{a numeric vector}
    \item{\code{R1S8}}{a numeric vector}
    \item{\code{R1S9}}{a numeric vector}
    \item{\code{R1S10}}{a numeric vector}
    \item{\code{R2S1}}{a numeric vector}
    \item{\code{R2S2}}{a numeric vector}
    \item{\code{R2S3}}{a numeric vector}
    \item{\code{R2S4}}{a numeric vector}
    \item{\code{R2S5}}{a numeric vector}
    \item{\code{R2S6}}{a numeric vector}
    \item{\code{R2S7}}{a numeric vector}
    \item{\code{R2S8}}{a numeric vector}
    \item{\code{R2S9}}{a numeric vector}
    \item{\code{R2S10}}{a numeric vector}
    \item{\code{R3S1}}{a numeric vector}
    \item{\code{R3S2}}{a numeric vector}
    \item{\code{R3S3}}{a numeric vector}
    \item{\code{R3S4}}{a numeric vector}
    \item{\code{R3S5}}{a numeric vector}
    \item{\code{R3S6}}{a numeric vector}
    \item{\code{R3S7}}{a numeric vector}
    \item{\code{R3S8}}{a numeric vector}
    \item{\code{R3S9}}{a numeric vector}
    \item{\code{R3S10}}{a numeric vector}
    \item{\code{R4S1}}{a numeric vector}
    \item{\code{R4S2}}{a numeric vector}
    \item{\code{R4S3}}{a numeric vector}
    \item{\code{R4S4}}{a numeric vector}
    \item{\code{R4S5}}{a numeric vector}
    \item{\code{R4S6}}{a numeric vector}
    \item{\code{R4S7}}{a numeric vector}
    \item{\code{R4S8}}{a numeric vector}
    \item{\code{R4S9}}{a numeric vector}
    \item{\code{R4S10}}{a numeric vector}
  }
}
\details{
An object of class "data.frame" of 60 rows x 40 columns, representing taxa and samples,
respectively.
}
\source{
Generated in silico to compare taxa of known distributions: D1 are generalists roughly
in equal abundance across the four environments. D2 decrease linearly across the four 
environments. D3 are specialists that decrease exponentially across the four environments.
D4 are specialists unique to one of the four environments. D5 are specialists in environments
two and four, with no relationship with an environmental gradient supplied in the vignette. 
Finally, D6 represent low, spurious counts haphazardly distributed across all samples often 
present in microbial count tables. 
}
\references{
NA
}
\examples{
data(df)
}
\keyword{datasets}

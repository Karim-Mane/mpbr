% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drop_samples.R
\name{remove_samples_from_vcf}
\alias{remove_samples_from_vcf}
\title{Remove samples from VCF file}
\usage{
remove_samples_from_vcf(vcf, samples_to_be_retained, path, index = 1L)
}
\arguments{
\item{vcf}{the input VCF file}

\item{samples_to_be_retained}{a vector of samples to be retained}

\item{path}{the path folder that contains the input data. This will be also
the folder where the output files will be stored.}

\item{index}{the index of the input VCF file.}
}
\value{
the path to filtered VCF file
}
\description{
Remove samples from VCF file
}
\keyword{internal}

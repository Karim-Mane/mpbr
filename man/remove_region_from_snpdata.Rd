% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drop_snps.R
\name{remove_region_from_snpdata}
\alias{remove_region_from_snpdata}
\title{Remove SNPs with a specified region of the genome.}
\usage{
remove_region_from_snpdata(snpdata, chrom, start, end)
}
\arguments{
\item{snpdata}{a \code{SNPdata} object}

\item{chrom}{the chromosome from which loci should be dropped}

\item{start}{the starting position of the region to be discarded}

\item{end}{the end position of the region to be discarded}
}
\description{
Remove SNPs with a specified region of the genome.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_genetics.R
\name{calculate_relatedness}
\alias{calculate_relatedness}
\title{Estimate the relatedness between pairs of isolates}
\usage{
calculate_relatedness(
  snpdata,
  mat_name = "Imputed",
  from = "Location",
  sweep_regions = NULL,
  groups = NULL
)
}
\arguments{
\item{snpdata}{SNPdata object}

\item{mat_name}{the name of the genotype table to be used. default="GT"}

\item{from}{the name of the column, in the metadata table, to be used
to represent the sample's population}

\item{sweep_regions}{a data frame with the genomic coordinates of the regions
of the genome to be discarded. The entire genome will be used if it is not
specified. This should contain the following 3 columns:
\enumerate{
\item Chrom: the chromosome ID
\item Start: the start position of the region on the chromosome
\item End:   the end position of the region on the chromosome
}}

\item{groups}{a vector of character. If specified, relatedness will be
estimated between those groups only}
}
\value{
\code{SNPdata} object with an extra field: \code{relatedness}. This will
contain the relatedness data frame of 3 columns and its correspondent
matrix
}
\description{
Estimate the relatedness between pairs of isolates
}
\details{
The relatedness calculation is based on the model developed by
Aimee R. Taylor and co-authors.
https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009101 # nolint: line_length_linter

\if{html}{\out{<div class="sourceCode">}}\preformatted{If the value for the `mat_name` argument is GT` or `Phased`, the missing
genotypes will be imputed prior to relatedness estimation. Similarly, if
`mat_name` is `Imputed`, the missing genotypes will be imputed first, if
the `Imputed` matrix does not exist in the SNPdata object.
}\if{html}{\out{</div>}}
}
\examples{
\dontrun{
  calculate_relatedness(
    snpdata,
    mat_name      = "GT",
    from          = "Location",
    sweep_regions = NULL,
    groups        = c("Chogen", "DongoroBa")
  )
 }

}

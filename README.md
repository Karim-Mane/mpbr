
<!-- README.md is generated from README.Rmd. Please edit that file. -->
<!-- The code to render this README is stored in .github/workflows/render-readme.yaml -->
<!-- Variables marked with double curly braces will be transformed beforehand: -->
<!-- `packagename` is extracted from the DESCRIPTION file -->
<!-- `gh_repo` is extracted via a special environment variable in GitHub Actions -->

# mpbr <img src="man/figures/logo.svg" align="right" width="120" />

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)
[![R-CMD-check](https://github.com/Karim-Mane/mpbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Karim-Mane/mpbr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Karim-Mane/mpbr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Karim-Mane/mpbr?branch=main)
[![lifecycle-concept](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-experimental.svg)](https://www.reconverse.org/lifecycle.html#concept)
<!-- badges: end -->

mpbr provides functions to perform population genomics analyses using
genome-wide Single Nucleotide Polymorphisms (SNPs) from the malaria
parasite *Plasmodium falciparum*.

<!-- This sentence is optional and can be removed -->

mpbr is developed at the [Medical Research Council, The Gambia Unit at
London School of Hygiene and Tropical
Medicine](https://www.lshtm.ac.uk/research/units/mrc-gambia) as part of
the [Malaria Population Genomics
(MPB)](https://data.org/initiatives/epiverse/) initiative.

## Installation

You can install the development version of mpbr from
[GitHub](https://github.com/) with:

``` r
pak::pak("Karim-Mane/mpbr")
library(mpbr)
```

``` r
library(mpbr)
```

## Create the `SNPdata` object

The functions in the mpbr require a **SNPdata** object. This is
generated with the `get_snpdata()` function.  
The function’s arguments and returned object are described in the
vignette and the function documentation.

``` r
snpdata <- get_snpdata(
  vcf_file   = system.file("extdata", "Input_Data.vcf.gz", package = "mpbr"), 
  meta_file  = system.file("extdata", "SampleMetadata.RDS", package = "mpbr"), 
  output_dir = tempdir(), 
  gof        = system.file("extdata", "pf_gene_ontology.RDS", package = "mpbr"), 
  gff        = system.file("extdata", "PlasmoDB-56_Pfalciparum3D7.RDS",
                           package = "mpbr")
)
```

## Package Vignettes

The vignette of the mpbr contains detailed illustrations about the use
of each function. This can be accessed by typing the command below:

``` r
# OPEN THE VIGNETTE WITHIN RSTUDIO
vignette("Summary statistics")
vignette("Population genetics analysis")

# OPEN THE VIGNETTE IN YOUR WEB BROWSER.
browseVignettes("mpbr")
```

## Development

### Lifecycle

This package is currently a *concept*, as defined by the [RECON software
lifecycle](https://www.reconverse.org/lifecycle.html). This means that
essential features and mechanisms are still being developed, and the
package is not ready for use outside of the development team.

### Contributions

Contributions are welcome via [pull
requests](https://github.com/Karim-Mane/mpbr/pulls).

### Code of Conduct

Please note that the mpbr project is released with a [Contributor Code
of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

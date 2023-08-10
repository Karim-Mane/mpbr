
<!-- README.md is generated from README.Rmd. Please edit that file. -->
<!-- The code to render this README is stored in .github/workflows/render-readme.yaml -->
<!-- Variables marked with double curly braces will be transformed beforehand: -->
<!-- `packagename` is extracted from the DESCRIPTION file -->
<!-- `gh_repo` is extracted via a special environment variable in GitHub Actions -->

# {{ packagename }} <img src="man/figures/logo.svg" align="right" width="120" />

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)
[![R-CMD-check](https://github.com/%7B%7B%20gh_repo%20%7D%7D/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/%7B%7B%20gh_repo%20%7D%7D/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/%7B%7B%20gh_repo%20%7D%7D/branch/main/graph/badge.svg)](https://app.codecov.io/gh/%7B%7B%20gh_repo%20%7D%7D?branch=main)
[![lifecycle-concept](https://raw.githubusercontent.com/reconverse/reconverse.github.io/master/images/badge-experimental.svg)](https://www.reconverse.org/lifecycle.html#concept)
<!-- badges: end -->

{{ packagename }} provides functions to perform population genomics
analyses using genome-wide Single Nucleotide Polymorphisms (SNPs) from
the malaria parasite *Plasmodium falciparum*.

<!-- This sentence is optional and can be removed -->

{{ packagename }} is developed at the [Medical Research Council, The
Gambia Unit at London School of Hygiene and Tropical
Medicine](https://www.lshtm.ac.uk/research/units/mrc-gambia) as part of
the [Malaria Population Genomics
(MPB)](https://data.org/initiatives/epiverse/) initiative.

## Installation

You can install the development version of {{ packagename }} from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("{{ gh_repo }}")
```

## Create the `SNPdata` object

The functions in the {{ packagename }} require a **SNPdata** object.
This can be generated with the `get_snpdata()` function.  
The function’s arguments and returned object are described in the
vignette and the function documentation.

``` r
snpdata <- get_snpdata(
  vcf_file = system.file("extdata", "file.vcf.gz", package = "mpbr"), 
  meta_file = system.file("extdata", "metadata.RDS", package = "mpbr"), 
  output_dir = "/path/to/output/dir", 
  gaf = system.file("extdata", "Pf_gene_ontology.RDS", package = "mpbr"), 
  gff = system.file("extdata", "PlasmoDB-56_Pfalciparum3D7.gff", 
                    package = "mpbr")
)
```

## Package Vignettes

The vignette of the {{ packagename }} contains detailed illustrations
about the use of each function. This can be accessed by typing the
command below:

``` r
# OPEN THE VIGNETTE WITHIN RSTUDIO
vignette("mpbr")

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
requests](https://github.com/%7B%7B%20gh_repo%20%7D%7D/pulls).

### Code of Conduct

Please note that the {{ packagename }} project is released with a
[Contributor Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

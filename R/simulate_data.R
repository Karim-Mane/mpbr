#' Simulate data given parameters, the genotype matrix fs and the genetic
#' distance between loci.
#'
#' @param fs a matrix of genotype data with N rows and M columns where N is the
#'      number of loci and M the number of possible allele types.
#'      Example: if 3 types are possible at position p, but 5 types are possible
#'      somewhere else, then fs should have five columns, and fs[p,] might look
#'      like [0.2, 0.3, 0.5, 0, 0]
#' @param gendist a vector of genetic distance between loci where gendist[p]
#'      contains the distance between position p and p+1 or equivalently,
#'      gendist[p-1] contains the distance between position p-1 and p, for p > 1
#' @param k the k
#' @param r the r
#' @param epsilon the epsilon
#' @param rho the rho
#'
#' @return a matrix of 2 columns and N rows where N is the number of loci in the
#'      genotype matrix.
#' @keywords internal
#'
simulate_data <- function(fs, gendist, k, r,
                          epsilon = 0.001,
                          rho     = 7.4 * 10^(-7)) { # nolint: implicit_integer_linter
  ndata      <- dim(fs)[[1L]]
  maxnstates <- dim(fs)[[2L]]
  nstates    <- 0L
  Ys         <- matrix(NA, nrow = ndata, ncol = 2L) # nolint: object_name_linter
  for (position in seq_len(ndata)) {
    if (position == 1L) {
      IBD_current <- (runif(1L) <= r) # nolint: object_name_linter
    } else {
      if (IBD_current) {
        IBD_current <- (runif(1L) < (1L - (1L - r) * (1 - exp(-k * rho * gendist[position - 1L])))) # nolint: line_length_linter
      } else {
        IBD_current <- (runif(1L) < r*(1L - exp(-k * rho * gendist[position - 1L]))) # nolint: line_length_linter
      }
    }
    # number of possible allele types at current position
    nstates <- 1L
    while ((nstates < maxnstates) && (fs[position, nstates] > 1e-10)) {
      nstates <- nstates + 1L
    }
    #
    Gi <- sample(x    = 0L:(nstates - 1L), # nolint: object_name_linter
                 size = 1L,
                 prob = fs[position, 1L:nstates])
    Gj <- NA # nolint: object_name_linter
    Yi <- NA # nolint: object_name_linter
    Yj <- NA # nolint: object_name_linter
    # generate Gi, Gj given IBD_current
    if (IBD_current) {
      Gj <- Gi # nolint: object_name_linter
    } else {
      Gj <- sample(x    = 0L:(nstates - 1L), # nolint: object_name_linter
                   size = 1L,
                   prob = fs[position, 1L:nstates])
    }
    # generate Yi, Yj given Gi, Gj,
    # i.e. genotyping error part of the model
    if (runif(1L) < (1L - (nstates - 1L) * epsilon)) {
      Yi <- Gi # nolint: object_name_linter
    } else {
      # sample uniformly on the other possible states
      otherstates <- setdiff(0L:(nstates - 1L), Gi)
      Yi <- sample(x = otherstates, size = 1L) # nolint: object_name_linter
    }
    if (runif(1L) < (1L - (nstates - 1L) * epsilon)) {
      Yj <- Gj # nolint: object_name_linter
    } else {
      otherstates <- setdiff(0L:(nstates - 1L), Gj)
      Yj <- sample(x = otherstates, size = 1L) # nolint: object_name_linter
    }
    Ys[position, ] <- c(Yi, Yj) # nolint: object_name_linter
  }
  return(Ys)
}

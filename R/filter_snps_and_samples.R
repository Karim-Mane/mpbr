#' Filter loci and samples
#'
#' This function filters the SNPs and samples based on the user-specified
#' conditions. It performs a soft filtering by default i.e. it allows users to
#' get the proportion of remaining SNPs and samples across different cut-offs.
#' Users can activate hard filtering by setting \code{soft = FALSE}.
#'
#' @param snpdata An object of class \code{SNPdata}
#' @param soft A boolean that specifies whether to perform iterative filtering
#'    or not. Iterative filtering consists in soft filtering the input object by
#'    varying the missingness cut-off within the specified range of the
#'    \code{cut_off} argument, and the MAF cut-off from 0.1 to 0.4.
#'    Default is \code{TRUE}.
#' @param missingness A numeric vector of three elements that specifies the
#'    cut-off range to use for the soft filtration. When provided, SNPs and
#'    samples will be filtered across this range. Default is
#'    \code{missingness = c(from = 0.1, to = 0.9, by = 0.1)}.
#' @param maf A numeric vector of three elements that specifies the
#'    cut-offs to use for the soft filtration on MAF. When provided, SNPs will
#'    be filtered across this range. Default is
#'    \code{maf = c(from = 0.01, to = 0.5, by = 0.01)}.
#' @param min_qual A numeric value representing the minimum call quality score
#'    below which a loci will be discarded. Default is 1000.
#' @param max_missing_sites A numeric value between `0` and `1` representing the
#'    maximum fraction of missing sites above which a sample should be
#'    discarded. Default is 0.2.
#' @param max_missing_samples A numeric value between `0` and `1` representing
#'    the maximum fraction of missing samples above which a loci should be
#'    discarded. Default is 0.2.
#' @param maf_cutoff A numeric value between `0` and `1` representing the MAF
#'    cut-off. Loci with \code{MAF < maf_cutoff} will be discarded.
#'    Default is 0.01.
#'
#' @return When \code{soft = TRUE}, the function returns a list with the
#'    following two elements of type data frame:
#'    \describe{
#'        \item{`missingness`}{A data frame with three columns showing the
#'        remaining proportion of samples and loci after every missingness
#'        cut-off.}
#'        \item{`maf`}{A data frame with two columns representing the remaining
#'        proportion of loci at every MAF cut-off}
#'    }
#' 
#'    When \code{soft = FALSE}, the function returns a subset of the input
#'    \code{SNPdata} object, where loci and samples that do not meet the
#'    specified criteria have been removed.
#'
#' @examples
#' \dontrun{
#'   # perform hard filtering
#'   snpdata <- filter(
#'     snpdata,
#'     soft = FALSE,
#'     min_qual = 10,
#'     max_missing_sites = 0.2,
#'     max_missing_samples = 0.2,
#'     maf_cutoff = 0.01
#'   )
#'  
#'   # perform soft filtering
#'   snpdata <- filter(
#'     snpdata,
#'     soft = TRUE,
#'     missingness = c(from = 0.1, to = 0.9, by = 0.1),
#'     maf = c(from = 0.01, to = 0.5, by = 0.01)
#'   )
#' }
#' @export
#'
filter <- function(snpdata,
                   soft = TRUE,
                   missingness = c(from = 0.1, to = 0.9, by = 0.1),
                   maf = c(from = 0.01, to = 0.5, by = 0.01),
                   min_qual = 1000L,
                   max_missing_sites = 0.2,
                   max_missing_samples = 0.2,
                   maf_cutoff = 0.01) {
  checkmate::assert_vector(missingness, len = 3, any.missing = FALSE,
                           null.ok = FALSE)
  checkmate::checkNames(missingness, type = "named",
                        identical.to = c("from", "to", "by"))
  checkmate::assert_vector(maf, len = 3, any.missing = FALSE,
                           null.ok = FALSE)
  checkmate::checkNames(maf, type = "named",
                        identical.to = c("from", "to", "by"))
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_qual, lower = 10L, any.missing = FALSE,
                            null.ok = TRUE)
  checkmate::assert_numeric(max_missing_sites, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = TRUE,
                            len = 1L)
  checkmate::assert_numeric(max_missing_samples, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = TRUE,
                            len = 1L)
  checkmate::assert_logical(soft, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(maf_cutoff, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = TRUE,
                            len = 1L)
  
  # allow for recursive filtering
  if (soft) {
    return(
      filter_soft(snpdata, missingness, maf)
    )
  } else {
    return(
      filter_hard(snpdata, min_qual, maf_cutoff, max_missing_samples,
                  max_missing_sites)
    )
  }
}


#' Soft filter on SNPs and samples
#'
#' @inheritParams filter
#' 
#' @return A list with the following two elements of type data frame:
#' \describe{
#'    \item{`missingness`}{A data frame with three columns showing the remaining
#'    proportion of samples and loci after every missingness cut-off.}
#'    \item{`maf`}{A data frame with two columns representing the remaining
#'    proportion of loci at every MAF cut-off}.
#' }
#'
#' @keywords internal
#'
filter_soft <- function(snpdata, missingness, maf) {
  # soft filtering snps and samples by varying their missingness cut-off from
  # 0.9 until 0.1 by a step of 0.1
  details <- snpdata[["details"]]
  meta <- snpdata[["meta"]]
  
  # create the table where to store the filtration outcome at every cut-off
  # we fill the first row with the number of samples and loci before filtration
  # starts
  missingness_filtration_res <- data.frame(
    missingness = character(),
    remaining_loci = numeric(),
    remaining_samples = numeric()
  )
  
  # we set the filtration cut-offs from the default or the one provided by the
  # user
  missingness_cutoff <- seq(
    missingness[["from"]], missingness[["to"]], missingness[["by"]]
  )

  # perform soft filtering on missingness
  for (c in missingness_cutoff) {
    remaining_snps <- round(
      length(which(details[["percentage_missing_samples"]] <= c)) /
        nrow(details),
      digits = 2
    )
    remaining_samples <- round(
      length(which(meta[["percentage_missing_sites"]] <= c)) / nrow(meta),
      digits = 2
    )
    temp_res <- data.frame(
      missingness = paste0("≥", c),
      remaining_loci = remaining_snps,
      remaining_samples = remaining_samples
    )
    missingness_filtration_res <- rbind(missingness_filtration_res, temp_res)
  }

  # perform soft filtering of snps on MAF
  # cut-off is varied from 0.01 to 0.5 by a step of 0.01 by default
  maf_cutoff <- seq(maf[["from"]], maf[["to"]], maf[["by"]])
  
  # the maf filtration outcome will be stored in the data frame below
  maf_filtration_res <- data.frame(
  maf = character(),
  remaining_loci = numeric()
  )
  
  # perform the MAF filtration
  for (c in maf_cutoff) {
    remaining_snps <- round(
      length(which(details[["MAF"]] >= c)) / nrow(details),
      digits = 2
    )
    temp_res <- data.frame(
      maf = paste0("≥", c),
      remaining_loci = remaining_snps
    )
    maf_filtration_res <- rbind(maf_filtration_res, temp_res)
  }
  
  # return a list with the 2 outputs
  # users might need these results for reporting purpose
  return(list(
    missingness = missingness_filtration_res,
    maf = maf_filtration_res
  ))
}



#' Hard filter on SNPs and samples
#'
#' @inheritParams filter
#'
#' @return A subset of the input \code{SNPdata} object where samples and SNPs
#'    that do not meet the defined thresholds have been remove.
#' @keywords internal
#'
filter_hard <- function(snpdata, min_qual, maf_cutoff, max_missing_samples,
                        max_missing_sites) {
  meta <- snpdata[["meta"]]
  details <- snpdata[["details"]]
  
  # determine loci that meet the criteria specified by the user
  passed_filters <- details[["Qual"]] >= min_qual &
    details[["percentage_missing_samples"]] <= max_missing_samples &
    details[["MAF"]] >= maf_cutoff

  # throw an error if there is SNP that passes the defined filters
  stopifnot("No locus in VCF file has satisfied the specified QC metrics." =
              sum(passed_filters) > 0L)
  
  # send a message if all SNPs have passed the defined filters
  # When there are some loci that do not meet the defined QC metrics, remove
  # them from the details table and all the genotype matrices.
  if (sum(passed_filters) == nrow(details)) {
    message("All loci have satisfied the specified QC metrics.")
  } else {
    fields <- c("GT", "Phased", "Phased_Imputed")
    details <- details[passed_filters, ]
    for (field in fields) {
      if (field %in% names(snpdata)) {
        snpdata[[field]] <- snpdata[[field]][passed_filters, ]
      }
    }
  }
  
  # filter out samples that do not meet the required criteria
  # or send a message if all samples have met the define criteria
  passed_filters <- meta[["percentage_missing_sites"]] <= max_missing_sites
  stopifnot("No sample has satisfied the specified QC metrics." =
              sum(passed_filters) > 0L)
  if (sum(passed_filters) == nrow(meta)) {
    message("All samples have satisfied the specified QC metrics.")
  } else {
    did_not_passed <- meta[["sample"]][!passed_filters]
    message("The following samples will be removed: ", toString(did_not_passed))
    meta <- meta[passed_filters, ]
    
    # filter out those samples from the genotype matrices
    fields <- c("GT", "Phased", "Phased_Imputed")
    for (field in fields) {
      if (field %in% names(snpdata)) {
        idx <- which(colnames(snpdata[[field]]) %in% did_not_passed)
        snpdata[[field]] <- snpdata[[field]][, -idx]
      }
    }
    
    # recalculate the missingness and MAF after filtration
    details[["MAF"]] <- NULL
    details[["percentage_missing_samples"]] <-
      rowSums(is.na(snpdata[["GT"]])) / ncol(snpdata[["GT"]])
    meta[["percentage_missing_sites"]] <-
      colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
    snpdata[["details"]] <- details
    snpdata[["meta"]] <- meta
    snpdata <- compute_maf(snpdata, include_het = FALSE, mat_name = "GT")
  }

  return(snpdata)
}


#' Plot the SNPs distribution across the genome
#'
#' @param snpdata An object of class \code{SNPdata}
#' @param window_size An integer that is used to specify the window size to
#'    consider. The chromosomes will be split by this window size and the
#'    number of SNPs per window is evaluated then plotted.
#'
#' @return Displays the distribution of the SNPs across the genome.
#' @export
#'
#' @examples
plot_distribution = function(snpdata, window_size = 100){
  details <- snpdata[["details"]]
  details[["Pos"]] <- as.numeric(details[["Pos"]])

  # set the window size at 10Kb
  window_size <- window_size
  
  # count SNPs per window for each chromosome
  snp_density <- details %>%
    dplyr::mutate(
      window = cut(
        Pos,
        breaks = seq(0, max(Pos), by = window_size),
        include.lowest = TRUE
      )
    ) %>%
    dplyr::group_by(Chrom, window) %>%
    dplyr::summarise(snp_count = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate( #get window start and stop
      window_start = as.numeric(
        gsub("[^0-9.e+-,.]", "",
             do.call(data.frame, strsplit(as.character(window), ","))[1,])
      ), 
      window_end = window_start + window_size
    ) 
  
  snp_density[["Chrom"]] <- as.numeric(substr(snp_density[["Chrom"]], 7, 8))
  
  # Create the SNP density plot
  p <- ggplot(snp_density, aes(x = window_start, y = snp_count)) +
    geom_bar(stat = "identity", position = "stack", fill = "steelblue",
             colour = "darkblue") +  
    facet_grid(rows = vars(Chrom), scales = "free_x") + # Correctly facet by row
    theme_minimal(base_size = 14) + # Increase font size for better readability
    labs(x = "Genomic Location (bp)", 
         y = "SNP Density", 
         title = "SNP Density Plot by Chromosome") +
    theme(
      # axis.text.x = element_blank(),  # Remove x-axis labels
      # axis.ticks.x = element_blank(),         # Remove x-axis ticks
      strip.text.y = element_text(size = 12), # Customize label size for y-axis
      axis.text.y = element_blank(),           # Remove y-axis labels
      axis.ticks.y = element_blank()           # Remove y-axis ticks
    )
  print(p)
}

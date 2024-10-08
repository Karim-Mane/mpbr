#' Filter loci and samples (requires **bcftools** and **tabix** to be installed)
#'
#' This function filters the SNPs and samples based on the specified conditions
#'
#' @param snpdata An object of class \code{SNPdata}
#' @param min_qual the minimum call quality score below which a loci will be
#'    discarded. default = 10
#' @param max_missing_sites the maximum fraction of missing sites above which
#'    a sample should be discarded. default = 0.2
#' @param max_missing_samples the maximum fraction of missing samples above
#'    which a loci should be discarded. default = 0.2
#'  @param soft A boolean that specify whether to perform iterative filtering
#'    or not. Iterative filtering consists in soft filtering the input object by
#'    varying the missingness cut-off from 0.8 to 0.1, and the MAF cut-off from
#'    0.4 to 0.1. Default is \code{FALSE}.
#'    @param hard A boolean that specify whether to perform hard iterative filtering
#'    or not. The hard Iterative filtering is done by filtering the input object by
#'    varying the missingness cut-off from 0.8 to 0.1, and filters out the inputted MAF and min_qual.
#'     Default is \code{FALSE}.
#' @param maf_cutoff the MAF cut-off. loci with a MAF < maf_cutoff will be
#'    discarded. default = 0.01
#'
#' @return a filtered SNPdata object
#' @examples
#' \dontrun{
#'  snpdata <- filter_snps_samples(
#'   snpdata,
#'   soft                = FALSE,
#'   hard                = FALSE,
#'   min_qual            = 10,
#'   max_missing_sites   = 0.2,
#'   max_missing_samples = 0.2,
#'   maf_cutoff          = 0.01
#'  )
#'  }
#' @export
#'
filter_snps_samples <- function(snpdata,
                                soft                = FALSE,
                                hard                = FALSE,
                                min_qual            = 10L,
                                max_missing_sites   = 0.2,
                                max_missing_samples = 0.2,
                                maf_cutoff          = 0.01) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_qual, lower = 10L, any.missing = FALSE,
                            null.ok = FALSE)
  checkmate::assert_numeric(max_missing_sites, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(max_missing_samples, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_logical(soft, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(maf_cutoff, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  
  # allow for recursive filtering
  if (soft) {
    filter_soft(snpdata[["details"]], snpdata[["meta"]])
  }else if (hard){
    message("Conducting hard filteration")
    filter_hard(snpdata, min_qual, maf_cutoff)
  } else {
    x      <- snpdata[["details"]]
    fields <- c("GT", "Phased", "Phased_Imputed")
    idx    <- which(x[["Qual"]] >= min_qual &
                      x[["percentage_missing_samples"]] <= max_missing_samples &
                      x[["MAF"]] >= maf_cutoff)
    stopifnot("\nNo locus in VCF file has satisfied the specified QC metrics." =
                length(idx) != 0L)
    if (all(length(idx) > 0L && length(idx) < nrow(snpdata[["details"]]))) {
      x    <- x[idx, ]
      snpdata[["details"]] <- x
      for (field in fields) {
        if (field %in% names(snpdata)) {
          snpdata[[field]] <- snpdata[[field]][idx, ]
        }
      }
      Chrom <- Pos <- NULL # nolint: object_name_linter
      f2c          <- x %>% dplyr::select(Chrom, Pos)
      output_dir   <- dirname(snpdata[["vcf"]])
      data.table::fwrite(f2c, file.path(output_dir, "loci_to_be_retained.txt"),
                         col.names = FALSE,
                         row.names = FALSE,
                         quote     = FALSE,
                         sep       = "\t",
                         nThread   = 4L)
      snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]],
                                               "loci_to_be_retained.txt",
                                               output_dir,
                                               index = snpdata[["index"]])
    } else if (length(idx) == nrow(snpdata[["details"]])) {
      message("All loci have satisfied the specified QC metrics.")
    }
    
    idx <- which(snpdata[["meta"]][["percentage_missing_sites"]] <=
                   max_missing_sites)
    stopifnot("\nNo sample in VCF file has satisfied the specified QC metrics." =
                length(idx) != 0L)
    if (all(length(idx) > 0L && length(idx) < nrow(snpdata[["meta"]]))) {
      message("\nthe following samples will be removed:\n",
              paste(snpdata[["meta"]][["sample"]][-idx], collapse = "\n"))
      snpdata[["meta"]] <- snpdata[["meta"]][idx, ]
      data.table::fwrite(snpdata[["meta"]][["sample"]][-idx],
                         file.path(output_dir, "samples_to_be_dropped.txt"),
                         col.names = FALSE, row.names = FALSE,
                         quote = FALSE, sep = "\t", nThread = 4L)
      snpdata[["vcf"]] <- remove_samples_from_vcf(snpdata[["vcf"]],
                                                  "samples_to_be_dropped.txt",
                                                  output_dir,
                                                  index = snpdata[["index"]])
    } else if (length(idx) == nrow(snpdata[["meta"]])) {
      message("\nAll samples have satisfied the specified QC metrics.")
    }
    
    snpdata[["index"]] <- snpdata[["index"]] + 1L
    plot_snpDistribution(snpdata[["details"]])
    
  }
  
  
  snpdata
}


#' Soft filter on SNPs and samples
#'
#' @param details The \code{SNPdata} object \code{details} data frame.
#' @param meta The \code{SNPdata} object \code{meta} data frame.
#'
#' @return Displays the number of SNPs and samples that will be left after
#'    applying the different cut-offs (from 0.8 to 0.1 for missingness, and
#'    0.4 to 0.1 for MAF).
#' @keywords internal
#'
filter_soft <- function(details, meta) {
  # soft filtering snps and samples by varying their missingness cut-off from
  # 0.8 until 0.1 by a step of 0.1
  num_sampleLoci =data.frame(matrix(rbind(c(nrow(meta), "samples", 1), c(nrow(details), "loci", 1)), ncol = 3))
  colnames(num_sampleLoci) = c("length", "variable", "missingness")
  missingness_cutoff <- seq(0.8, 0.1, -0.1)
  message("soft filter on missingness")
 
  for (c in missingness_cutoff) {
    idx_snps <- length(which(details[["percentage_missing_samples"]] <= c))
    idx_sample <- length(which(meta[["percentage_missing_sites"]] <= c))
    num_sampleLoci = rbind.data.frame(num_sampleLoci, c(idx_sample, "samples", c), c(idx_snps, "loci", c))
    cli::cli_alert_info("cut-off = {c}: remaining {idx_snps} SNPs and \\
                        {idx_sample} samples.")
    
  }
  
  # Change the colors manually
  p <- ggplot(data = num_sampleLoci, aes(x = missingness, y = length, fill = variable)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    theme_minimal() +
    geom_text(aes(label = length), hjust = -0.5, color = "black",
              position = position_dodge(0.9), size = 3) +
    scale_fill_manual(values = c('#999999', 'grey25')) +
    labs(x = "Percentage Missingness", 
         y = "Number of Variables", 
         title = "Horizontal Bar Plot of Length by Percentage Missingness") +
    coord_flip()  # Flip the axes to create a horizontal bar plot
    #theme_classic()
  # Use custom colors
  print(p)
  
  # soft filtering snps by varying their MAF cut-off from
  # 0.4 until 0.1 by a step of 0.1
  num_sampleLoci1 =data.frame(matrix(rbind(c(nrow(details), 0)), ncol = 2))
  colnames(num_sampleLoci1) = c("length", "MAF")
  message("\nsoft filter on MAF")
  maf_cutoff <- seq(0.05, 0.01, -0.01)
  for (c in maf_cutoff) {
    idx_snps <- length(which(details[["MAF"]] >= c))
    num_sampleLoci1 = rbind.data.frame(num_sampleLoci1, c(idx_snps, c))
    
    cli::cli_alert_info("cut-off = {c}: remaining {idx_snps} SNPs.")
  }
  
  pp = ggplot(num_sampleLoci1, aes(x=MAF, y=length)) +
    geom_bar(stat="identity")+theme_minimal()+
    geom_text(aes(label=length), vjust=-0.5, color="darkblue", size=4, position="stack")+
    labs(x = "Minor Allele Frequency (MAF)", 
         y = "number of loci", 
         title = "Bar Plot of Length by Minor Allele Frequency") #+  # Add title
    #scale_x_discrete(breaks = num_sampleLoci1$MAF)
  print(pp)
}


#fwrite(as.data.frame(snpdata[["details"]]), paste0(path, "/details.txt"), sep = "\t")


#
filter_hard <- function(snpdata, min_qual, maf_cutoff) {
  # hard filtering snps and samples by varying their missingness cut-off from
  # 0.8 until 0.1 by a step of 0.1
  meta = snpdata[["meta"]]
  data = snpdata[["GT"]]
  details = as.data.frame(snpdata[["details"]])
  vcf_file = snpdata[["vcf"]]
  pos_seq = paste0(details[,1], ":",details[,2])
  rownames(data) = pos_seq
  
  missingness_cutoff <- seq(0.8, 0.1, -0.1)
  
 for (i in missingness_cutoff) {
   lociMissing = calculate_missingness(data)
   data = data[lociMissing < i,]
   sampleMissing = calculate_missingness(t(data))
   data = data[,sampleMissing < i] 
   }
  
  
  
  details = details[(pos_seq %in% rownames(data) & details[["Qual"]] > min_qual & details[["MAF"]] > maf_cutoff),] 
  snpdata[["details"]] = details
  snpdata[["meta"]] = meta[meta$sample %in% colnames(data),]
  snpdata[["GT"]] = data[(rownames(data) %in% paste0(details[["Chrom"]], ":", details[["Pos"]])),]
  #**filter snps and position from the vcf
  data.table::fwrite(as.data.frame(snpdata[["details"]]), paste0(path, "/details.txt"))
  plot_snpDistribution(details)
  #apply filteration on list
}

plot_snpDistribution = function(details){
  library(ggplot2)
  details[["Pos"]] = as.numeric(details[["Pos"]])
  # Define window size (10Kb)
  window_size <- 1000
  # Count SNPs per window for each chromosome
  snp_density <- details %>%
    mutate(Window = cut(Pos, breaks = seq(0, max(Pos), by = window_size), include.lowest = TRUE)) %>%
    group_by(Chrom, Window) %>%
    summarise(SNP_count = n(), .groups = 'drop') %>%
    mutate(Window_start = as.numeric(gsub("[^0-9.e+-,.]", "", do.call(data.frame, strsplit(as.character(Window), ","))[1,])), 
           Window_end = Window_start + window_size) #get window start and stop
  
  snp_density[["Chrom"]] = as.numeric(substr(snp_density[["Chrom"]], 7, 8))
  # Create the SNP density plot
 p= ggplot(snp_density, aes(x = Window_start, y = SNP_count)) +
    geom_bar(stat = "identity", position = "stack", fill = "steelblue", colour = "darkblue") +  
    facet_grid(rows = vars(Chrom), scales = "free_y") +  # Correctly facet by rows
    theme_minimal(base_size = 14) +         # Increase base font size for better readability
    labs(x = "Genomic Location (bp)", 
         y = "SNP Density", 
         title = "SNP Density Plot by Chromosome") +
    theme(
      axis.text.x = element_blank(),          # Remove x-axis labels
      axis.ticks.x = element_blank(),         # Remove x-axis ticks
      strip.text.y = element_text(size = 12), # Customize facet label size for y-axis
      axis.text.y = element_blank(),           # Remove y-axis labels
      axis.ticks.y = element_blank()           # Remove y-axis ticks
    )
 print(p)
}

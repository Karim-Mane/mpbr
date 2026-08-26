#' Check whether a VCF file contains SNPs overlapping a set of genes of interest
#'
#' @param snpdata A <SNPdata> object obtained from `get_snpdata()` function
#' @param regions A character with the path to the `BED` file with the
#'    genomic coordinates of the regions of interest. It must contain at least
#'    the CHROM, START, and END columns
#'
#' @returns A list of two elements: the variants overlapping each region and the
#'    region's variant counts.
#' @export
#'
#' @examples
#' \dontrun{
#'   res <- get_region_variants(
#'     snpdata = snpdata,
#'     regions = file.path("data", "genes_cordinates_on_v68.txt")
#'   )
#'   region_variants <- res[["region_variants"]]
#'   variant_count <- res[["variant_count"]]
#' }
get_region_variants <- function(snpdata, regions) {
  checkmate::assert_class(snpdata, classes = "SNPdata", null.ok = FALSE)
  checkmate::assert_file_exists(regions, access = "r")
  
  # extract the SNPs genomic coordinates
  details <- snpdata[["details"]]
  genomic_coordinates <- details |>
    dplyr::select(Chrom, Pos)
  
  # read in the region file
  regions <- data.table::fread(regions)
  
  # prepare a data frame to find overlap between the genomic coordinates and the
  # region annotation
  pos_dt <- data.table::data.table(
    start = genomic_coordinates[["Pos"]],
    end = genomic_coordinates[["Pos"]]
  )
  
  # find the overlaps between the regions and the SNPs genomic coordinates
  data.table::setkey(regions, start, end)
  result <- data.table::foverlaps(
    pos_dt, regions, type = "within", nomatch = NULL
  )
  result <- result |>
    dplyr::select(chrom, i.start, gene_id, gene_name, gene_description)
  names(result)[[2]] <- "pos"
  data.table::setkeyv(result, c("chrom", "pos"))
  
  # join the genes to their annotations
  annotation <- details |>
    dplyr::select(Chrom, Pos, effect, impact, hgvs_c, hgvs_p)
  names(annotation)[1:2] <- c("chrom", "pos")
  annotation <- data.table::as.data.table(annotation)
  data.table::setkeyv(annotation, c("chrom", "pos"))
  result <- result |>
    dplyr::left_join(annotation, by = c("chrom", "pos"))
  
  # get the count
  counts <- data.table::setDT(result)[, .N, by = .(gene_id)]
  names(counts)[[2]] <- "num_variants"
  
  # reorder the genes by chromosome
  result[["chrom"]] <- factor(
    result[["chrom"]],
    levels = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3",
               "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3",
               "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3",
               "Pf3D7_13_v3", "Pf3D7_14_v3", "Pf3D7_MIT_v3")
  )
  result <- result |> dplyr::arrange(chrom)
  
  return(list(
    region_variants <- result,
    variant_count <- counts
  ))
}

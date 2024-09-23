
test_that("compute_maf works as expected", {
  skip_on_ci() 
  skip_on_covr()
  # build the SNPdata object
  snpdata <- get_snpdata(
    vcf_file   = system.file("extdata", "Input_Data.vcf.gz", package = "mpbr"), 
    meta_file  = system.file("extdata", "SampleMetadata.RDS", package = "mpbr"), 
    output_dir = tempdir(), 
    gof        = system.file("extdata", "pf_gene_ontology.RDS", package = "mpbr"), 
    gff        = system.file("extdata", "PlasmoDB-56_Pfalciparum3D7.RDS",
                             package = "mpbr"),
    num_threads = 2L
  )
  
  # test if compute_maf works fine when include_het = FALSE
  test_maf <- compute_maf(snpdata,
                          include_het = FALSE,
                          mat_name    = "GT")
  expect_true(inherits(test_maf, "SNPdata"))
  expect_identical(names(test_maf), c("meta", "details", "GT", "vcf"))
  expect_identical(test_maf[["meta"]], snpdata[["meta"]])
  expect_true(ncol(test_maf[["details"]]) > ncol(snpdata[["details"]]))
  expect_true(all(c("MAF", "MAF_allele") %in% names(test_maf[["details"]])))
  expect_identical(test_maf[["GT"]], snpdata[["GT"]])
  expect_identical(test_maf[["vcf"]], snpdata[["vcf"]])
  expect_identical(test_maf[["index"]], snpdata[["index"]])
  expect_identical(unique(test_maf[["details"]][["MAF_allele"]]),
                   c("1", "0", "0=1"))
  
  # test if compute_maf works fine when include_het = TRUE
  test_maf <- compute_maf(snpdata,
                          include_het = TRUE,
                          mat_name    = "GT")
  expect_true(inherits(test_maf, "SNPdata"))
  expect_identical(names(test_maf), c("meta", "details", "GT", "vcf"))
  expect_identical(test_maf[["meta"]], snpdata[["meta"]])
  expect_true(ncol(test_maf[["details"]]) > ncol(snpdata[["details"]]))
  expect_true(all(c("MAF", "MAF_allele") %in% names(test_maf[["details"]])))
  expect_identical(test_maf[["GT"]], snpdata[["GT"]])
  expect_identical(test_maf[["vcf"]], snpdata[["vcf"]])
  expect_identical(test_maf[["index"]], snpdata[["index"]])
  expect_identical(unique(test_maf[["details"]][["MAF_allele"]]),
                   c("1", "0", "0/1", "0=1=2"))
})

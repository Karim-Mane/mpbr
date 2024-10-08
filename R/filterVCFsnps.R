calculate_missingness = function(gt_data){
  missingness = apply(gt_data, 1, function(x){
    length(x[is.na(x)])/length(x)
  })
}

read_vcf = function(vcf, path){
  
  samples      <- file.path(path, "Samples.txt")
  header       <- file.path(path, "Header.txt")
  body         <- file.path(path, "Body.txt")
  ## Read VCF file
  os_type <- Sys.info()["sysname"]
  
  if (os_type == "Windows") {
    # Use gzip for Windows systems
    # get list of sample ids in vcf
    system(sprintf("gzip -dc %s | findstr /R \"^#CHROM\" | cut -d\"\\t\" -f 10- | tr \"\\t\" \"\\n\" > %s", vcf_file, samples))
    
    # get header of vcf
    system(sprintf("gzip -dc %s | findstr /V \"^#\" > %s", vcf_file, body))
    
    # get the body of the vcf
    system(sprintf("gzip -dc %s | findstr /R \"^#\" > %s", vcf_file, header))
  } else {
    # Use gzcat for Unix-like systems (macOS, Linux)
    # get list of sample ids in vcf
    system(sprintf("gzcat %s | grep '^#CHROM' | cut -d$'\t' -f 10- | paste -sd '\t' - | tr '\t' '\n' > %s", vcf_file, samples))
    
    #get header of vcf
    system(sprintf("gzcat %s | grep -v '^#' > %s", vcf_file, body))
    
    # get the body of the vcf
    system(sprintf("gzcat %s | grep '^#' > %s", vcf_file, header))
    
  }
  return(list(samples=samples, header=header, body=body))
}


#Remove snps from vcf
remove_snps_from_vcf <- function(vcf, loci_to_be_retained, path, index = 1L) {
  checkmate::assert_file_exists(vcf)
  #checkmate::assert_character(loci_to_be_retained, any.missing = FALSE,
  #                           null.ok = FALSE)
  checkmate::assert_directory_exists(path)
  checkmate::assert_numeric(index, lower = 0L, len = 1L, null.ok = FALSE,
                            any.missing = FALSE)
  #target_loci  <- file.path(path, loci_to_be_retained)
  
  #correct_rows <- file.path(path, "Good_snps.txt")
  filtered_vcf <- file.path(path, paste0("Filtered_snps_", index, ".vcf"))
  bodyFilter   <- file.path(path, "BodyFilter.txt")
  loci_to_be_retained = file.path(path, "details.txt")
  
  vcfData = read_vcf(vcf_file, path)
  
  #retain snp in list
  if (os_type == "Windows") {
    # Extract common rows based on columns 1 and 2
    system(sprintf("powershell -Command \"Get-Content '%s' | ForEach-Object { $line = $_ -split '\\t'; $key = $line[0] + $line[1]; $key } | Sort-Object | Get-Unique | Out-File '%s'\"", 
                   loci_to_be_retained, bodyFilter))
    
    # Filter vcf_data_file based on keys from bodyFilter
    system(sprintf("powershell -Command \"Get-Content '%s' | ForEach-Object { $line = $_ -split '\\t'; $key = $line[0] + $line[1]; if (Test-Path '%s' -and (Select-String -Pattern $key -Path '%s')) { $_ } } | Out-File '%s'\"", 
                   vcf_data_file, bodyFilter, bodyFilter, filtered_vcf))
    
    # Concatenate header and bodyFilter into filtered_vcf
    system(sprintf("powershell -Command \"Get-Content '%s', '%s' | Set-Content '%s'\"", 
                   header, bodyFilter, filtered_vcf))
    
    # Compress the filtered VCF file using bgzip (make sure bgzip is accessible)
    system(sprintf("powershell -Command \"& 'path\\to\\bgzip.exe' '%s'\"", 
                   filtered_vcf))
    
    # Remove temporary files using PowerShell
    system(sprintf("powershell -Command \"Remove-Item '%s', '%s', '%s' -Force\"", 
                   header, bodyFilter, correct_rows))
    
    # Define the final output path for the compressed VCF
    filtered_vcf_compressed <- file.path(path, paste0("Filtered_snps_", index, ".vcf.gz"))
    
  }
  else{
    
    #system(sprintf("awk 'NR==FNR {a[$1]; next} $2 in a' %s %s > %s", vcf_data[[body]], loci_to_be_retained, bodyFilter))
    #Extract common rows based on columns 1 and 2 of the vcf body and loci_to_be_retained
    system(sprintf("awk -F'\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' %s %s > %s",
                   loci_to_be_retained, vcfData[["body"]], bodyFilter))
    #concatenate body and header to form vcf
    system(sprintf("cat %s %s > %s", header, bodyFilter, filtered_vcf))
    
    #zip vcf file
    system(sprintf("bgzip %s", filtered_vcf))
    
    #remove pre produced files to clear hard disk
    system(sprintf("rm -f %s %s %s", header, body, correct_rows))
    
    #get path to zipped vcf
    filtered_vcf <- file.path(path, paste0("Filtered_snps_", index, ".vcf.gz"))
  }
  #return vcf path
  as.character(filtered_vcf)
}

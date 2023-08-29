#' Calculate Weir & Cockerham's Fst
#' 
#' @param snpdata `SNPdata` object
#' @param groups a vector of population names. Every sample in the metadata file
#'    is associated to its population of origin. The differentiation index will
#'    be estimated between these groups.
#' @param from the metadata column that contains the sample's population
#'    information.
#'    
#' @return a `SNPdata` object with an extra field named as **Fst**. This is a
#'    `list` of data frames with the Fst values for each pair of comparison.
#'    Every data frame contains N rows (where N is the number of loci) and the
#'    following 7 columns:
#'    \enumerate{
#'      \item the chromosome ID
#'      \item the SNPs positions
#'      \item the allele frequency in the first group
#'      \item the allele frequency in the second group
#'      \item the resulting Fst values
#'      \item the p-values associated with the Fst results
#'      \item the p-values corrected for multiple testing using the 
#'            Benjamini-Hochberg method
#'    }
#'    
#' @examples
#' \dontrun{
#'   snpdata <- calculate_wcFst(
#'     snpdata,
#'     groups = c("Senegal", "Gambia"),
#'     from   = "Country"
#'   )
#' }
#' @export
#' 
calculate_wcFst <- function(snpdata, from, groups) {
  checkmate::assert_vector(groups, any.missing = FALSE, min.len = 0,
                           null.ok = TRUE)
  checkmate::assert_character(from, len = 1, any.missing = FALSE,
                              null.ok = FALSE)
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  
  if (is.null(groups) & !is.null(from)) {
    groups <- as.character(unique(snpdata$meta[[from]]))
  } else if(all(!is.null(groups) &&
                !is.null(from) &&
                (any(!(groups %in% unique(snpdata$meta[[from]])))))) {
        stop(sprintf("Not all specified groups belong to the %s column of 
             the metadata table", from))
  }
  system(sprintf("tabix %s", snpdata$vcf))
  Fst <- list()
  for(i in 1:(length(groups)-1)) {
      idx1 <- which(snpdata$meta[[from]] == groups[i])
      idx1 <- paste(idx1, collapse = ",")
      for(j in (i+1):length(groups)) {
          idx2 <- which(snpdata$meta[[from]] == groups[j])
          idx2 <- paste(idx2, collapse = ",")
          out  <- file.path(dirname(snpdata$vcf), "out.wc.fst")
          pout <- file.path(dirname(snpdata$vcf), "out.wc.fst.pvalues")
          system(sprintf("wcFst --target %s --background %s --file %s \
                         --deltaaf 0 --type GT > %s",
                         idx1, idx2, snpdata$vcf, out))
          system(sprintf("pFst --target %s --background %s --file %s \
                         --deltaaf 0 --type GT > %s",
                         idx1, idx2, snpdata$vcf, pout))
          out        <- data.table::fread(out)
          pout       <- data.table::fread(pout)
          data.table::setkeyv(out,  c("V1","V2"))
          data.table::setkeyv(pout, c("V1","V2"))
          tmp        <- out[pout, nomatch = 0]
          names(tmp) <- c("Chrom", "Pos", "AF_in_target", "AF_in_background",
                          "wcFst", "wcFst_pvalue")
          idx        <- which(tmp$wcFst < 0)
          if (length(idx) > 0) {
              tmp$wcFst[idx] <- 0
          }
          tmp$wcFst_Adj_pvalue_BH <- p.adjust(tmp$wcFst_pvalue, method = "BH")
          Fst[[paste0(groups[i], "_vs_", groups[j])]] <- tmp
      }
  }
  snpdata$Fst <- Fst
  snpdata
}

#' Calculate LD R^2 between pairs of loci
#' 
#' @param snpdata a `SNPdata` object
#' @param min_r2 the minimum r2 value below which the LD value is not reported
#' @param inter_chrom a `logical` that determines whether to calculate the
#'    inter-chromosomal LD or not. Default is `FALSE`.
#' @param chroms a character vector of chromosome names. If provided, LD will be
#'    calculated only across these chromosomes.
#' 
#' @return a `SNPdata` object with an extra field named as **LD**
#' @examples
#' \dontrun{
#'   snpdata <- calculate_LD(snpdata, min_r2=0.2, inter_chrom=FALSE,
#'     chroms=c("Pf3D7_04_v3","Pf3D7_05_v3"))
#'  }
#' @details The output file from LD calculation could be large. In order to
#'     reduce the size of that file, we recommend to specify the list of
#'     chromosomes for which LD should be calculated using the `chroms` option.
#' @export
calculate_LD <- function(snpdata, min_r2 = 0.2, inter_chrom = FALSE,
                         chroms = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_r2, lower = 0, upper = 1, finite = TRUE,
                            any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_logical(inter_chrom, any.missing = FALSE, len = 1,
                            null.ok = FALSE)
    out = paste0(dirname(snpdata$vcf),"/tmp_ld")
    if(inter_chrom){
        cat("inter-chromosomal LD will be calculated between sites on",paste(chroms, collapse = ","))
        system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --interchrom-geno-r2", snpdata$vcf, out, min_r2))
    }else{
        system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --geno-r2", snpdata$vcf, out, min_r2))
    }
    out = paste0(out,".geno.ld")
    system(sprintf("bgzip %s", out))
    out = paste0(out,".gz")
    ld = fread(out, nThread = 4)
    if(!is.null(chroms)){
        idx1 = which(ld$CHR1 %in% chroms)
        idx2 = which(ld$CHR2 %in% chroms)
        idx = unique(c(idx1, idx2))
        ld = ld[idx,]
    }

    snpdata$LD=ld
}

#' Generate dissimilarity matrix (1-IBS) between all pairs of isolates
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="GT"
#' @return SNPdata object with an extra field: IBS
#' @examples
#' \dontrun{
#'   snpdata = calculate_IBS(snpdata, mat.name="GT")
#'  }
#'  
#' @export
calculate_IBS = function(snpdata, mat.name="GT"){
    if(!(mat.name %in% names(snpdata))){
        stop("specified genotype matrix does not exist")
    }
    X = t(snpdata[[mat.name]])
    y = matrix(NA, nrow = nrow(X), ncol = nrow(X))
    # colnames(y) = rownames(y) = rownames(X)
    pb = txtProgressBar(min = 0, max = nrow(X), initial = 0,style = 3, char = "*")
    for(i in 1:nrow(X)){
        for(j in 1:nrow(X)){
            m = X[i,]-X[j,]
            y[i,j] = 1-(length(which(m==0))/ncol(X))
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    y[upper.tri(y)]=NA
    colnames(y) = rownames(X); rownames(y) = rownames(X)
    snpdata$IBS=y
    snpdata
}


#' Calculate iR index to detect loci with excess of IBD sharing
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="Phased"
#' @param family the name of the column, in the metadata table, to be used to represent the sample's population
#' @param number.cores the number of cores to be used. default=4
#' @return SNPdata object with an extra field: iR
#' @examples
#' \dontrun{
#'   snpdata = calculate_iR(snpdata, mat.name="Phased", family="Location", 
#'   number.cores=4)
#'  }
#' @export
calculate_iR = function(snpdata, mat.name="Phased", family="Location", number.cores=4){
    if(!(family %in% names(snpdata$meta))){
        stop("No column name ",family," in the metadata table")
    }
    if(mat.name=="GT" | mat.name=="Phased"){
        cat("phasing the mixed genotypes\n")
        snpdata = phase_mixed_genotypes(snpdata, nsim=10)
    }
    ped = make_ped(snpdata[[mat.name]], snpdata$meta, family)
    ped$sex=as.numeric(ped$sex)
    map = make_map(snpdata$details)
    ped.map = list(ped,map)
    my.geno = getGenotypes(ped.map, reference.ped.map=NULL, maf=0.01, isolate.max.missing=0.2, snp.max.missing=0.2, chromosomes=NULL, input.map.distance="cM", reference.map.distance="cM")
    my.param = getIBDparameters(ped.genotypes = my.geno, number.cores = number.cores)
    my.ibd = getIBDsegments(ped.genotypes = my.geno,parameters = my.param, number.cores = number.cores, minimum.snps = 20, minimum.length.bp = 50000,error = 0.001)
    my.matrix = getIBDmatrix(ped.genotypes = my.geno, ibd.segments = my.ibd)
    my.iR = getIBDiR(ped.genotypes = my.geno,ibd.matrix = my.matrix,groups = NULL)
    pvalues = as.numeric(as.character(lapply(my.iR$log10_pvalue, get.pvalue)))
    my.iR$log10_pvalue = pvalues
    names(my.iR)[8] = "pvalues"
    my.iR$adj_pvalue_BH = p.adjust(my.iR$pvalues, method = "BH")
    if(!("iR" %in% names(snpdata))){
        snpdata$iR=list()
    }
    groups = unique(snpdata$meta[[family]])
    snpdata$iR[[paste0(groups[1],"_vs_",groups[2])]] = my.iR
    snpdata
}

get.pvalue = function(x){10^-x}

make_ped = function(mat, metadata, family){
    mat[mat==1]=2
    mat[mat==0]=1
    mat[is.na(mat)]=0
    mat=t(mat)

    new.genotype = matrix(-9 , nrow = dim(mat)[1], ncol = ((dim(mat)[2])*2) )
    j=1
    k=1
    while(j<=ncol(mat)){
        # print(j)
        new.genotype[,k] = mat[,j]
        new.genotype[,(k+1)] = mat[,j]
        j=j+1
        k=k+2
    }

    PED6 = data.frame(cbind(metadata[[family]],metadata$sample, fatherID = 0, motherID = 0, sex = metadata$COI, phenotype = -9), stringsAsFactors = FALSE)
    names(PED6)[1:2] = c("pop","sample")
    ped = data.frame(cbind(PED6, data.frame(new.genotype)))
    ped
}

make_map = function(details){
    c1 = details$Chrom
    c4 = details$Pos
    get_xme = function(x){as.character(unlist(strsplit(x,"_"))[2])}
    xme = as.numeric(as.character(lapply(c1, get_xme)))
    c2 = paste0(LETTERS[xme], c4)
    c3 = c4/12000
    map = data.frame(chrom = c1, post = c2, gd = c3, pos = c4, stringsAsFactors = FALSE)
    map
}

haplotypeToGenotype = function(haplotypes, moi){
    genotypes = matrix(-9, dim(haplotypes)[1], dim(haplotypes)[2]/2)
    for(i in 1:dim(haplotypes)[1]){
        #print(i)
        j=1
        k=1
        while(j <= dim(haplotypes)[2]){
            if(haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 1) genotypes[i,k] = 0
            if(haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 2) genotypes[i,k] = 2
            if(haplotypes[i,j] == 0 && haplotypes[i,(j+1)] == 0) genotypes[i,k] = -1
            if((moi[i] == 1) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = -1
            if((moi[i] == 2) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = 1
            if(haplotypes[i,j] != 0 && haplotypes[i,j] != 1 && haplotypes[i,j] != 2) genotypes[i,k] = -1
            if(haplotypes[i,(j+1)] != 0 && haplotypes[i,(j+1)] != 1 && haplotypes[i,(j+1)] != 2) genotypes[i,k] = -1
            j=j+2
            k=k+1
        }
    }
    return(t(genotypes))
}

calculatePopAlleleFreq = function(genotypes, moi){
    pop_allele_freqs = vector('numeric',dim(genotypes)[1])
    number_isolates = dim(genotypes)[2]
    number_snps = dim(genotypes)[1]

    for (t in 1:number_snps){
        #print(t)
        A = 0
        B = 0
        for (i in 1:number_isolates)
        {
            if (genotypes[t,i] == 0 && moi[i] == 1)  A=A+1
            if (genotypes[t,i] == 0 && moi[i] == 2)  A=A+2
            if (genotypes[t,i] == 1 && moi[i] == 2){ A=A+1; B=B+1 }
            if (genotypes[t,i] == 2 && moi[i] == 1)  B=B+1
            if (genotypes[t,i] == 2 && moi[i] == 2)  B=B+1
        }
        if (A + B == 0) pop_allele_freqs[t] = -1
        else pop_allele_freqs[t] = A/(A+B)
    }
    return(pop_allele_freqs)
}

calculateMissingness = function(genotypes) {
    proportion_missing=vector('numeric',dim(genotypes)[2])
    number_snps = dim(genotypes)[2]
    number_isolates = dim(genotypes)[1]
    number_snps_1 = dim(genotypes)[1]

    for (i in 1:number_snps){
        #print(i)
        number_missing = 0.0
        for (j in 1:number_isolates) {
            if(genotypes[j,i] == -1 ) number_missing = number_missing+1;
        }
        proportion_missing[i] = number_missing/number_snps_1;
    }
    return(proportion_missing)
}

getGenotypes = function(ped.map, reference.ped.map=NULL, maf=0.01, isolate.max.missing=0.1, snp.max.missing=0.1, chromosomes=NULL, input.map.distance="cM", reference.map.distance="cM"){
    # check input PED and MAP files
    if (!is.list(ped.map) | length(ped.map) != 2) stop ("'ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
    input.ped <- as.data.frame(ped.map[[1]])
    input.map <- as.data.frame(ped.map[[2]])
    if (!is.data.frame(input.ped)) stop ("'ped.map' has incorrect format - PED is not a data.frame")
    if (!is.data.frame(input.map)) stop ("'ped.map' has incorrect format - MAP is not a data.frame")

    # check the PED and MAP files have the same number of SNPs
    if (ncol(input.ped) != (2*nrow(input.map)+6)) stop ("'ped.map' has incorrect format - PED and MAP must have the same number of SNPs")

    # check the MAP file has 4 coloumns
    if (ncol(input.map) != 4) stop ("'ped.map' has incorrect format - MAP must have four columns")
    colnames(input.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
    if (!is.numeric(input.map[,"pos_M"])) stop ("'ped.map' has incorrect format - genetic-map positions in MAP file are non-numeric")
    if (!is.numeric(input.map[,"pos_bp"])) stop ("'ped.map' has incorrect format - base-pair positions in MAP file are non-numeric")
    if (is.factor(input.map[,"chr"]))
        input.map[,"chr"] <- as.character(input.map[,"chr"])
    if (is.factor(input.map[,"snp_id"]))
        input.map[,"snp_id"] <- as.character(input.map[,"snp_id"])


    # check PED for factors
    # input.ped = as.matrix(input.ped)
    # input.map = as.matrix(input.map)
    if (is.factor(input.ped[,1]))
        input.ped[,1] <- as.character(input.ped[,1])
    if (is.factor(input.ped[,2]))
        input.ped[,2] <- as.character(input.ped[,2])


    # check reference data
    if (!is.null(reference.ped.map)) {
        if (!is.list(reference.ped.map) | length(reference.ped.map) != 2) stop ("'reference.ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
        reference.ped <- reference.ped.map[[1]]
        reference.map <- reference.ped.map[[2]]
        if (!is.data.frame(reference.ped)) stop ("'reference.ped.map' has incorrect format - PED is not a data.frame")
        if (!is.data.frame(reference.map)) stop ("'reference.ped.map' has incorrect format - MAP is not a data.frame")

        # check the PED and MAP files have the same number of SNPs
        if (ncol(reference.ped) != (2*nrow(reference.map)+6)) stop ("'reference.ped.map' has incorrect format - PED and MAP must have the same number of SNPs")

        # check the MAP file has 4 coloumns
        if (ncol(reference.map) != 4) stop ("'reference.ped.map' has incorrect format - MAP must have four columns")
        colnames(reference.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
        if (!is.numeric(reference.map[,"pos_M"])) stop ("'reference.ped.map' has incorrect format - genetic-map positions in reference MAP file are non-numeric")
        if (!is.numeric(reference.map[,"pos_bp"])) stop ("'reference.ped.map' has incorrect format - base-pair positions in reference MAP file are non-numeric")
        if (is.factor(reference.map[,"chr"]))
            reference.map[,"chr"] <- as.character(reference.map[,"chr"])
        if (is.factor(reference.map[,"snp_id"]))
            reference.map[,"snp_id"] <- as.character(reference.map[,"snp_id"])

        # check PED for factors
        if (is.factor(reference.ped[,1]))
            reference.ped[,1] <- as.character(reference.ped[,1])
        if (is.factor(reference.ped[,2]))
            reference.ped[,2] <- as.character(reference.ped[,2])
    }

    # check maf
    if (!is.vector(maf)) stop ("'maf' has incorrect format - must be a vector")
    if (!is.numeric(maf)) stop ("'maf' has incorrect format - must be numeric")
    if (length(maf) != 1) stop ("'maf' has incorrect format - must be a single numeric value")
    if (maf > 1 | maf < 0) stop ("'maf' has incorrect format - must be between 0 and 1")

    # check isolate.max.missing
    if (!is.vector(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be numeric")
    if (length(isolate.max.missing) != 1) stop ("'isolate.max.missing' has incorrect format - must be a single numeric value")
    if (isolate.max.missing > 1 | isolate.max.missing < 0) stop ("'isolate.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")

    # check snp.max.missing
    if (!is.vector(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be numeric")
    if (length(snp.max.missing) != 1) stop ("'snp.max.missing' has incorrect format - must be a single numeric value")
    if (snp.max.missing > 1 | snp.max.missing < 0) stop ("'snp.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")

    # check chromosomes
    if (!is.null(chromosomes)) {
        if (!is.vector(chromosomes)) stop ("'chromosomes' has incorrect format - must be a vector")
        if(!all(chromosomes %in% input.map[,"chr"]))
            stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% input.map[,"chr"])]," not in 'ped.map'\n")))
        if (!is.null(reference.ped.map)) {
            if(!all(chromosomes %in% reference.map[,"chr"]))
                stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% reference.map[,"chr"])]," not in 'reference.ped.map'\n")))
        }
    } else
        chromosomes <- unique(as.character(input.map[,"chr"]))

    # check input map distance
    if (!is.vector(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be a vector")
    if (!is.character(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (length(input.map.distance) != 1) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (input.map.distance != "M" & input.map.distance != "cM")
        stop ("'input.map.distance' has incorrect format - must be either M or cM")
    if (input.map.distance == "cM") {
        input.map[,"pos_M"] <- input.map[,"pos_M"]/100
    }

    # check reference map distance
    if (!is.vector(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be a vector")
    if (!is.character(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (length(reference.map.distance) != 1) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (reference.map.distance != "M" & reference.map.distance != "cM")
        stop ("'reference.map.distance' has incorrect format - must be either M or cM")
    if (reference.map.distance == "cM" & !is.null(reference.ped.map)) {
        reference.map[,"pos_M"] <- reference.map[,"pos_M"]/100
    }

    # begin data filtering

    # create new isolate IDs from PED FIDs and IIDs
    isolate.names <- paste(input.ped[,1], input.ped[,2], sep="/")
    if (any(duplicated(isolate.names)))
        stop ("duplicate sample IDs found")

    # merge input data with reference data
    if (!is.null(reference.ped.map)) {
        input.map.v1      <- cbind(1:nrow(input.map), input.map)
        reference.map.v1  <- cbind(1:nrow(reference.map), reference.map)
        input.map.v1      <- merge(input.map.v1, reference.map.v1, by.x="snp_id", by.y="snp_id")
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map'")
        input.map.v1  <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
        if (!is.null(chromosomes))
            input.map.v1 <- input.map.v1[input.map.v1[,"chr.x"] %in% chromosomes,]
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map' for selected chromosomes")
        input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
        input.ped.columns <- input.ped.columns[order(input.ped.columns)]
        input.ped.v1      <- input.ped[,input.ped.columns]
        reference.ped.columns  <- c(1:6, 2*input.map.v1[,"1:nrow(reference.map)"] + 5, 2*input.map.v1[,"1:nrow(reference.map)"] + 6)
        reference.ped.columns  <- reference.ped.columns[order(reference.ped.columns)]
        reference.ped.v1       <- reference.ped[,reference.ped.columns]
        input.map.v2           <- input.map.v1[,c("chr.x", "snp_id", "pos_M.x", "pos_bp.x")]
        colnames(input.map.v2) <- c("chr", "snp_id", "pos_M","pos_bp")
    } else {
        if (!is.null(chromosomes)) {
            input.map.v1 <- cbind(1:nrow(input.map), input.map)
            input.map.v1 <- input.map.v1[input.map.v1[,"chr"] %in% chromosomes,]
            if (nrow(input.map.v1) == 0)
                stop ("no SNPs remaining after subsetting 'ped.map'by selected chromosomes")
            input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
            input.ped.columns <- input.ped.columns[order(input.ped.columns)]
            input.ped.v1      <- input.ped[,input.ped.columns]
            input.map.v2      <- input.map.v1[,c("chr", "snp_id", "pos_M", "pos_bp")]
        } else {
            input.map.v2 <- input.map
            input.ped.v1 <- input.ped
        }
    }

    # call genotypes
    input.matrix        <- as.matrix(input.ped.v1[,7:ncol(input.ped.v1)])
    input.genders       <- input.ped.v1[,5]
    input.genotypes.v0  <- cbind(input.map.v2, haplotypeToGenotype(input.matrix, input.genders))
    if (!is.null(reference.ped.map)) {
        reference.matrix       <- as.matrix(reference.ped.v1[,7:ncol(reference.ped.v1)])
        reference.genders      <- reference.ped.v1[,5]
        reference.genotypes.v0 <- cbind(input.map.v2, haplotypeToGenotype(reference.matrix, reference.genders))
    }


    # calculate allele frequencies form reference data
    if (is.null(reference.ped.map)) {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(input.genotypes.v0[,5:ncol(input.genotypes.v0)]), input.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,5:ncol(input.genotypes.v0)])
    } else {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(reference.genotypes.v0[,5:ncol(reference.genotypes.v0)]), reference.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,c(5:ncol(input.genotypes.v0))])
    }
    colnames(input.genotypes.v1) <- c("chr", "snp_id", "pos_M","pos_bp", "freq", isolate.names)
    cat(paste("Begin filtering of",length(isolate.names),"isolates and ",nrow(input.genotypes.v1),"SNPs...\n",sep=""))


    # remove SNPs with low population MAF
    input.genotypes.v2 <- subset(input.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf) #removing SNPs with AF>0.99 and AF<0.01
    if (nrow(input.genotypes.v2) == 0)
        stop("0 SNPs remain after MAF removal")
    cat(paste(nrow(input.genotypes.v2),"SNPs remain after MAF removal...\n",sep=""))


    # remove snps with high missingness
    snp.missingness    <- calculateMissingness(as.matrix(t(input.genotypes.v2[,6:ncol(input.genotypes.v2)])))
    input.genotypes.v3 <- input.genotypes.v2[snp.missingness <= snp.max.missing,]
    if (nrow(input.genotypes.v3) == 0)
        stop("0 SNPs remain after missingness removal")
    cat(paste(nrow(input.genotypes.v3),"SNPs remain after missingness removal...\n",sep=""))


    # remove samples with high missingness
    isolate.missingness <- round(calculateMissingness(as.matrix(input.genotypes.v3[,6:ncol(input.genotypes.v3)])),digits=3)
    if (length(isolate.names[isolate.missingness > isolate.max.missing]) > 0) {
        my.remove <- isolate.names[isolate.missingness > isolate.max.missing]
        warning("isolates removed due to genotype missingness: ",paste(my.remove, collapse=", "))
        sample.keep        <- input.ped.v1[isolate.missingness <= isolate.max.missing,1:6]
        input.genotypes.v4 <- input.genotypes.v3[,c(1:5, which(isolate.missingness <= isolate.max.missing) + 5)]
        if(nrow(sample.keep) < 1) stop(paste("All isolates removed with missingness > ",isolate.max.missing*100,"%. No isolates remaining.",sep=""))
    } else {
        sample.keep        <- input.ped.v1[,1:6]
        input.genotypes.v4 <- input.genotypes.v3
    }
    colnames(sample.keep) <- c("fid", "iid", "pid", "mid", "moi", "aff")
    if ((ncol(input.genotypes.v4)-5) == 0) {
        stop("0 samples remain after missingness removal")
    }
    cat(paste(ncol(input.genotypes.v4)-5,"isolates remain after missingness removal...\n",sep=""))


    return.genotypes <- list(sample.keep, input.genotypes.v4)
    names(return.genotypes) <- c("pedigree", "genotypes")
    return(return.genotypes)

}

#' Estimate the relatedness between pairs of isolates
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="GT"
#' @param from the name of the column, in the metadata table, to be used to represent the sample's population
#' @param sweepRegions a data frame with the genomic coordinates of the regions of the genome to be discarded. This should contain the following 3 columns:
#' \enumerate{
#' \item Chrom: the chromosome ID
#' \item Start: the start position of the region on the chromosome
#' \item End: the end position of the region on the chromosome
#' }
#' @param groups a vector of character. If specified, relatedness will be generated between these groups
#' @return SNPdata object with an extra field: relatedness. This will contain the relatedness data frame of 3 columns and its correspondent matrix
#' @examples
#' \dontrun{
#'   calculate_relatedness(snpdata, mat.name="GT", family="Location", 
#'   sweepRegions=NULL, groups=c("Chogen","DongoroBa"))
#'  }
#' @details The relatedness calculation is based on the model developed by Aimee R. Taylor and co-authors. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009101
#' @export
calculate_relatedness = function(snpdata, mat.name="Imputed", from="Location", sweepRegions=NULL, groups=NULL){
    # sourceCpp("src/hmmloglikelihood.cpp")
    details = snpdata$details
    metadata = snpdata$meta
    if(mat.name=="GT" | mat.name=="Phased"){
        cat("Imputing the missing genotypes\n")
        snpdata = impute_missing_genotypes(snpdata, genotype=mat.name, nsim=10)
    }
    if((mat.name=="Imputed") & (!("Imputed" %in% names(snpdata)))){
        mat.name="GT"
        cat("Imputing the missing genotypes\n")
        snpdata = impute_missing_genotypes(snpdata, genotype=mat.name, nsim=10)
    }

    # mat.name = "Imputed"
    mat = snpdata[["Imputed"]]
    if(!is.null(groups) & all(groups %in% unique(metadata[[from]]))){
        pops = groups
    }else{
        pops = unique(metadata[[from]])
    }
    cat("creating the genotype data\n")
    genotypes = constructGenotypeFile(pops, details, mat, metadata, from, sweepRegions)
    sites = names(genotypes)
    dir = dirname(snpdata$vcf)
    ibd = NULL
    for(ii in 1:length(sites)){
        site1 = sites[ii]
        for(jj in ii:length(sites)){
            site2 = sites[jj]
            cat("calculating the relatedness between",site1,"and",site2,"\n")
            ibd = rbind(ibd, gen_mles(genotypes, site1, site2, f=0.3, dir))
        }
    }
    ibd = as.data.table(ibd)
    names(ibd) = c('iid1','iid2','k','relatedness')
    ibd = subset(ibd, select = -3)
    ibd$relatedness = as.numeric(ibd$relatedness)
    cat("creating the relatedness matrix\n")
    rmatrix = createRelatednessMatrix(ibd, metadata, pops, from)
    snpdata$relatedness = list()
    snpdata$relatedness[["df"]] = ibd
    snpdata$relatedness[["matrix"]] = rmatrix
    system(sprintf("rm -rf %s", paste(dir,"/ibd")))
    snpdata
}

createRelatednessMatrix = function(ibd, metadata, pops, from){
    idx = NULL
    for(pop in pops){
        idx = c(idx, which(metadata[[from]]==pop))
    }
    idx = unique(idx)
    metadata = metadata[idx,]
    modifier = function(x){gsub('-','.',x)}
    metadata$sample = as.character(lapply(metadata$sample,modifier))
    relatednessMatrix = matrix(NA,nrow = length(metadata$sample), ncol = length(metadata$sample))
    rownames(relatednessMatrix) = metadata$sample
    colnames(relatednessMatrix) = metadata$sample
    surLigne = unique(ibd$iid1)
    for(i in 1:length(surLigne)){
        # print(paste0('i=',i))
        ligne = surLigne[i]
        l = match(ligne,rownames(relatednessMatrix))
        target = ibd[which(ibd$iid1==ligne),]
        t = unique(target$iid2)
        for(j in 1:length(t)){
            tt = target[which(target$iid2==t[j]),]
            k = match(t[j],colnames(relatednessMatrix))
            if(nrow(tt) == 0)
                relatednessMatrix[l,k] = 0
            else if(nrow(tt)==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$relatedness), digits = 5)
            else if(nrow(tt)>1 & length(unique(tt$iid1))==1 & length(unique(tt$iid2))==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$relatedness[1]), digits = 5)
            else if(nrow(tt)>1 & length(unique(tt$iid1))>1 | length(unique(tt$iid2))>1)
                relatednessMatrix[l,k] = mean(round(as.numeric(tt$relatedness), digits = 5), na.rm = TRUE)
        }
    }
    relatednessMatrix
}

constructGenotypeFile = function(pops, details, mat, metadata, from, sweepRegions){
    if(!is.null(sweepRegions)){
        selectiveRegions = fread(sweepRegions)
        rtd = NULL
        for(j in 1:nrow(selectiveRegions))
            rtd = c(rtd, which(details$Chrom==selectiveRegions$Chrom[j] & (details$Pos>=selectiveRegions$Start[j] & details$Pos<=selectiveRegions$End[j])))
        rtd = unique(rtd)
        details = details[-rtd,]
        mat = mat[-rtd,]
    }
    res = list()
    # pops = unique(metadata[[from]])
    for(pop in pops){
        idx = which(metadata[[from]]==pop)
        s = metadata$sample[idx]
        m = match(s, colnames(mat))
        X=mat[,m]
        chroms=details$Chrom; pos=details$Pos; samps=metadata$sample[idx]
        L = list(X, chroms, pos, samps)
        names(L)=c('X','chroms','pos','samps')
        res[[pop]] = L
    }
    res
}

gen_mles = function(res, site1, site2, f=0.3, dir){
    nproc=700
    epsilon=0.001
    nboot=100
    Ps=c(0.025,0.975)
    outFilePrefix = paste0(site1,'_',site2)
    countryA = res[[site1]]
    countryB = res[[site2]]
    outputDir = paste0(dir,'/ibd')
    system(sprintf("mkdir -p %s", outputDir))
    if (site1==site2){
        # within country comparison
        L=run_country(countryA,f)
        data_set=L$data_set
        individual_names=L$individual_names
        nindividuals=L$nindividuals
        chrom=L$chrom
        pos=L$pos
        name_combinations = matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
        Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
        k=0
        for ( i in 1 : (nindividuals-1)){
            j=(1+k):(k+nindividuals-i)
            name_combinations[j,1]=rep(individual_names[i],each=length(j))
            name_combinations[j,2]=individual_names[(i+1):nindividuals]
            k=k+length(j)      # within country comparison
        }
    }else{
        # within country comparison
        LA=run_country(countryA,f)
        individual_names_A=LA$individual_names
        nindividuals_A=LA$nindividuals
        # # within country comparison
        LB=run_country(countryB,f)
        individual_names_B=LB$individual_names
        nindividuals_B=LB$nindividuals
        c_offset = dim(LA$data_set)[2]
        L=run_2country(countryA,countryB,f)
        data_set=L$data_set
        individual_names=L$individual_names
        nindividuals=L$nindividuals
        chrom=L$chrom
        pos=L$pos
        name_combinations <- matrix(nrow = nindividuals_A*nindividuals_B, ncol = 2)
        Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
        name_combinations[,1]=rep(individual_names_A,each=nindividuals_B)
        name_combinations[,2]=rep(individual_names_B,nindividuals_A)
    }
    X=as.matrix(data_set)
    data_set$fs = rowMeans(data_set, na.rm = TRUE) # Calculate frequencies
    data_set$pos =pos
    data_set$chrom=chrom
    data_set$dt <- c(diff(data_set$pos), Inf)
    pos_change_chrom <- 1 + which(diff(data_set$chrom) != 0) # find places where chromosome changes
    data_set$dt[pos_change_chrom-1] <- Inf
    # note NA result is undocumented - could change
    a0=Rfast::rowCountValues(X, rep(0,dim(X)[1]))  #getting the count of 0 on each row (SNPs)
    a1=Rfast::rowCountValues(X, rep(1,dim(X)[1]))  #getting the count of 1 on each row (SNPs)
    a2=Rfast::rowCountValues(X, rep(2,dim(X)[1]))  #getting the count of 2 on each row (SNPs)
    ana=rowSums(is.na(X))
    frequencies=cbind(a0,a1,a2)/(dim(X)[2]-ana)
    if (all.equal(rowSums(frequencies),rep(1,dim(X)[1]))!=TRUE){
        cat(paste0("frequency ERROR"))
        return()
    }
    # if (iloop < nproc){
    #     N=floor(dim(name_combinations)[1]/nproc)
    #     starty=(iloop-1)*N+1
    #     endy=iloop*N
    # }else{
    #     N=dim(name_combinations)[1]-(nproc-1)*floor(dim(name_combinations)[1]/nproc)
    #     starty = 1+(nproc-1)*floor(dim(name_combinations)[1]/nproc)
    #     endy=dim(name_combinations)[1]
    # }
    starty = 1
    endy = dim(name_combinations)[1]
    X=matrix(data=NA,nrow =endy-starty+1,ncol = 4 )
    for (icombination in starty:endy){
        #cat(paste0("icombination=",icombination,"\n"))
        individual1 <- name_combinations[icombination,1]
        individual2 <- name_combinations[icombination,2]
        if (site1==site2){
            # Indices of pair
            i1 = which(individual1 == names(data_set))  #index of ind1 on the data frame
            i2 = which(individual2 == names(data_set)) #index of ind2 on the data frame
        }else{
            i1 = which(individual1 == individual_names_A) #index of ind1 on the data frame
            i2 = which(individual2 == individual_names_B) + c_offset #index of ind2 on the data frame
        }
        # Extract data
        subdata <- cbind(data_set[,c("fs","dt")],data_set[,c(i1,i2)])
        names(subdata) <- c("fs","dt","Yi","Yj") # note fs not used
        krhat_hmm <- compute_rhat_hmm(frequencies, subdata$dt, cbind(subdata$Yi, subdata$Yj), epsilon)
        X[icombination-starty+1,1]=individual1
        X[icombination-starty+1,2]=individual2
        X[icombination-starty+1,3:4]=krhat_hmm

    }
    saveRDS(X, file = paste0(outputDir,'/',outFilePrefix,"_",starty,"_",endy,"_",gsub("\\.","p",sprintf("%.4f",f)),".RDS"))
    X
}

compute_rhat_hmm = function(frequencies, distances, Ys, epsilon){
    ll <- function(k, r) loglikelihood(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7)) #loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
    optimization <- optim(par = c(50, 0.5), fn = function(x) - ll(x[1], x[2]))
    rhat <- optimization$par
    return(rhat)
}

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
    Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
    return(Ys)
}


run_country=function(countryA,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    L=countryA  #readRDS(countryA)
    q=data.frame(L$X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(L$chroms[-j],matchy)))
    pos=L$pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}

run_2country=function(countryA,countryB,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    Z=countryA  #readRDS(countryA)
    M=countryB  #readRDS(countryB)
    X=cbind(Z$X,M$X)
    samps=c(Z$samps,M$samps)
    chroms=c(Z$chroms)
    pos=c(Z$pos)#
    q=data.frame(X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(chroms[-j],matchy)))
    pos=pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}

loglikelihood = function(k, r, Ys, f, gendist, epsilon, rho = 7.4 * 10^(-7)){
    loglikelihood_value = 0
    if (r < 0 | r > 1 | k < 0){
        return(log(0))
    }
    current_predictive = numeric(length = 2)
    current_predictive[1] = 1 - r
    current_predictive[2] = r
    current_filter = numeric(length = 2)
    ndata = nrow(Ys)
    maxnstates = ncol(f)
    for (idata in 1:ndata){
        nstates = 1
        # print(paste0("idata=",idata))
        while((nstates <= maxnstates) && (f[idata,nstates] > 1e-20)){
            nstates = nstates+1
        }
        lk0 = 0
        incr = 0
        # print(paste0("nstates=",nstates))
        # print(Ys[idata,])
        if(nstates>maxnstates) nstates=maxnstates
        for (g in 1:nstates){
            for (gprime in 1:nstates){
                if(gprime<=nstates){
                    incr = f[idata, g] * f[idata, gprime]
                    if (Ys[idata,1] == g){
                        incr = incr * (1 - (nstates - 1) * epsilon)
                    } else {
                        incr = incr*epsilon
                    }
                    if (Ys[idata,2] == gprime){
                        incr = incr * (1 - (nstates - 1) * epsilon)
                    } else {
                        incr = incr*epsilon
                    }
                    lk0 = lk0 + incr
                }
            }
        }
        lk1 = 0
        incr = 0

        for (g in 1:nstates){
          incr = f[idata, g]
          if (Ys[idata,1] == g){
            incr = incr * (1 - (nstates - 1) * epsilon)
          } else {
            incr = incr*epsilon
          }
          if (Ys[idata,2] == g){
            incr = incr * (1 - (nstates - 1) * epsilon)
          } else {
            incr = incr*epsilon
          }
          lk1 = lk1+incr
        }
        current_filter[1] = current_predictive[1] * lk0
        current_filter[2] = current_predictive[2] * lk1
        l_idata = current_filter[1] + current_filter[2]
        loglikelihood_value = loglikelihood_value + log(l_idata)
        if (idata < ndata-1){
          current_filter = current_filter / l_idata
          exp_ = exp(- k * rho * gendist[idata])
          a01 = r * (1 - exp_)
          a11 = r + (1-r) * exp_
          current_predictive[2] = current_filter[1] * a01 + current_filter[2] * a11
          current_predictive[1] = 1 - current_predictive[2]
        }
    }
    return(loglikelihood_value)
}





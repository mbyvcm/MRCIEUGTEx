
#' generate set of "query" snps using MRBase
#'
#' @param outcomes MRBase outcome index.
#' @param p Pvalue threshold on which to select SNPs.
#' @return Data.frame of SNPs, betas and alleles.
#' @export
get_query_snps_mrbase <- function(outcomes,p) {
  return(TwoSampleMR::extract_instruments(outcomes = outcomes, p1 = p)[,c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure')])
}

# given a set of dbSNP rsids, get genome coordinates (hg19) and extract these ranges from GTEx VCF
#' harmonise GTEx and query data and calculate polygenic risk score
#'
#' @param rsids Vector of dbSNP IDs.
#' @param gtex_vcf_dif Path to GTEx VCF file.
#' @return Data.frame: polygenic risk score calculated in GTEX samples.
#' @export
extract_query_snps_gtex <- function(rsids, gtex_vcf_dir) {

  message(paste0("reading ",length(rsids)," rsids..."))

  query_snps_granges <- update_query_snp_genome_coordinats(rsids)
  genome(query_snps_granges) <- "hg19"

  tabix_file <- TabixFile(gtex_vcf_dir)
  vcf <- readVcf(tabix_file, "hg19", query_snps_granges)

  rsids <- mcols(query_snps_granges)$RefSNP_id[subjectHits(findOverlaps(rowRanges(vcf), query_snps_granges))]

  dose_alt_allele <- geno(vcf)$DS
  type_alt_allele <- as.character(CharacterList(alt(vcf)))
  type_ref_allele <- as.character(ref(vcf))

  df <- data.frame(rsids, type_ref_allele, type_alt_allele, dose_alt_allele, stringsAsFactors = F)
  message(paste0(length(df$rsids), " rsids available for PRS"))
  return(df)
}

# harmonise GTEx and query data and calculate polygenic risk score
update_query_snp_genome_coordinats <- function(rsids) {
  message("getting genome coordinates for supplied rsids...")
  snps_filter <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, as.character(rsids), ifnotfound = 'drop')
  snps_filter <- renameSeqlevels(snps_filter, mapSeqlevels(seqlevels(snps_filter),"NCBI"))
  return(snps_filter)
}


# harmonise GTEx and query data and calculate polygenic risk score
calculate_prs_gtex <- function(query, gtex) {

  message('calculating polygenic risk score...')

  ## reformat so all betas are positive
  # positive beta, keep column order
  #ea <- query[query$beta.exposure > 0,]
  # negative beta, reverse alleles and invert beta
  #oa <- query[query$beta.exposure < 0,c('SNP','other_allele.exposure','effect_allele.exposure','beta.exposure')]
  #oa$beta.exposure <- abs(oa$beta.exposure)
  #names(oa) <- c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure')
  #query <- rbind(ea, oa)

  ## hamonise gtex data with query
  gtex_meta <- gtex[,1:3]
  mer <- merge(gtex_meta, query, all.x = T, by.x = 'rsids', by.y = 'SNP')

  mer <- strand_issues(mer)

  alt <- mer[mer$effect_allele.exposure == mer$type_alt_allele, c('rsids','beta.exposure')]
  ref <- mer[mer$effect_allele.exposure == mer$type_ref_allele, c('rsids','beta.exposure')]
  alt$allele <- "alt"
  ref$allele <- "ref"
  df <- rbind(alt,ref)
  df <- merge(x = df, y = gtex, by.x = "rsids")

  # reformat sampleid
  is_sample <- grep(names(df), pattern = '^GTEX')
  names(df)[is_sample] <- sapply(strsplit(names(df)[is_sample], split = '.', fixed = T), function(x) paste(x[1],x[2], sep = "_"))

  # get prs
  prs <- sapply(is_sample, calculatePRS, df = df)
  u <- mean(prs)
  sd <- sd(prs)
  prs <- (prs - u)/sd

  return(data.frame("IID" = names(df)[is_sample], "PRS" = prs, stringsAsFactors = F))
}


# identify SNPs where the alleles in the query are not both found in GTEx. The order (ref, alt) does not matter.
strand_issues <- function(mer) {

  exp <- mer$effect_allele.exposure == mer$type_alt_allele | mer$effect_allele.exposure == mer$type_ref_allele
  other <- mer$other_allele.exposure == mer$type_alt_allele | mer$other_allele.exposure == mer$type_ref_allele
  return(mer[exp&other,])
}


# Calculate PRS
calculatePRS <- function(x, df) {
  geno   <- df[,x]
  beta   <- df[,'beta.exposure']
  allele <- df[,'allele']
  prs <- 1/length(geno) * (
    sum(sum(geno[allele == "alt"] * beta[allele == "alt"]),
        sum((2-geno[allele == "ref"]) * beta[allele == "ref"]))
    )
  return(prs)
}

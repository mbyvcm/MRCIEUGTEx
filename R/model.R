
#' Run linear models
#'
#' @param geno Data.frame (or list of data.frames) of GTEx polygenic risk scores.
#' @param tissue Vector of tissues to be tested.
#' @param tx Vector of transcripts to be tested.
#' @param restrict.coding if TRUE (default) ignores non-coding transcripts.
#' @return List of data.frames by tissue; each dataframe contains summary statistics
#'   (p,se,b) per transcript.
#' @export
run_eqtl <- function(geno, tissue, tx = NULL, restrict_coding = T) {

  if (is.null(tissue)) {error("No tissue argument was included!"  )}
  if (is.null(geno))   {error("No genotype argument was included!")}

  res <- lapply(tissue, run_eqtl2, geno = geno, tx = tx, restrict_coding = restrict_coding)
  names(res) <- tissue
  return(res)
}


run_eqtl2 <- function(x, expression, geno, tx, restrict_coding = T) {

  message('')
  message(paste0(x))

  exp_name <- paste0('../../resources/',x,'_Analysis.expr.txt')
  cov_name <- paste0('../../resources/',x,'_Analysis.covariates.txt')

  # unzip expression & covar files from tar fot this tissue if not already
  if(!(file.exists(exp_name))) {
    message("  untaring..")
    untar(tarfile = expression_matrix_tar, files = paste0(x,"_Analysis.expr.txt"), exdir = '../../resources/')
  }

  if(!(file.exists(cov_name))) {
    message("  untaring..")
    untar(tarfile = covariate_matrix_tar, files = paste0(x,"_Analysis.covariates.txt"), exdir = '../../resources/')
  }

  # read files
  exp <- as.data.frame(read.table(exp_name, header = T, row.names = 1))
  cov <- as.data.frame(read.table(cov_name, header = T, row.names = 1))

  # if tx argument included, filter transcripts
  #if (!(is.null(tx))) {exp <- exp[,tx]}

  # if restricted argument TRUE, logical vector of which transcripts are
  # known protein coding
  #if (restrict_coding == T) {
  #  message("  Restricting analysis to known, protein coding genes")
  #  exp <- exp[,restrict_protein_coding(gtf_path = gtf_path, tx = names(exp))]
  #}

  # report basic stats
  message(paste0("  Number of transcripts: ",dim(exp)[1]))
  Sys.sleep(1)
  message(paste0("  Expression matrix samples: ",dim(exp)[2]))
  Sys.sleep(1)
  message(paste0("  Covariate matrix samples: ",dim(cov)[2]))

  # ensure expression and covariate samples are in the same order
  if (length(names(exp)) != length(names(cov))) { stop("covariate and expression  matricies have different number of samples!")}
  if (all(names(exp) != names(cov))) { stop("covariate and expression matricies have different sample order!")}

  # Load covariate data matrixEQTL
  covariates_file_name = cov_name
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = "\t"
  cvrt$fileOmitCharacters = "NA"
  cvrt$fileSkipRows = 1
  cvrt$fileSkipColumns = 1
  cvrt$fileSliceSize = 2000
  cvrt$LoadFile(covariates_file_name)

  # Load expression data matrixEQTL
  expression_file_name = exp_name
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000
  gene$LoadFile(expression_file_name)

  ret <- lapply(names(geno), function(prs) {

    geno <- geno[[prs]]

    # geno samples
    genoIID <- gsub(geno$IID, pattern = '_', replacement = '.')

    # write geno for samples in cov/exp
    g <- geno[match(names(cov), genoIID),]
    write.table(t(g),'./snps.txt', quote = F, row.names = F, col.names = F, sep = "\t")

    df <- run_matrixEQTL(exp_name, cov_name)
    return(df)
    })

  names(ret) <- names(geno)
  return(ret)
}


run_matrixEQTL <- function(gene,cvrt) {

  # SNP data
  snps_file_name = "./snps.txt"
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 0
  snps$fileSliceSize = 2000
  snps$LoadFile( snps_file_name )

  output_file_name = tempfile()
  useModel = modelLINEAR
  pvOutputThreshold = 1
  errorCovariance = numeric()

  me <- Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  tx  <- me$all$eqtls$gene
  p   <- me$all$eqtls$pvalue
  b   <- me$all$eqtls$beta
  fdr <- me$all$eqtls$FDR

  return(data.frame(tx, b, p, fdr, stringsAsFactors = F, row.names = 1))
}

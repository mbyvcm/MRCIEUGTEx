
#' Run linear models
#'
#' @param geno Data.frame of GTEx polygenic risk scores.
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
  exp <- as.data.frame(t(read.table(exp_name, header = T, row.names = 1)))
  cov <- as.data.frame(t(read.table(cov_name, header = T, row.names = 1)))

  # if tx argument included, filter transcripts
  if (!(is.null(tx))) {exp <- exp[,tx]}

  # if restricted argument TRUE, logical vector of which transcripts are
  # known protein coding
  if (restrict_coding == T) {
    message("  Restricting analysis to known, protein coding genes")
    exp <- exp[,restrict_protein_coding(gtf_path = gtf_path, tx = names(exp))]
  }

  # report basic stats
  message(paste0("  Number of transcripts: ",dim(exp)[2]))
  Sys.sleep(1)
  message(paste0("  Expression matrix samples: ",dim(exp)[1]))
  Sys.sleep(1)
  message(paste0("  Covariate matrix samples: ",dim(cov)[1]))

  # exp and cov sampling should be the same
  e <- rownames(exp)
  c <- rownames(cov)

  # merge PRS to covariates
  cov$IID <- c
  geno$IID <- gsub(geno$IID, pattern = '_', replacement = '.')
  geno_cov <- merge(x = geno, y = cov, by = "IID")

  # ensure expression and covariate samples are in the same order
  exp <- exp[e %in% geno_cov$IID,]
  geno_cov <- geno_cov[geno_cov$IID %in% e,]
  exp <- exp[order(rownames(exp)),]
  geno_cov <- geno_cov[order(geno_cov$IID),]

  # error out if the order and/or lengty of samples in covariate and expression data is different
  if (sum(rownames(exp) == geno_cov$IID) != length(rownames(exp))) {warning("expression and geno matrix mismatch!")}

  message(paste0("  Genotypes matrix samples: ",dim(geno_cov)[1]))

  # run regression
  message("")
  message("  Running model")

  res <- pblapply(exp, runlm, geno_cov)
  res <- do.call(rbind.data.frame, res)
  names(res) <- c('b','se','p')
  return(res)
}


runlm <- function(x, g) {

  model <- fastLm(as.matrix(g[,-1]), x)
  sum   <- summary(model)

  pvalue <- coefficients(sum)["PRS","p.value"]
  beta   <- coefficients(sum)["PRS","Estimate"]
  se     <- coefficients(sum)["PRS","StdErr"]

  return(c(beta,se,pvalue))
}

#' Returns list of tissues in GTEX available for analysis.
#'
#' @param expression_matrix_tar Path to tar file containing expression matricies by tissue.
#' @return Data.frame: polygenic risk score calculated in GTEX samples.
#' @export
available_tissues <- function(expression_matrix_tar = expression_matrix_tar) {

  file <- "../../resources/tissue_index.txt"

  if(file.exists(file)) {
    tissues <- scan(file = file, what = 'character')
    names(tissues) <- seq_len(length(tissues))

  } else {
    tissues <- gsub(
      x = untar(expression_matrix_tar, list = T),
      pattern = '_Analysis.expr.txt',
      replacement = '')
    names(tissues) <- seq_len(length(tissues))
    write(tissues, file)
  }
  return(tissues)
}


restrict_protein_coding <- function(gtf_path = gtf_path, tx) {

  gtf <- read.table(gtf_path, skip = 5, sep = '\t', stringsAsFactors = F)
  geneid <- gsub(stringr::str_extract(string = gtf[,9], pattern = "gene_id ENSG\\d+.\\d+"), pattern = "gene_id ", replacement = "")
  genetype <- gsub(stringr::str_extract(string = gtf[,9], pattern = "gene_type \\w+"), pattern = "gene_type ", replacement = "")
  genestatus <- gsub(stringr::str_extract(string = gtf[,9], pattern = "gene_status \\w+"), pattern = "gene_status ", replacement = "")

  return(tx %in% geneid[genetype == 'protein_coding' & genestatus == 'KNOWN'])
}

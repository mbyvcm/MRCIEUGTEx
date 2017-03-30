
#' Generate annotated list of top hits
#'
#' @param x List object returned by "run_eQTL".
#' @param fdr Should FDR be used.
#' @param p.thresh Threshold for annotation gene names.
#' @return Data.frame of genes passing threshold, annotated with gene names.
#' @export
extract_top_hits <- function(x, fdr = T, p.thresh = 0.05) {

  if (fdr != 1) {
    th <- lapply(x, function(x) {
      x[x[,'p'] <= p.thresh,]
    })
  } else {
    th <- lapply(x, function(x) {
      x[x[,'fdr'] <= p.thresh,]
    })
  }

  th <- th[unlist(lapply(th, function(x) dim(x)[1] > 0))]

  output <- lapply(names(th), function(x) {

    tissue <- x
    df     <- th[[tissue]]
    tx     <- rownames(df)
    txrfm  <- gsub(tx, pattern = '\\.\\d+', replacement = '')
    symbol <- try(add_gene_names(ids = txrfm), silent = T)

    # reorder symbol to conform to input
    if (class(symbol) != "try-error") {
      symbol <- symbol[match(txrfm, names(symbol))]
    } else {symbol <- rep("NA", length(txrfm)); names(symbol) <- txrfm}

    p      <- as.numeric(df[,'p'])
    fdr    <- as.numeric(df[,'fdr'])
    b      <- as.numeric(df[,'b'])

    data.frame(tissue, tx, symbol, b, fdr, p, row.names = NULL)
  })

  return(do.call(rbind, output))
}


#' Render volcanoe plot for each element in output list.
#'
#' @param x List object returned by "run_eQTL"
#' @param fdr Should FDR be used.
#' @param p.thresh Threshold for annotation gene names.
#' @param outdir Output directory for plots.
#' @return volcanoe plots
#' @export
volcanoeplot <- function(x, fdr = T, p.thresh = 0.05, outdir = '../../output/') {

  tissues <- names(x)

  lapply(tissues, function(tissue) {

    message(paste0("volcanoe plot ", tissue))

    df  <- x[[tissue]]
    row.names(df) <- gsub(row.names(df), pattern = '\\.\\d+', replacement = '')

    if (fdr != 1) {
      sig_ids   <- row.names(df[df$p <= p.thresh,]); pass <- df$p <= p.thresh
    } else {sig_ids <- row.names(df[df$fdr <= p.thresh,]); pass <- df$fdr <= p.thresh}


    if (length(sig_ids) > 0) {
      try(symbols   <- add_gene_names(sig_ids), silent = T)
      if (class(symbols) == "try-error") {
        df$symbol <- NA
      } else {df$symbol <- symbols[match(x = rownames(df), table = names(symbols))]}
    } else {df$symbol <- NA }

    plot <- ggplot(df, aes(x= b, y= -log10(p), colour = pass, label=as.character(symbol))) +
      geom_point(alpha=0.4, size=1.75) +
      geom_text(aes(label=ifelse(pass,as.character(symbol),'')),hjust=0,vjust=0) +
      theme(legend.position="none") +
      ggtitle(paste0(tissue,":",p.thresh))

    ggsave(filename = paste0(outdir,"_",tissue,"_VOL.pdf"), plot = plot, device = 'pdf')
    return(plot)
  })
}

add_gene_names <- function(ids) {
  # one-to-many output
  output <- select(x = org.Hs.eg.db, columns = 'SYMBOL', keys = as.character(ids), keytype = 'ENSEMBL')
  # make one-to-one
  return(unlist(lapply(split(output, output$ENSEMBL), function(x) paste(x[,'SYMBOL'], collapse = ','))))
}

#' Returns transcripts surpassing given FDR threshold in each tissue
#'
#' @param x List object returned by "run_eQTL"
#' @param fdr FDR threshold
#' @return Data.frame of top hits
#' @export
extract_top_hits <- function(x, fdr) {

  th <- lapply(x, function(x) {
    x[p.adjust(x[,'p'], method = 'BH') <= fdr,]
  })

  th <- th[unlist(lapply(th, function(x) dim(x)[1] > 0))]

  output <- lapply(names(th), function(x) {

    tissue <- x
    df     <- th[[tissue]]
    tx     <- rownames(df)
    txrfm  <- gsub(tx, pattern = '\\.\\d+', replacement = '')
    symbol <- try(add_gene_names(ids = txrfm), silent = T)

    # reorder symbol to conform to input
    if (class(symbol) != "try-error") {
      symbol <- symbol[match(txrfm, names(symbol))]
    } else {symbol <- rep("NA", length(txrfm)); names(symbol) <- txrfm}

    p      <- as.numeric(df[,'p'])
    se     <- as.numeric(df[,'se'])
    b      <- as.numeric(df[,'b'])

    data.frame(tissue, tx, symbol, b, se, p, row.names = NULL)
  })

  return(do.call(rbind, output))

}

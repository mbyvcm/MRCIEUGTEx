qc <- function(x, outdir = './') {
  names <- names(x)
  lapply(names, plotqq, list = x, outdir = outdir)
}


plotqq <- function(name,list, outdir) {

  df <- list[[name]]

  pdf(file = paste0(outdir,name,"_QQ.pdf"), paper = 'a4')
  qqplot(pnorm(rnorm(length(df$p), mean = 0, sd = 1)),df$p, xlab = "expected", ylab = "observed")
  abline(coef= c(0,1), col = 'red')
  dev.off()
}

#' Render volcanoe plot for each element in output list.
#'
#' @param x List object returned by "run_eQTL"
#' @return volcanoe plots
#' @export
volcanoeplot <- function(x, fdr, outdir = '../../output/') {

  tissues <- names(x)

  lapply(tissues, function(tissue) {

    message(paste0("volcanoe plot ", tissue))

    df  <- x[[tissue]]
    row.names(df) <- gsub(row.names(df), pattern = '\\.\\d+', replacement = '')
    sig_ids   <- row.names(df[p.adjust(df$p, method = "BH") < fdr,])

    if (length(sig_ids) > 0) {
      try(symbols   <- add_gene_names(sig_ids), silent = T)
      if (class(symbols) == "try-error") {
        df$symbol <- NA
      } else {df$symbol <- symbols[match(x = rownames(df), table = names(symbols))]}
    } else {df$symbol <- NA }

    plot <- ggplot(df, aes(x= b, y= -log10(p), colour = p.adjust(p, method = "BH") <= fdr, label=as.character(symbol))) +
      geom_point(alpha=0.4, size=1.75) +
      geom_text(aes(label=ifelse(p.adjust(p, method = "BH")<= fdr,as.character(symbol),'')),hjust=0,vjust=0) +
      theme(legend.position="none") +
      ggtitle(paste0(tissue,"_FDR:",fdr))

    ggsave(filename = paste0(outdir,"_",tissue,"_VOL.pdf"), plot = plot, device = 'pdf')
    return(plot)
  })
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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

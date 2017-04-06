#' Plot distances between transcript profiles
#'
#' @param x List object returned by "run_eQTL".
#' @param collapse If true all tissues x traits are unlisted and plotted together
#' @return List of plots
#' @export
distanceByPRS <- function(all, collapse = F) {

  if (collapse == 1) {
    l <- unlist(all[['models']], recursive = F)
    return(distanceByPRS2(l, collapse = T))
  } else {

    tissues <- names(all[['models']])
    out <- lapply(tissues, function(tissue) {
      message(tissue)
      return(distanceByPRS2(all[['models']][[tissue]], collapse = F))
    })
    names(out) <- tissues
    return(out)
  }
}

distanceByPRS2 <- function(l, collapse) {

  if (length(l) < 3) {stop('n must be > 2 to calculate eigen')}

  df <- Reduce(
    function(dtf1, dtf2) merge(dtf1, dtf2, by = 'iid', all = TRUE),
    lapply(l, function(trait) {data.frame('iid' = rownames(trait),'b' = trait$b)}
    ))

  rownames(df) <- df$iid
  mat <- t(as.matrix(df[,-1]))
  rownames(mat) <- names(l)

  dist <- dist(mat, diag = T, method = 'e')
  mds <- cmdscale(dist, eig = TRUE, k = 2)
  df <- data.frame("C1" = mds$points[,1], "C2" = mds$points[,2])

  if (collapse == 1) {

    Tissue <-     gsub(rownames(df), pattern = '^(\\S+)\\.(\\S+)', replacement = '\\1')
    PRS <-  gsub(rownames(df), pattern = '^(\\S+)\\.(\\S+)', replacement = '\\2')

    p <- ggplot(df, aes(x = C1, y = C2, color = PRS, label = Tissue)) +
      geom_text(aes(label=Tissue))

    return(p)
  } else {

    PRS <- rownames(df)

    p <- ggplot(df, aes(x = C1, y = C2, color = PRS, label = PRS)) +
      geom_text(aes(label=PRS))

    return(p)

  }
}
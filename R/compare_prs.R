#' Plot distances between transcript profiles
#'
#' @param x List object returned by "run_eQTL".
#' @param collapse If true all tissues x traits are unlisted and plotted together
#' @return List of plots
#' @export
distanceByPRS <- function(all, collapse = F, plot.eig = F) {

  if (collapse == 1) {
    l <- unlist(all[['models']], recursive = F)
    return(distanceByPRS2(l, collapse = T, plot.eig = plot.eig))
  } else {

    tissues <- names(all[['models']])
    out <- lapply(tissues, function(tissue) {
      message(tissue)
      return(distanceByPRS2(all[['models']][[tissue]], collapse = F, plot.eig = plot.eig, tissue = tissue))
    })
    names(out) <- tissues
    return(out)
  }
}

#' Plot distances between transcript profiles from unlisted
#' 
#' @param l
#' @param collapse
#' @param plot.eig
#' @param tissue
#' @return List of plots
#' @export
distanceByPRS2 <- function(l, collapse, plot.eig, tissue = NULL) {

  if (length(l) < 3) {stop('n must be > 2 to calculate eigen')}

  df <- Reduce(
    function(dtf1, dtf2) merge(dtf1, dtf2, by = 'iid', all = TRUE),
    lapply(l, function(trait) {data.frame('iid' = rownames(trait),'b' = trait[,'b'])}
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

    isBrain <- as.factor(as.numeric(grepl(Tissue, pattern = 'Brain')))
    
    p <- ggplot(df, aes(x = C1, y = C2, color = PRS, label = Tissue)) +
      geom_point(aes(shape = isBrain)) +
      geom_text(aes(label=Tissue), size = 2, colour = 'black', nudge_y = -0.2)
    

    if (plot.eig == 1) {
      p2 <- autoplot(prcomp(x = mat), loadings = T, loadings.label = T, label = T)
      return(list(p,p2))
    } else {
      return(p)
    }



  } else {

    PRS <- rownames(df)

    p <- ggplot(df, aes(x = C1, y = C2, color = PRS, label = PRS)) +
      geom_point() +
      geom_text(aes(label=PRS))

    if (plot.eig == 1) {
      p2 <- autoplot(prcomp(x = mat), loadings = T, loadings.label = T, label = T, main = tissue)
      return(list(p,p2))
    } else {
      return(p)
    }

  }
}


#' Generate GO object
#'
#' @param GOTERMS Vector of GO Terms (GO:#######).
#' @return List Ensembl gene IDS by GO Term
#' @export
generateSetListGO <- function(GOTERMS) {

  x <- select(org.Hs.eg.db, columns = "ENSEMBL", keys = GOTERMS, keytype = 'GO')
  x <- x[x$GO %in% GOTERMS,c('GO','ENSEMBL')]
  x <- lapply(split(x,x$GO), function(x) x[,'ENSEMBL'])

  if (all(lapply(x, length) > 10) != 1) {warning('Some sets contain fewer than ten genes!')}

  return(x)
}


#' Run Pathways Analysis
#'
#' @param x List object returned by "run_eQTL".
#' @param setList Generated using included function
#' @return List pathways scores
#' @export
runPathwayAnalysis <- function(x, setList) {

  out <- lapply(output[['models']], function(tissue) {
    out <- lapply(tissue, function(trait) {

      gene <- gsub(rownames(trait), pattern = '\\.\\d+$', replacement = '')
      z    <- qnorm(trait[,'p'])
      return(runPathwayAnalysis2(gene,z, setList))
    })
    return(out)
  })
  return(list('models' = out))
}


runPathwayAnalysis2 <- function(gene, z, setList) {

  out <- lapply(setList, function(term) {

    # 1 if gene is in catagory, 0 otherwise
    catVar <- as.numeric(gene %in% term)
    sum    <- summary(lm(z ~ catVar))
    beta   <- coefficients(sum)['catVar','Estimate']
    p      <- coefficients(sum)['catVar','Pr(>|t|)']
    return(c("n_cat" = dim(term)[1], "n_avail" = sum(catVar), 'b' = beta, 'p' = p))

  })
  return(t(as.data.frame(out)))
}

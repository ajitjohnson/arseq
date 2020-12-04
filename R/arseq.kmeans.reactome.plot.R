#' @title Reactome Pathway Enrichment of Kmeans clustering output
#' @description Reactome Pathway Enrichment of Kmeans clustering output
#' @param clusters List of gene clusters returned by Kmeans clustering algorithm
#' @return DotPlot of Reactome Enrichment
#' @import ReactomePA
#' @import clusterProfiler
#' @examples
#' arseq.kmeans.reactome.plot (clusters)
#' @export

arseq.kmeans.reactome.plot <- function(clusters, save.plot=FALSE, save.dir=getwd()){

  # clusters from Kmeans
  clusters <- data.frame(clusters)
  gene2entrez <- function(cluster){
    require(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    genes = row.names(clusters[clusters[,1] == cluster,,drop=F])
    ENT <- AnnotationDbi::select (hs, keys = genes,columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    ENT <- ENT[complete.cases(ENT), ]
    entid <- ENT$ENTREZID
    return(entid)
  }
  # Covert all gene symbols to entrezid
  gn <- lapply(unique(clusters[,1]), gene2entrez)
  names(gn) <- unique(clusters[,1])
  # Perform reactome enrichment
  require(ReactomePA)
  res <- compareCluster(gn, fun="enrichPathway")

  # Save or plot the DotPlot
  print("Plotting Reactome Pathway Enrichment")
  if (isTRUE(save.plot)){
    pdf(paste(save.dir,"Reactome Enrichment Plot.pdf",sep = ""),width=16,height=8,paper='special')
    plot(dotplot(res))
    dev.off()
  } else {
    dotplot(res)
  }
}



#' @title Reactome enrichment Analysis
#' @description Performing reactome enrichment analysis on the significantly differentially expressed genes.
#' @param deg Differentially expressed genes dataframe returned by DESeq2 analysis
#' @param Padj Numeric. Adjusted P Value cut off to define the differentially expressed genes. Default = 0.05
#' @return Enrichment dataframe
#' @import ReactomePA
#' @examples
#' \dontrun{
#' reactome.enrich <- arseq.go.enrich (example_deg)
#' }
#' @export

arseq.reactome.enrich <- function(deg,Padj=0.05,plot=FALSE,save.plot=FALSE,save.dir=getwd()){
  print("Performing Reactome enrichment analysis between the constrast groups")

  ranked_list <- data.frame(deg[which(deg$padj <= Padj), ])
  r_list <- ranked_list$log2FoldChange
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  ranked.list <- r_list[is.finite(r_list)]

  # Convert gene names to entrez ID
  require(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  ENT <- select (hs, keys = names(ranked.list), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  ENT <- ENT[complete.cases(ENT), ]
  entid <- ENT$ENTREZID

  # Rename genenames to entrez id
  ranked.list.entrez <- ranked.list
  names(ranked.list.entrez) <- entid

  # Perfoem enrichment
  arseq.reactome.enrich <- enrichPathway(gene=names(ranked.list.entrez), pvalueCutoff=1, readable=T)
  arseq.reactome.enrich.df <- as.data.frame(arseq.reactome.enrich)
  arseq.reactome.enrich.df <- arseq.reactome.enrich.df[,-1]

  if (isTRUE(plot)){
    dotplot(arseq.reactome.enrich, showCategory=20)
    if (isTRUE(save.plot)){
      # Dot Plot
      pdf(paste(save.dir,"Reactome Enrichment Dot Plot.pdf",sep = ""),width=16,height=8,paper='special')
      plot(dotplot(arseq.reactome.enrich, showCategory=20))
      dev.off()
      # emapplot
      pdf(paste(save.dir,"Reactome Enrichment EMAP Plot.pdf",sep = ""),width=16,height=8,paper='special')
      plot(emapplot(arseq.reactome.enrich))
      dev.off()
      # cnetplot
      pdf(paste(save.dir,"Reactome Enrichment CNET Plot.pdf",sep = ""),width=16,height=8,paper='special')
      plot(cnetplot(arseq.reactome.enrich, categorySize="pvalue", foldChange=ranked.list.entrez))
      dev.off()

    }
  }

  # Return the results
  return(arseq.reactome.enrich.df)
}

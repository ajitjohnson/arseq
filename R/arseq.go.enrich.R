#' @title Gene Ontology (GO) term enrichment
#' @description Performing GO term enrichment analysis on the significantly differentially expressed genes.
#' @param deg Differentially expressed genes dataframe returned by DESeq2 analysis
#' @param Padj Numeric. Adjusted P Value cut off to define the differentially expressed genes. Default = 0.05
#' @return Enrichment dataframe
#' @import goseq
#' @examples
#' \dontrun{
#' go.enrich <- arseq.go.enrich (example_deg)
#' }
#' @export

arseq.go.enrich <- function(deg,Padj=0.05){
  print("Performing GO enrichment analysis between the constrast groups")
  diff_genes <- data.frame(deg[which(deg$padj <= Padj), ])
  gene.vector <- as.integer(row.names(data.frame(deg))%in%row.names(diff_genes))
  names(gene.vector) <- row.names(data.frame(deg))
  # Prepare for goseq (gene length data for the DEG's)
  pwf <- suppressWarnings(nullp(gene.vector,"hg19","geneSymbol"))
  # Using the Wallenius approximation
  go.enrich=goseq(pwf,"hg19","geneSymbol")
  return(go.enrich)
}

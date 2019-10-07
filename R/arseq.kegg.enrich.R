#' @title KEGG pathway enrichment
#' @description Performing KEGG pathway enrichment analysis on the significantly differentially expressed genes.
#' @param deg Differentially expressed genes dataframe returned by DESeq2 analysis
#' @param kegg.compare Character, which comparison scheme to be used: 'paired', 'unpaired', '1ongroup', 'as.group'. 'paired' is the default, ref and samp are of equal length and one-on-one paired by the original experimental design; 'as.group', group-on-group comparison between ref and samp; 'unpaired' (used to be '1on1'), one-on-one comparison between all possible ref and samp combinations, although the original experimental design may not be one-on-one paired; '1ongroup', comparison between one samp column at a time vs the average of all ref columns.
#' @return Enrichment dataframe
#' @importFrom AnnotationDbi mapIds
#' @import gage
#' @examples
#' \dontrun{
#' kegg.enrich <- arseq.kegg.enrich (example_deg)
#' }
#' @export

arseq.kegg.enrich <- function(deg, kegg.compare="as.group"){
  print("Performing KEGG pathway analysis between the constrast groups")
  # Find Entrez ID
  deg$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(deg),
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  # Intialize the gage function
  foldchanges <- deg$log2FoldChange
  names(foldchanges) <- deg$entrez
  # Get the results
  kegg.sets.hs <- kegg.sets.hs
  keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE, compare=kegg.compare)
  # List to store results
  pathways <- list()
  pathways$up_regulated_pathways <- data.frame(keggres[["greater"]])
  pathways$down_regulated_pathways <- data.frame(keggres[["less"]])
  pathways$foldchanges <- foldchanges
  return(pathways)
}

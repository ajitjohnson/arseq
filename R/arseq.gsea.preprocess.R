#' @title GSEA pre-processing step
#' @description Function takes in the differentially expressed genes from DESEq2 and produces a ranked list for GSEA.
#' @param deg Differentially expressed genes dataframe returned by DESeq2 analysis.
#' @return ranked list of genes.
#' @examples
#' ranked.list <- arseq.gsea.preprocess (deg=example_deg)
#' @export

arseq.gsea.preprocess <- function(deg){
  # Create a ranked list of genes from deg
  ranked_list <- data.frame(deg[!is.na(deg$padj),])
  ranked_list$score <- -log(ranked_list$padj)*sign(ranked_list$log2FoldChange)
  ranked_list <- ranked_list[order(-ranked_list$score),]
  r_list <- ranked_list$score
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  ranked.list = r_list[is.finite(r_list)]
  return(ranked.list)
}

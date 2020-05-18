#' @title Most Variable Genes
#' @description Identifying the most variable genes in the dataset
#' @param data DESeq2 object or normalized expression dataframe
#' @param variable.genes numeric: The number of most variable genes to be identified. By default, the program identifies the top 1000 most variable genes.
#' @param dds.object Logical parameter indicating if the data is a DESeq2 object. Default=TRUE
#' @return The most variable genes
#' @import DESeq2
#' @importFrom stats var
#' @examples
#' mvg <- arseq.mvg (example_dds)
#' @export

arseq.mvg <- function(data,variable.genes=variable.genes,dds.object=TRUE){
  # Normalize the data if the data is a DESeq2 object
  if (isTRUE(dds.object)){data <- log2(counts(data, normalized=TRUE)+1)}
  # MVG
  print("Identifying the most variable genes")
  data.var <- apply(data, 1, stats::var)
  arseq.mvg <- data[order(data.var, decreasing = TRUE)[1:variable.genes],]
  return(arseq.mvg)
}

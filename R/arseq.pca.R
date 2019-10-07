#' @title Principal Component Analysis
#' @description Plots any two specified principal components. By default it plots the first two principal components.
#' @param dds DESeq2 object
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @return Eigen vectors of the Principal components
#' @import DESeq2
#' @importFrom matrixStats rowVars
#' @importFrom stats prcomp
#' @examples
#' \dontrun{
#' pca.ev <- arseq.pca (example_dds)
#' }
#' @export

arseq.pca = function(dds, ntop=500){
  print("Performing a PCA analysis")
  # Normalize data
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  # calculate the variance for each gene
  rv <- rowVars(vsd@assays$data[[1]])
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(vsd@assays$data[[1]][select,])
  pca <- data.frame(pca$x)
  arseq.pca <- pca[order(-pca[,1]),]

  # return the dataframe
  return(arseq.pca)
}

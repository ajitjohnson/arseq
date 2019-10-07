#' @title Euclidean distance
#' @description Calculate the Euclidean distance between samples.
#' @param dds DESeq2 object
#' @return Euclidean distance matrix
#' @import DESeq2
#' @importFrom stats dist
#' @examples
#' euclid.dist <- arseq.euclid.dist (example_dds)
#' @export

# Sample to sample distance
arseq.euclid.dist <- function (dds){
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  print("Calculating the Euclidean distance")
  euclid.dist <- dist(t(vsd@assays$data[[1]]))
  return(euclid.dist)
}

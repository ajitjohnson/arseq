#' @title Poisson distance
#' @description Calculate the Poisson distance between samples.
#' @param dds DESeq2 object
#' @return Poisson distance matrix
#' @import DESeq2
#' @import PoiClaClu
#' @examples
#' poisd.dist <- arseq.poisd.dist (example_dds)
#' @export

arseq.poisd.dist <- function (dds){
  print("Calculating the Poisson distance")
  poisd <- PoissonDistance(t(counts(dds)))
  poisd.dist <- poisd$dd
  attr(poisd.dist, "Labels") <- colnames(dds)
  return(poisd.dist)
}

#' @title Euclidean distance heatmap
#' @description Plots a heatmap of the distance between samples using pheatmap.
#' @param euclid.dist Euclidean distance matrix.
#' @return Euclidean distance heatmap.
#' @import pheatmap
#' @import viridis
#' @examples
#' euclid.dist <- arseq.euclid.dist (example_dds)
#' euclid.plot <- arseq.euclid.dist.plot (euclid.dist)
#' @export

arseq.euclid.dist.plot <- function(euclid.dist){
  euclid.dist.plot <- pheatmap(as.matrix(euclid.dist),
           clustering_distance_rows = euclid.dist,
           clustering_distance_cols = euclid.dist,
           border_color = NA,
           color = inferno(50),
           main = "Euclidean distance between samples")
  return(euclid.dist.plot)
}


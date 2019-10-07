#' @title Poisson distance heatmap
#' @description Plots a heatmap of the Poisson distance between samples using pheatmap.
#' @param poisd.dist Poisson distance matrix.
#' @return Poisson distance heatmap.
#' @import pheatmap
#' @import viridis
#' @examples
#' poisd.dist <- arseq.poisd.dist (example_dds)
#' poisd.plot <- arseq.poisd.dist.plot (poisd.dist)
#' @export

arseq.poisd.dist.plot <- function(poisd.dist){
  poisd.dist.plot <- pheatmap(as.matrix(poisd.dist),
                               clustering_distance_rows = poisd.dist,
                               clustering_distance_cols = poisd.dist,
                               border_color = NA,
                               color = inferno(50),
                               main = "Poisson distance between samples")
  return(poisd.dist.plot)
}

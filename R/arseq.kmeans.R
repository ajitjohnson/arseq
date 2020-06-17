#' @title K Means clustering
#' @description K means clustering to identify gene modules
#' @param data Normalized expression dataframe
#' @param kmeans integer: The number of clusters to be returned
#' @return Kmeans clusters
#' @importFrom stats kmeans
#' @examples
#' clusters <- arseq.kmeans (data, kmeans=10)
#' @export

arseq.kmeans <- function(data, kmeans=kmeans){
  print("Clustering the genes using the Kmeans clustering algorithm")
  clusters <- kmeans(data, centers=kmeans)
  arseq.kmeans <- clusters$cluster
  return(arseq.kmeans)
}

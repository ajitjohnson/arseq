#' @title K Means clustering Heatmap
#' @description Heatmap of the Most Variable genes divided into clusters
#' @param data Normalized expression dataframe
#' @param clusters List of gene clusters returned by Kmeans clustering algorithm
#' @import ComplexHeatmap
#' @import circlize
#' @import RColorBrewer
#' @import ReactomePA
#' @return Heatmap
#' @examples
#' arseq.kmeans.plot (data,clusters=kmeansclusters,metadata= example_meta,intgroup="treatment")
#' @export

arseq.kmeans.plot <- function(data, clusters, metadata,
                              intgroup="arseq.group", save.plot=FALSE, save.dir=getwd()){

  # Scale the data
  scaled_data <- t(scale(t(data)))
  clusters <- data.frame(clusters)

  # Set the colors for heatmap
  h1_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#F8F6B8","#F5A15B","#B11E4B"))

  # Add Row annotation (clusters)
  row_annotation <- clusters
  colnames(row_annotation) <- c('cluster')

  # Add column annotation (Sample group)
  col_ann <- data.frame(as.character(metadata[,intgroup]))
  colnames(col_ann) <- c('Groups')
  # Create color pallete
  colourCount <- length(unique(col_ann[,1]))
  getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
  col_color <- getPalette(colourCount)
  #col_color <- brewer.pal(length(unique(col_ann[,1])),"Dark2") #BrBG
  names(col_color) <- unique(col_ann[,1])
  col_color <- list(Groups = col_color)
  col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))

  # Plot the heatmap
  h1 <- Heatmap(as.matrix(scaled_data), col = h1_col, cluster_columns = F, cluster_rows = T, show_row_names = FALSE,
                name='z score', top_annotation=col_Ann,
                row_split = factor(row_annotation$cluster, levels = unique(row_annotation$cluster)))

  # Save or plot the heatmap
  if (isTRUE(save.plot)){
    pdf(paste(save.dir,"Heatmap.pdf",sep = ""),width=6,height=6,paper='special')
    draw(h1)
    dev.off()
  } else {
    draw(h1)
  }
}

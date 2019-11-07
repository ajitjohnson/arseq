#' @title Multidimensional scaling (MDS) Analysis
#' @description Plots the Multidimensional scaling (MDS) factors.
#' @param dds DESeq2 object
#' @param intgroup interesting groups: a character vector of names in colData(x) to use for grouping.
#' @return ggplot of the two selected principal components
#' @import DESeq2
#' @import ggrepel
#' @import ggplot2
#' @importFrom stats dist cmdscale
#' @importFrom dplyr tbl_df "%>%"
#' @examples
#' \dontrun{
#' mds.plot <- arseq.mds.plot (example_dds, intgroup="treatment")
#' }
#' @export

arseq.mds.plot = function(dds,intgroup="arseq.group"){
  # normalize data
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  # Calculate MDS
  print("Performing a multidimensional scaling (MDS) analysis")
  mds <- as.data.frame(vsd@colData)  %>% cbind(cmdscale(as.matrix(dist(t(vsd@assays$data[[1]])))))
  # Plot MDS
  ggplot(mds, aes_string(x = mds$`1`, y = mds$`2`, color = intgroup)) +
    geom_point(size = 3) + geom_text_repel(aes(label = rownames(mds)),size = 3)+
    theme_classic()+ coord_fixed()+
    xlab ("MDS_1")+ ylab ("MDS_2") + ggtitle("Multidimensional scaling (MDS) Plot")+
    theme(plot.title = element_text(hjust = 0.5))
}

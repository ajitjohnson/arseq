#' @title Most Variable Genes
#' @description Identifying the most variable genes in the dataset
#' @param mvg Expression matrix of the most variable genes identified using the function 'arseq.mvg'
#' @param metadata A CSV file with information regarding the samples.
#' @param intgroup interesting groups: a character vector of names in colData(x) to use for grouping.
#' @param save.plot Logical. Argument to state if the plot needs to be saved in disk. Default: FALSE
#' @param save.dir Location to save the plot when 'save.plot=TRUE'. Default: Working Directory.
#' @return Heatmap of the most variable genes
#' @import RColorBrewer
#' @import grDevices
#' @import gplots
#' @importFrom stats hclust
#' @examples
#' mvg <- arseq.mvg (example_dds)
#' mvg.plot <- arseq.mvg.plot (mvg[1:100,],metadata= example_meta,intgroup="treatment")
#' @export

arseq.mvg.plot <- function(mvg,metadata,intgroup="arseq.group", save.plot=FALSE, save.dir=getwd()){
  # Heatmap of the most variable genes
  my_group = factor(metadata[,intgroup])
  my_col = brewer.pal(length(levels(factor(metadata[,intgroup]))), "Set1")[my_group]
  coul = colorRampPalette(brewer.pal(9,"OrRd"))(10)
  # heatmap
  if (isTRUE(save.plot)){
    pdf(paste(save.dir,"Heatmap of the most variable genes.pdf",sep = ""),width=6,height=6,paper='special')
    heatmap.2(as.matrix(mvg), Colv = NA, scale="row",
              col =coul, hclustfun = hclust,trace="none",
              dendrogram= "row",ColSideColors=my_col,
              margins = c(10,5))
    dev.off()
  }else{
    heatmap.2(as.matrix(mvg), Colv = NA, scale="row",
              col =coul, hclustfun = hclust,trace="none",
              dendrogram= "row",ColSideColors=my_col,
              margins = c(10,5))
  }
}

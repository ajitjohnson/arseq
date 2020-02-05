#' @title Heatmap of the differentially expressed genes
#' @description Generating a heatmap of the differentially expressed genes
#' @param deg Differtially expressed genes calculated using DESeq2
#' @param data Normalized count's matrix
#' @param dds DeSes2 object
#' @param save.plot Logical. Argument to state if the plot needs to be saved in disk. Default: FALSE
#' @param save.dir Location to save the plot when 'save.plot=TRUE'. Default: Working Directory.
#' @return Heatmap of the differentially expressed genes.
#' @import RColorBrewer
#' @import grDevices
#' @import gplots

arseq.deg.plot <- function(gene_count,deg, data, dds, save.plot=FALSE,save.dir=getwd()){
  # Subset the significant genes from deg and subset the rows of n_data based on those genes and columns based on goi
  deg_heatmap <- data.frame(deg[which(deg$padj <= 0.05), ])
  # Identify the top 500 genes for the heatmap
  if (nrow(deg_heatmap) > gene_count){
    deg_heatmap <- deg_heatmap[order(deg_heatmap$padj),][1:gene_count,]
  }
  # Order the dataframe based on fold change
  deg_heatmap <- deg_heatmap[order(deg_heatmap$log2FoldChange),]
  # Get the matrix
  deg_heatmap <- data[row.names(deg_heatmap), ]
  # Create the grouping for coloring
  my_group = factor(dds$arseq.group)
  my_col = suppressWarnings(brewer.pal(length(levels(factor(dds$arseq.group))), "Set1")[my_group])
  coul = colorRampPalette(brewer.pal(9,"OrRd"))(10)
  # Save the heatmap
  if (isTRUE(save.plot)){
    print("Generating heatmap of the diferentially expressed genes")
    pdf(paste(save.dir,"DEG Heatmap.pdf",sep = ""),width=20,height=20,paper='special')
    heatmap.2(as.matrix(deg_heatmap), Rowv = NA, Colv = NA, scale="row",
              col =coul, trace="none",main="Top Differentially Expressed Genes",
              dendrogram= "none",ColSideColors=my_col,cexCol = 2,
              margins = c(10,8))
    dev.off()
  }else{
    print("Generating heatmap of the diferentially expressed genes")
    heatmap.2(as.matrix(deg_heatmap), Rowv = NA, Colv = NA, scale="row",
              col =coul, trace="none",main="Top Differentially Expressed Genes",
              dendrogram= "none",ColSideColors=my_col,cexCol = 2,
              margins = c(10,8))
  }

}

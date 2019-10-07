#' @title Principal Component Analysis
#' @description Plots any two specified principal components. By default it plots the first two principal components.
#' @param dds DESeq2 object
#' @param intgroup interesting groups: a character vector of names in colData(x) to use for grouping.
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @param pc.a numeric: Defines the first principal component. Default:1
#' @param pc.b numeric: Defines the second principal component.Default:2
#' @param returnData should the function only return the data.frame of 'pc.a' and 'pc.b' with intgroup covariates for custom plotting (default is FALSE)
#' @return ggplot of the two selected principal components
#' @import DESeq2
#' @import ggrepel
#' @import ggplot2
#' @importFrom stats prcomp
#' @examples
#' pca.plot <- arseq.pca.plot (example_dds, intgroup="treatment")
#' @export

arseq.pca.plot = function(dds, intgroup="arseq.group", ntop=500, pc.a= 1, pc.b = 2, returnData=FALSE){
  print("Performing a PCA analysis")
  # Normalize data
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  # calculate the variance for each gene
  rv <- rowVars(vsd@assays$data[[1]])

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(vsd@assays$data[[1]][select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(vsd@colData))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(vsd@colData[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    vsd@colData[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(pca$x[,pc.a], pca$x[,pc.b], group=group, intgroup.df, name=colnames(vsd))
  names(d)[1] <- paste("PC",pc.a,sep = "")
  names(d)[2] <- paste("PC",pc.b,sep = "")


  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.a:pc.b]
    return(d)
  }

  ggplot(data=d, aes_string(x=paste("PC",pc.a,sep = ""), y=paste("PC",pc.b,sep = ""), color="group")) + geom_point(size=3) +
    xlab(paste0(paste("PC",pc.a,sep = ""),": ",round(percentVar[pc.a] * 100),"% variance")) +
    ylab(paste0(paste("PC",pc.b,sep = ""),": ",round(percentVar[pc.b] * 100),"% variance")) +
    theme_classic()+
    geom_text_repel(aes(label = .data$name),size = 3) +
    coord_fixed() + ggtitle("Principal component analysis (PCA) Plot")+
    theme(plot.title = element_text(hjust = 0.5))
}

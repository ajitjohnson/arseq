#' @title Pheatmap plotting function
#' @description The function saves the pheatmap as a pdf.
#' @param x Expression matrix or any numeric matrix to be visualized as a heatmap.
#' @param filename The pdf file name.
#' @param width Width of the pdf file.
#' @param height Height of the pdf file.
#' @param save.dir Location to save the heatmap
#' @return Does not return anything. The function saves a pdf file to the working directory.
#' @import grid

arseq.pheatmap <- function(x, filename, width=7, height=7, save.dir=getwd()) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(paste(save.dir,"/",filename, sep=""), width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

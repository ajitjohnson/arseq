#' @title Save plots
#' @description Function to save plots as a pdf file in user defined directory
#' @param arseq.plot The plot to be saved
#' @param filename The pdf file name.
#' @param width Width of the pdf file.
#' @param height Height of the pdf file.
#' @param save.dir Location to save the plot
#' @return Does not return an object. Function saves plots in directory
#' @import grDevices
#' @import graphics

arseq.plot.save <- function(arseq.plot, filename, width=7, height=7, save.dir=getwd()){
  pdf(paste(save.dir,"/",filename, sep=""),width=width,height=height,paper='special')
  plot(arseq.plot)
  dev.off()
}

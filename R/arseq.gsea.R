#' @title Gene Set Enrichment Analysis (GSEA)
#' @description Performing GSEA enrichment analysis on a ranked list of genes.
#' @param ranked.list Ranked list of genes for GSEA analysis. Use 'arseq.gsea.preprocess' to generate ranked list of genes from DESeq2 differential gene expression analysis.
#' @param f.name Folder name to save the results when 'save=TRUE'
#' @param fsub.name File name of the results when 'save=TRUE'
#' @param category collection, such as H, C1, C2, C3, C4, C5, C6, C7
#' @param subcategory sub-collection, such as CGP, MIR, BP, etc
#' @param save Logical. Argument to indicate if to save the results in disk (pdf and csv files). Default=FALSE
#' @param save.dir Directory to save the results when 'save=TRUE'. Default= Working Directory.
#' @param custom.gsea User defined gene list to perform GSEA. File need to be supplied as a dataframe with each row as a gene list
#' @return List containing a dataframe and a ggplot object.
#' @import fgsea
#' @import msigdbr
#' @import utils
#' @importFrom dplyr "%>%"
#' @importFrom stats na.omit
#' @import stringr
#' @examples
#' ranked.list <- arseq.gsea.preprocess (deg=example_deg)
#' gsea.output <- arseq.gsea (ranked.list, category = "H")
#' @export

arseq.gsea <- function(ranked.list, custom.gsea=NULL, f.name="GeneSet", fsub.name="enrichment", category= NULL, subcategory= NULL, save=FALSE, save.dir= getwd()){

  # Get the geneset of interest
  if (is.null(custom.gsea)){
    geneset = msigdbr(species = "Homo sapiens", category = category, subcategory=subcategory)
    geneset = geneset %>% split(x = .$gene_symbol, f = .$gs_name)
  }

  # Create a geneset for custom GSEA analysis
  if (!is.null(custom.gsea)){
    if (is.data.frame(custom.gsea)){
      geneset <- as.list(as.data.frame(t(custom.gsea)))
      geneset <- lapply(custom.gsea, as.character)
      geneset <- lapply(custom.gsea, unique)
      custom.gsea <- geneset
    }
    if (is.list(custom.gsea)){
      print("Performing GSEA Custom GeneSet")
    } else{print("Please provide custom geneset in the right format")}
  }

  # Table
  gsea <- suppressWarnings(fgsea(pathways = geneset, stats = ranked.list, nperm=1000))
  fgsea_result<- gsea[with(gsea, order(padj, -NES)), ]
  # Leading edge conversion
  leadingEdge  <- do.call(rbind.data.frame, fgsea_result$leadingEdge)
  colnames(leadingEdge) <- seq(1:ncol(leadingEdge))
  fgsea_table <- as(fgsea_result, "data.frame")[,-ncol(as(fgsea_result, "data.frame"))]
  fgsea_table <- cbind(fgsea_table,leadingEdge)
  # Plot
  fgsea_result <- fgsea_result[order(fgsea_result$NES),]
  fgsea_result <- na.omit(fgsea_result)
  fgsea_sorted <- rbind(fgsea_result[1:10,], fgsea_result[(nrow(fgsea_result)-9):nrow(fgsea_result),])
  fgsea_sorted <- fgsea_sorted[order(fgsea_sorted$NES),]
  fgsea_sorted$pathway <- factor(fgsea_sorted$pathway , levels = fgsea_sorted$pathway )
  gsea.plot <- ggplot(data=fgsea_sorted, aes(x=.data$pathway, y=.data$NES)) +
    geom_col(aes(fill=.data$padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Top 10 most and least enriched sets") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
    theme_bw()

  # Save the results to file
  if (isTRUE(save)){
    # Create a folder to save results
    suppressWarnings(dir.create(paste(save.dir,"/GSEA analysis/",sep = "")))
    suppressWarnings(dir.create(paste(save.dir,"/GSEA analysis/",f.name,sep = "")))

    # Save plot
    pdf(paste(save.dir,"/GSEA analysis/",f.name,"/",fsub.name,".pdf",sep = ""),width=10,height=12,paper='special')
    plot(gsea.plot)
    dev.off()
    # Save table
    write.csv(fgsea_table, file = paste(save.dir,"/GSEA analysis/",f.name,"/",fsub.name,".csv",sep = ""))
  }

  # Result list
  gsea.res <- list()
  gsea.res$table <- fgsea_table
  gsea.res$plot <- gsea.plot

  # Return the result
  return(gsea.res)
}

#' @title Bulk Gene Set Enrichment Analysis (GSEA)
#' @description Performing GSEA enrichment analysis on a ranked list of genes for all msigDB.
#' @param ranked.list Ranked list of genes for GSEA analysis. Use 'arseq.gsea.preprocess' to generate ranked list of genes from DESeq2 differential gene expression analysis.
#' @param save Logical. Argument to indicate if to save the results in disk (pdf and csv files). Default=FALSE
#' @param save.dir Directory to save the results when 'save=TRUE'. Default= Working Directory.
#' @param custom.gsea User defined gene list to perform GSEA. File need to be supplied as a dataframe with each row as a gene list.
#' @return List containing all the GSEA dataframe and ggplot objects.
#' @examples
#' \dontrun{
#' ranked.list <- arseq.gsea.preprocess (deg=example_deg)
#' gsea.output <- arseq.gsea.runall (ranked.list)
#' }
#' @export

arseq.gsea.runall <- function(ranked.list, save=FALSE, custom.gsea=NULL, save.dir= getwd()){

  # Create an empty list to hold all results
  gsea.results <- list()

  # Custom Geneset
  if (!is.null(custom.gsea)){
    print("Performing GSEA with Custom Gene Set")
    gsea.results$customgeneset <- arseq.gsea (ranked.list=ranked.list, custom.gsea=custom.gsea,
                                              f.name="Custom gene sets",fsub.name="custom_geneset",
                                              save=save, save.dir=save.dir)
  }

  # H1
  print("Performing GSEA with Hallmark gene sets")
  gsea.results$h1 <- arseq.gsea (ranked.list=ranked.list,
                                     f.name="Hallmark gene sets",fsub.name="H_geneset",
                                     category = "H",
                                     save=save, save.dir=save.dir)
  # C1
  print("Performing GSEA with C1: Positional gene sets")
  gsea.results$C1 <- arseq.gsea (ranked.list=ranked.list,
                                     f.name="Positional gene sets",fsub.name="C1_geneset",
                                     category = "C1",
                                     save=save, save.dir=save.dir)
  # C2
  print("Performing GSEA with C2: Curated gene sets")
  gsea.results$C1.CGP <- arseq.gsea (ranked.list=ranked.list,
                                     f.name="Curated gene sets",fsub.name="CGP",
                                     category = "C2",subcategory="CGP",
                                     save=save, save.dir=save.dir)

  gsea.results$C1.CP <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="Curated gene sets",fsub.name="CP",
                                    category = "C2",subcategory="CP",
                                    save=save, save.dir=save.dir)

  gsea.results$C1.CP_BIOCARTA <- arseq.gsea (ranked.list=ranked.list,
                                             f.name="Curated gene sets",fsub.name="CP_BIOCARTA",
                                             category = "C2",subcategory="CP:BIOCARTA",
                                             save=save, save.dir=save.dir)

  gsea.results$C1.CP_KEGG <- arseq.gsea (ranked.list=ranked.list,
                                         f.name="Curated gene sets",fsub.name="CP_KEGG",
                                         category = "C2",subcategory="CP:KEGG",
                                         save=save, save.dir=save.dir)

  gsea.results$C1.CP_REACTOME <- arseq.gsea (ranked.list=ranked.list,
                                             f.name="Curated gene sets",fsub.name="CP_REACTOME",
                                             category = "C2",subcategory="CP:REACTOME",
                                             save=save, save.dir=save.dir)
  # C3
  print("Performing GSEA with C3: Motif gene sets")
  gsea.results$C3.MIR <- arseq.gsea (ranked.list=ranked.list,
                                     f.name="Motif gene sets",fsub.name="MIR",
                                     category = "C3",subcategory="MIR:MIRDB",
                                     save=save, save.dir=save.dir)

  gsea.results$C3.TFT <- arseq.gsea (ranked.list=ranked.list,
                                     f.name="Motif gene sets",fsub.name="TFT",
                                     category = "C3",subcategory="TFT:GTRD",
                                     save=save, save.dir=save.dir)
  # C4
  print("Performing GSEA with C4: Computational gene sets")
  gsea.results$C4.CGN <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="Computational gene sets",fsub.name="CGN",
                                    category = "C4",subcategory="CGN",
                                     save=save, save.dir=save.dir)

  gsea.results$C4.CM <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="Computational gene sets",fsub.name="CM",
                                    category = "C4",subcategory="CM",
                                    save=save, save.dir=save.dir)
  # C5
  print("Performing GSEA with C5: GO gene sets")
  gsea.results$C5.BP <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="GO gene sets",fsub.name="BP",
                                    category = "C5",subcategory="BP",
                                    save=save, save.dir=save.dir)

  gsea.results$C5.CC <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="GO gene sets",fsub.name="CC",
                                    category = "C5",subcategory="CC",
                                    save=save, save.dir=save.dir)

  gsea.results$C5.MF <- arseq.gsea (ranked.list=ranked.list,
                                    f.name="GO gene sets",fsub.name="MF",
                                    category = "C5",subcategory="MF",
                                    save=save, save.dir=save.dir)
  # C6
  print("Performing GSEA with C6: Oncogenic signatures")
  gsea.results$C6 <- arseq.gsea (ranked.list=ranked.list,
                                 f.name="Oncogenic signatures",fsub.name="C6_geneset",
                                 category = "C6",
                                 save=save, save.dir=save.dir)
  # C7
  print("Performing GSEA with C7: Immunologic signatures")
  gsea.results$C7 <- arseq.gsea (ranked.list=ranked.list,
                                 f.name="Immunologic signatures",fsub.name="C7_geneset",
                                 category = "C7",
                                 save=save, save.dir=save.dir)

  # Return the final result
  return(gsea.results)

}

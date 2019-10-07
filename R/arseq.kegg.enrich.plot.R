#' @title Gene Ontology (GO) term enrichment
#' @description Plotting the top 5 enriched term from each category
#' @param kegg.enrich Enrichment dataframe. Output of the function 'arseq.kegg.enrich'
#' @param save.dir User defined directory to save the plots in. Default is working directory.
#' @param foldchanges Foldchange information to overlay on the pathway plots.
#' @param pathway.plots Logical.Parameter to indicate if the pathways figures need to be drawn.
#' @return GO enrichment plot
#' @import ggplot2
#' @importFrom dplyr tbl_df filter row_number "%>%"
#' @import pathview
#' @import stringr
#' @examples
#' \dontrun{
#' kegg.enrich <- arseq.kegg.enrich (example_deg)
#' kegg.plot <- arseq.kegg.enrich.plot (kegg.enrich, foldchanges=kegg.enrich$foldchanges)
#' }
#' @export

arseq.kegg.enrich.plot <- function(kegg.enrich,foldchanges,save.dir=getwd(), pathway.plots=TRUE){

  if (isTRUE(pathway.plots)){
    # Individual pathway figures for the top 5 pathways
    keggrespathways_up = data.frame(id=rownames(kegg.enrich$up_regulated_pathways), kegg.enrich$up_regulated_pathways) %>% tbl_df() %>%
      dplyr::filter(.data$q.val<=0.05) %>% dplyr::filter(row_number()<=5) %>% .data$id %>% as.character()

    keggrespathways_down = data.frame(id=rownames(kegg.enrich$down_regulated_pathways), kegg.enrich$down_regulated_pathways) %>% tbl_df() %>%
      dplyr::filter(.data$q.val<=0.05) %>% dplyr::filter(row_number()<=5) %>% .data$id %>% as.character()
    # collate the pathways
    keggrespathways = c(keggrespathways_up, keggrespathways_down)
    # Get the IDs.
    keggresids = substr(keggrespathways, start=1, stop=8)
    # Define plotting function for applying later
    plot_pathway <- function(dir, gene.data, pathway.id, species) {
      old_wd <- getwd()
      on.exit(setwd(old_wd))
      setwd(dir)
      pathview(gene.data=foldchanges, pathway.id=pathway.id, species="hsa")
    }

    if (length(keggresids) >0 ){
      suppressWarnings(dir.create(paste(save.dir,"/top 5 pathways (up and down)",sep = "")))
      tmp = sapply(keggresids, function (x) plot_pathway(dir=paste(save.dir,"/top 5 pathways (up and down)/",sep=""),
                                                         gene.data=foldchanges, pathway.id=x, species="hsa"))
    }
  }
  # Prepare the KEGG dataframe for generating a figure
  print("Generating a figure for the top 10 enriched terms")
  kegg_top <- kegg.enrich$up_regulated_pathways[1:10,]
  kegg_top$adjp.val <- -log(kegg_top$q.val)
  kegg_top$direction <- "Up regulated"
  kegg_bottom <- kegg.enrich$down_regulated_pathways[1:10,]
  kegg_bottom$adjp.val <- -log(kegg_bottom$q.val)
  kegg_bottom$direction <- "Down regulated"
  kegg.plot <- rbind(kegg_top,kegg_bottom)

  # Plot the kegg.plot data
  kegg.plot$Name <- as.character(row.names(kegg.plot))
  kegg.plot$Name <- substring(kegg.plot$Name, 10)
  kegg.plot$Name <- factor(kegg.plot$Name, levels = kegg.plot$Name) # to avoid re-ordering

  ggplot(data=kegg.plot, aes(x=.data$Name, y=.data$adjp.val,fill=.data$direction)) +
    geom_bar(stat="identity")+ geom_hline(yintercept = 3)+
    coord_flip()+ theme_bw()+ theme(legend.position="bottom",legend.title = element_blank(), axis.text=element_text(size=12)) +
    ggtitle("Top 10 KEGG Pathways")+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 25))+
    labs(y = "-log(P value)")
 }

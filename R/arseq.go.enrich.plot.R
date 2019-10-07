#' @title Gene Ontology (GO) term enrichment
#' @description Plotting the top 5 enriched term from each category
#' @param go.enrich Enrichment dataframe. Output of the function 'arseq.go.enrich'
#' @return GO enrichment plot
#' @import ggplot2
#' @import gridExtra
#' @examples
#' \dontrun{
#' go.enrich <- arseq.go.enrich (example_deg)
#' go.plot <- arseq.go.enrich.plot (go.enrich)
#' }
#' @export

arseq.go.enrich.plot <- function(go.enrich){
  GO.plot <- data.frame()
  # Prepare the GO terms
  for (i in levels(as.factor(go.enrich$ontology))){
    GO.wall_subset <- go.enrich[go.enrich$ontology %in% i,]
    GO.wall_top <- GO.wall_subset[rownames(GO.wall_subset)[order(GO.wall_subset$over_represented_pvalue, decreasing=FALSE)][1:5],]
    GO.wall_top$direction <- "over represented"
    GO.wall_top$p.val <- -log(GO.wall_top$over_represented_pvalue)
    GO.wall_bottom <- GO.wall_subset[rownames(GO.wall_subset)[order(GO.wall_subset$under_represented_pvalue, decreasing=FALSE)][1:5],]
    GO.wall_bottom$direction <- "under represented"
    GO.wall_bottom$p.val <- -log(GO.wall_bottom$under_represented_pvalue)
    GO_subset <- rbind(GO.wall_top,GO.wall_bottom)
    GO.plot <- rbind(GO.plot,GO_subset)
  }
  # Plot the GO.plot data
  print("Generating a figure for the top 5 enriched terms in each category")
  GO.plot$term <- factor(GO.plot$term, levels = GO.plot$term) # to avoid re-ordering
  for (i in levels(as.factor(GO.plot$ontology))){
    GO.plot_subset <- GO.plot[GO.plot$ontology %in% i,]
    p<-ggplot(data=GO.plot_subset, aes(x=.data$term, y=.data$p.val,fill=.data$direction)) +
      geom_bar(stat="identity")+ geom_hline(yintercept = 3)+
      coord_flip()+ theme_bw()+ theme(legend.position="bottom",legend.title = element_blank(), axis.text=element_text(size=12)) +
      ggtitle(paste(i, "- Top 5"))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 25))+
      labs(y = "-log(P value)")
    my_plots[[i]] <- p
  }
  grid.arrange(grobs = my_plots, ncol = 1, padding = unit(2, "line"), top = textGrob("Top 5 enriched terms",gp=gpar(fontsize=20,font=3)))
}

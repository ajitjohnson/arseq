#' @title Automated RNASeq Analysis Pipeline
#' @description An easy to use pipeline for analysing RNA Seq data. The package currently supports the following analysis- Differential gene expression analysis using DESeq2, Calculate the most variable genes, PCA analysis, GO enrichment of the differentially expressed genes, KEGG pathway enrichment of the differentially expressed genes, GSEA analysis.
#' @param data Un-normalized counts matrix (please note that you should NOT pass in normalized data).  The counts' table should contain unique gene names as the first column. ENSEMBL ID's are also allowed but no other form of ID's are currently supported. Check example- head(example_data): \code{\link{example_data}}.
#' @param meta A CSV file with information regarding the samples. It is absolutely critical that the columns of the counts' matrix and the rows of the metadata are in the same order. The function will not make guesses as to which column of the count matrix belongs to which row of the metadata, these must be provided in a consistent order. Check example- head(example_meta): \code{\link{example_meta}}.
#' @param design The formula that expresses how the counts for each gene depend on your metadata (used to calculate the necessary data for differential gene expression analysis). Check DESeq2 documentation for designing the formula. In general, you pass in a column name (e.g. treatment) of your metadata file or a combination of column names (e.g. treatment + cell_type).
#' @param contrast Information regarding the groups between which you would like to perform differential gene expression analysis. It could be between two groups or between multiple groups and needs to follow the following format: contrast = list(A = c(" "), B= c(" ")). If you are comparing two groups (e.g. control vs treatment), the constrast argument should look like the following: contrast = list(A = c("control"), B= c("treatment")). In situations where you have multiple groups to compare- (e.g. control vs treatment1 and treatment2), you should do the following- contrast = list(A = c("control"), B= c("treatment1", "treatment2")).
#' @param variable.genes numeric: The number of most variable genes to be identified. By default, the program identifies the top 1000 most variable genes.
#' @param general.stats TRUE/FALSE. When passed in TRUE, the program would run the general stat modules on the entire dataset. If you are planning to perform multiple comparisons using the contrast argument, run  general.stats = TRUE for the first time and then change it to general.stats = FALSE for the subsequent comparisons to speed up the analysis.
#' @return The program does not return anything to the R environment. All analysis results will be stored under a folder named "arseq" within your working directory.
#' @import DESeq2
#' @import goseq
#' @import plyr
#' @importFrom dplyr tbl_df filter row_number "%>%"
#' @import ggplot2
#' @import gridExtra
#' @import stringr
#' @import KEGGREST
#' @import fgsea
#' @import msigdbr
#' @importFrom biomaRt useEnsembl getBM
#' @import RColorBrewer
#' @import pheatmap
#' @import viridis
#' @import PoiClaClu
#' @import ggrepel
#' @import grid
#' @import gplots
#' @import EnhancedVolcano
#' @import pathview
#' @import gage
#' @import utils
#' @import grDevices
#' @import graphics
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats aggregate as.formula cmdscale complete.cases dist hclust na.omit prcomp
#' @examples
#' \dontrun{
#' contrast = list(A = c("control"), B= c("drug_A"))
#' arseq (data = example_data,meta = example_meta, design = "treatment", contrast = contrast)
#' }
#' @export

arseq <- function(data,meta,design, contrast, general.stats= TRUE, variable.genes=1000){

  # Resolve the contrasts and design
  goi <-  c(paste(unlist(contrast[[1]]), collapse="_"),paste(unlist(contrast[[2]]), collapse="_"))
  groups <- c(unlist(contrast[[1]]),contrast[[2]])

  # Figure out the contrast column programatically
  colum.names <- c()
  for (i in 1: length(groups)){
    colum.names <- c(colum.names,colnames(meta[, which(meta == groups[i], arr.ind=T)[1,][[2]], drop=FALSE]))
  }
  if (all(colum.names == colum.names[1])){
    contrast.by = colum.names[1]
    # Remname the column with contrast information to arseq.group
    colnames(meta)[which(colnames(meta) == contrast.by)] <- "arseq.group"
  }else{stop('All of your contrast groups for Differntial expression analysis need to be listed under the same column within your metadata file')}

  # resolving design
  design <- gsub(contrast.by, "arseq.group", design)
  design <- paste("~",design,sep = "")
  design <- as.formula(design)

  # Creating a folder to save results
  print("Creating an 'ARSeq' folder to store your analysis results")
  suppressWarnings(dir.create("ARSeq"))
  suppressWarnings(dir.create("ARSeq/General stats"))
  suppressWarnings(dir.create(paste("ARSeq/",goi[1], " vs ", goi[2],sep = "")))
  location <- paste("ARSeq/",goi[1], " vs ", goi[2],"/",sep = "")
  suppressWarnings(dir.create(paste(location,"Differential expression/",sep = "")))
  suppressWarnings(dir.create(paste(location,"Sample relation plots",sep = "")))
  suppressWarnings(dir.create(paste(location,"GO enrichment analysis/",sep = "")))
  suppressWarnings(dir.create(paste(location,"KEGG pathway enrichment/",sep = "")))
  suppressWarnings(dir.create(paste(location,"GSEA analysis/",sep = "")))

  # If data is based of ENSEMBL ID convert it to gene names
  if (all(substr(row.names(data),1,3) == "ENS")){
    print("Converting ENSEMBL ID's to gene names")
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)
    data_m <- merge(data, genes, by.x="row.names", by.y= "ensembl_gene_id")[,-1]
    data_m <- data_m[!(data_m$hgnc_symbol==""), ]
    # Collapse to one gene
    data_m$hgnc_symbol <- as.factor(data_m$hgnc_symbol)
    data <- aggregate(.~hgnc_symbol, data=data_m, max)
    rownames(data) <- data[,1]
    data <- data[,-1]}

  # If data has duplicates collapse it into one gene
  # Remove genes that are completely not expressed
  print("Removing genes that are not expressed: if any")
  data <- data[rowSums(data) > 0, ]

  # Generate the DESEq2 data object
  print("Generating the DESeq object")
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = design)
  dds <- DESeq(dds)

  # Get the normalized matrix for heatmaps
  print("Normalizing data and saving the normalized data")
  n_data <- log2(counts(dds, normalized=TRUE)+1)
  write.csv(n_data, file = paste("ARSeq/normalized_data",".csv",sep = ""))

  # Common function
  save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }

  if (general.stats == TRUE){
    # General structure of the dataset
    print("Computing general structure of the entire dataset")
    # Sample to sample distance
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    print("Calculating the Euclidean distance between samples")
    sampleDists <- dist(t(vsd@assays$data[[1]]))
    sampleDistMatrix <- as.matrix( sampleDists )
    colnames(sampleDistMatrix) <- NULL
    # save plot
    euclid_dist <- pheatmap(sampleDistMatrix,
                            clustering_distance_rows = sampleDists,
                            clustering_distance_cols = sampleDists,
                            border_color = NA,
                            color = inferno(50),
                            main = "Euclidean distance between samples")

    save_pheatmap_pdf(euclid_dist, "ARSeq/General stats/Euclidean distance.pdf")

    # Poisson Distance
    print("Calculating the Poisson distance between samples")
    poisd <- PoissonDistance(t(counts(dds)))
    samplePoisDistMatrix <- as.matrix( poisd$dd )
    rownames(samplePoisDistMatrix) <- colnames(dds)
    colnames(samplePoisDistMatrix) <- NULL
    poi_dist <- pheatmap(samplePoisDistMatrix,
                         clustering_distance_rows = poisd$dd,
                         clustering_distance_cols = poisd$dd,
                         border_color = NA,
                         color = inferno(50),
                         main = "Poisson distance between samples")
    save_pheatmap_pdf(poi_dist, "ARSeq/General stats/Poisson distance.pdf")

    #PCA plot
    print("Performing a PCA analysis")
    pca.plot <- plotPCA (vsd,intgroup = "arseq.group") + geom_text_repel(aes(label = name),size = 3)
    pdf("ARSeq/General stats/PCA plot.pdf",width=6,height=6,paper='special')
    plot(pca.plot)
    dev.off()

    # Save the PC's
    pca = prcomp(n_data, center=TRUE, scale=TRUE)
    pca <- data.frame(pca$x)
    pca <- pca[order(-pca[,1]),]
    write.csv(pca, file = "ARSeq/General stats/PCA Eigenvectors.csv")

    # MSD plot
    print("Performing a multidimensional scaling (MDS) analysis")
    mds <- as.data.frame(vsd@colData)  %>% cbind(cmdscale(sampleDistMatrix))
    mds.plot <- ggplot(mds, aes_string(x = mds$`1`, y = mds$`2`, color = "arseq.group" )) +
      geom_point(size = 3) + coord_fixed()+ geom_text_repel(aes(label = rownames(mds)),size = 3)+
      theme(legend.title = element_blank())
    pdf("ARSeq/General stats/MDS plot.pdf",width=6,height=6,paper='special')
    plot(mds.plot)
    dev.off()

    # Most variable genes
    print("Identifying the most variable genes in your dataset")
    data.var <- apply(n_data, 1, stats::var)
    var.genes <- n_data[order(data.var, decreasing = TRUE)[1:variable.genes],]
    write.csv(var.genes, file = "ARSeq/General stats/most variable genes.csv")
    # Heatmap of the most variable genes
    my_group = factor(meta[,"arseq.group"])
    my_col = brewer.pal(length(levels(factor(meta[,"arseq.group"]))), "Set1")[my_group]
    coul = colorRampPalette(brewer.pal(9,"OrRd"))(10)
    pdf("ARSeq/General stats/Heatmap of the most variable genes.pdf",width=6,height=6,paper='special')
    heatmap.2(as.matrix(var.genes), Colv = NA, scale="row",
              col =coul, hclustfun = hclust,trace="none",
              dendrogram= "row",ColSideColors=my_col,
              margins = c(10,5))
    dev.off()
  }

  # Manipulate the dds groups to look at the contrasts of interest
  dds$arseq.group <- mapvalues(dds$arseq.group, from=contrast[[1]], to=rep(paste(unlist(contrast[[1]]), collapse="_"),length(contrast[[1]])))
  dds$arseq.group <- mapvalues(dds$arseq.group, from=contrast[[2]], to=rep(paste(unlist(contrast[[2]]), collapse="_"),length(contrast[[2]])))

  # Subset the normalized matrix and dds based on the groups of interest
  colums_to_subset <- row.names(meta[meta$arseq.group %in% groups,,drop=FALSE]) # Changed
  n_data_goi <- data.frame(n_data[,colums_to_subset])
  #n_data_goi <- data.frame(n_data[,unique(grep(paste(groups,collapse="|"), colnames(n_data), value=TRUE))])
  dds_subset <- dds[ , dds$arseq.group %in% goi ]
  dds_subset$arseq.group <- factor(dds_subset$arseq.group)

  # DEG
  print("Performing Differential gene expression analysis between the  constrast groups")

  dds <- DESeq(dds)
  deg <- results(dds, contrast=c("arseq.group",goi[1],goi[2]))
  deg <- deg[order(deg$padj),]
  # Save the DEG matrix
  print("Saving the differentially expressed genes")
  write.csv(deg, file = paste(location,"Differential expression/",goi[1], " vs ", goi[2],".csv",sep = ""))
  # Subset the significant genes from deg and subset the rows of n_data based on those genes and columns based on goi
  deg_heatmap <- data.frame(deg[which(deg$padj <= 0.05), ])
  # Identify the top 500 genes for the heatmap
  if (nrow(deg_heatmap) > 500){
    deg_heatmap <- deg_heatmap[order(deg_heatmap$padj),][1:500,]
  }
  # Order the dataframe based on fold change
  deg_heatmap <- deg_heatmap[order(deg_heatmap$log2FoldChange),]
  # Get the matrix
  deg_heatmap <- n_data_goi[row.names(deg_heatmap[which(deg_heatmap$padj <= 0.05), ]), ]
  # Create the grouping for coloring
  my_group = factor(dds_subset$arseq.group)
  my_col = suppressWarnings(brewer.pal(length(levels(factor(dds_subset$arseq.group))), "Set1")[my_group])
  coul = colorRampPalette(brewer.pal(9,"OrRd"))(10)
  # Save the heatmap
  print("Generating heatmap of the diferentially expressed genes")
  pdf(paste(location,"Differential expression/","DEG Heatmap- ",goi[1], " vs ", goi[2],".pdf",sep = ""),width=6,height=6,paper='special')
  heatmap.2(as.matrix(deg_heatmap), Rowv = NA, Colv = NA, scale="row",
            col =coul, trace="none",main="Top 500 Differentially Expressed Genes",
            dendrogram= "none",ColSideColors=my_col,
            margins = c(10,5))
  dev.off()

  # Euclidean distance
  print("Calculating the Euclidean distance between constrast groups")

  vsd_subset <- varianceStabilizingTransformation(dds_subset, blind = FALSE)

  sampleDists <- dist(t(vsd_subset@assays$data[[1]]))
  sampleDistMatrix <- as.matrix( sampleDists )
  colnames(sampleDistMatrix) <- NULL
  # save plot
  euclid_dist <- pheatmap(sampleDistMatrix,
                          clustering_distance_rows = sampleDists,
                          clustering_distance_cols = sampleDists,
                          border_color = NA,
                          color = inferno(50),
                          main = "Euclidean distance between samples")
  print("Saving the Euclidean distance plot")
  save_pheatmap_pdf(euclid_dist, paste(location,"Sample relation plots/","Euclidean distance- ",goi[1], " vs ", goi[2],".pdf",sep = ""))


  # Poisson Distance
  print("Calculating the Poisson distance between constrast groups")

  poisd <- PoissonDistance(t(counts(dds_subset)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  rownames(samplePoisDistMatrix) <- colnames(dds_subset)
  colnames(samplePoisDistMatrix) <- NULL
  poi_dist <- pheatmap(samplePoisDistMatrix,
                       clustering_distance_rows = poisd$dd,
                       clustering_distance_cols = poisd$dd,
                       border_color = NA,
                       color = inferno(50),
                       main = "Poisson distance between samples")
  print("Saving the Poisson distance plot")
  save_pheatmap_pdf(poi_dist, paste(location,"Sample relation plots/","Poisson distance- ",goi[1], " vs ", goi[2],".pdf",sep = ""))

  # PCA
  print("Performing Principal Component Analysis (PCA) between the constrast groups")

  pca.plot <- plotPCA (vsd_subset,intgroup = "arseq.group") + geom_text_repel(aes(label = name),size = 3)
  pdf(paste(location,"Sample relation plots/","PCA- ",goi[1], " vs ", goi[2],".pdf",sep = ""),width=6,height=6,paper='special')
  plot(pca.plot)
  dev.off()
  print("Saving the PCA plot")
  # Save the PC's
  pca = prcomp(n_data, center=TRUE, scale=TRUE)
  pca <- data.frame(pca$x)
  pca <- pca[order(-pca[,1]),]
  write.csv(pca, file = paste(location,"Sample relation plots/",goi[1], " vs ", goi[2]," PCA Eigenvectors.csv",sep = ""))

  # MSD plot
  print("Performing multidimensional scaling (MDS) analysis between the constrast groups")

  mds <- data.frame(vsd_subset@colData)  %>% cbind(cmdscale(sampleDistMatrix))
  mds.plot <- ggplot(mds, aes_string(x = mds$`1`, y = mds$`2`, color = "arseq.group" )) +
    geom_point(size = 3) + coord_fixed()+ geom_text_repel(aes(label = rownames(mds)),size = 3)+
    theme(legend.title = element_blank())
  print("Saving the MDS plot")
  pdf(paste(location,"Sample relation plots/","MDS- ",goi[1], " vs ", goi[2],".pdf",sep = ""),width=6,height=6,paper='special')
  plot(mds.plot)
  dev.off()

  # Volcano plot
  print("Generating a volcano plot between the constrast groups")

  deg_volcano <- data.frame(deg)[complete.cases(data.frame(deg)),]
  volcano.plot <- EnhancedVolcano(deg_volcano,
                  lab = rownames(deg_volcano),
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  pCutoff = 0.01,
                  FCcutoff = 2.0,
                  transcriptLabSize = 4.0,
                  transcriptPointSize = 2.0,
                  colAlpha = 1,
                  legend=c('NS','Log2 FC','Adjusted p-value',
                           'Adjusted p-value (<0.01) & Log2 FC'),
                  legendPosition = 'bottom',
                  legendIconSize = 3.0)
  print("Saving the Volcano plot")
  pdf(paste(location,"Sample relation plots/","Volcano Plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),width=8,height=8,paper='special')
  plot(volcano.plot)
  dev.off()

  # GO enrichment analysis for the DEG's
  print("Performing GO enrichment analysis between the constrast groups")

  gene.vector=as.integer(row.names(data.frame(deg))%in%row.names(deg_heatmap))
  names(gene.vector)=row.names(data.frame(deg))
  # Prepare for goseq (gene length data for the DEG's)
  pwf=suppressWarnings(nullp(gene.vector,"hg19","geneSymbol"))
  # Using the Wallenius approximation
  GO.wall=goseq(pwf,"hg19","geneSymbol")
  # save the goterms
  write.csv(GO.wall, file = paste(location,"GO enrichment analysis/","GO Enrichment- ",goi[1], " vs ", goi[2],".csv", sep = ""))

  # Prepare the GO dataframe for generating a figure
  GO.plot <- data.frame()
  for (i in levels(as.factor(GO.wall$ontology))){
    GO.wall_subset <- GO.wall[GO.wall$ontology %in% i,]
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
    p<-ggplot(data=GO.plot_subset, aes(x=term, y=p.val,fill=direction)) +
      geom_bar(stat="identity")+ geom_hline(yintercept = 3)+
      coord_flip()+ theme_bw()+ theme(legend.position="bottom",legend.title = element_blank(), axis.text=element_text(size=12)) +
      ggtitle(paste(i, "- Top 5"))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 25))+
      labs(y = "-log(P value)")
    my_plots[[i]] <- p
  }
  pdf(paste(location,"GO enrichment analysis/","GO Enrichment- ",goi[1], " vs ", goi[2],".pdf",sep = ""),width=6,height=20,paper='special')
  grid.arrange(grobs = my_plots, ncol = 1, padding = unit(2, "line"), top = textGrob("Top 5 enriched terms",gp=gpar(fontsize=20,font=3)))
  dev.off()

  # Kegg Pathway Analysis for the DEG's
  print("Performing KEGG pathway analysis between the constrast groups")

  # Find Entrez ID
  deg$entrez = mapIds(org.Hs.eg.db,
                      keys=row.names(deg),
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  # Intialize the gage function
  foldchanges = deg$log2FoldChange
  names(foldchanges) = deg$entrez
  # Get the results
  kegg.sets.hs <- kegg.sets.hs
  keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
  up_regulated_pathways <- data.frame(keggres[["greater"]])
  down_regulated_pathways <- data.frame(keggres[["less"]])
  # write the reults
  write.csv(up_regulated_pathways, file = paste(location,"KEGG pathway enrichment/","KEGG Up Regulated Pathways- ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(down_regulated_pathways, file = paste(location,"KEGG pathway enrichment/","KEGG Down Regulated Pathways- ",goi[1], " vs ", goi[2],".csv",sep = ""))

  # Figure
  # Individual pathway figures for the top 5 pathways
  # Get the pathways
  keggrespathways_up = data.frame(id=rownames(keggres$greater), keggres$greater) %>% tbl_df() %>%
    dplyr::filter(q.val<=0.05) %>% dplyr::filter(row_number()<=5) %>% .$id %>% as.character()

  keggrespathways_down = data.frame(id=rownames(keggres$less), keggres$less) %>% tbl_df() %>%
    dplyr::filter(q.val<=0.05) %>% dplyr::filter(row_number()<=5) %>% .$id %>% as.character()
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
    suppressWarnings(dir.create(paste(location,"KEGG pathway enrichment/top 5 pathways (up and down)",sep = "")))
    tmp = sapply(keggresids, function (x) plot_pathway(dir =paste(location,"KEGG pathway enrichment/top 5 pathways (up and down)/",sep=""),
                                                       gene.data=foldchanges, pathway.id=x, species="hsa"))
  }
  # Prepare the KEGG dataframe for generating a figure
  print("Generating a figure for the top 10 enriched terms")
  kegg_top <- up_regulated_pathways[1:10,]
  kegg_top$adjp.val <- -log(kegg_top$q.val)
  kegg_top$direction <- "Up regulated"
  kegg_bottom <- down_regulated_pathways[1:10,]
  kegg_bottom$adjp.val <- -log(kegg_bottom$q.val)
  kegg_bottom$direction <- "Down regulated"
  kegg.plot <- rbind(kegg_top,kegg_bottom)

  # Plot the kegg.plot data
  kegg.plot$Name <- as.character(row.names(kegg.plot))
  kegg.plot$Name <- substring(kegg.plot$Name, 10)
  kegg.plot$Name <- factor(kegg.plot$Name, levels = kegg.plot$Name) # to avoid re-ordering

  kegg_plot<-ggplot(data=kegg.plot, aes(x=Name, y=adjp.val,fill=direction)) +
    geom_bar(stat="identity")+ geom_hline(yintercept = 3)+
    coord_flip()+ theme_bw()+ theme(legend.position="bottom",legend.title = element_blank(), axis.text=element_text(size=12)) +
    ggtitle("Top 10 KEGG Pathways")+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 25))+
    labs(y = "-log(P value)")
  # Plot
  pdf(paste(location,"KEGG pathway enrichment/","KEGG Pathway Analysis- ",goi[1], " vs ", goi[2],".pdf",sep=""),width=8,height=12,paper='special')
  plot(kegg_plot)
  dev.off()

  # GSEA Analysis
  print("GSEA analysis between the constrast groups")
  print("Preparing the GSEA components")
  # Create a ranked list of genes
  ranked_list <- data.frame(deg[!is.na(deg$padj),])
  ranked_list$score <- -log(ranked_list$padj)*sign(ranked_list$log2FoldChange)
  ranked_list <- ranked_list[order(-ranked_list$score),]
  r_list <- ranked_list$score
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  r_list = r_list[is.finite(r_list)]

  # Download the necessary genesets
  geneset_h = msigdbr(species = "Homo sapiens", category = "H")
  geneset_c1 = msigdbr(species = "Homo sapiens", category = "C1")
  geneset_c2 = msigdbr(species = "Homo sapiens", category = "C2")
  geneset_c3 = msigdbr(species = "Homo sapiens", category = "C3")
  geneset_c4 = msigdbr(species = "Homo sapiens", category = "C4")
  geneset_c5 = msigdbr(species = "Homo sapiens", category = "C5")
  geneset_c6 = msigdbr(species = "Homo sapiens", category = "C6")
  geneset_c7 = msigdbr(species = "Homo sapiens", category = "C7")

  # Create a plot of the top pathways
  fgsea_plot <- function(fgsea_result,geneset,r_list,analysis){
    fgsea_result <- fgsea_result[order(fgsea_result$NES),]
    fgsea_result <- na.omit(fgsea_result)
    fgsea_sorted <- rbind(fgsea_result[1:10,], fgsea_result[(nrow(fgsea_result)-9):nrow(fgsea_result),])
    fgsea_sorted <- fgsea_sorted[order(fgsea_sorted$NES),]
    fgsea_sorted$pathway <- factor(fgsea_sorted$pathway , levels = fgsea_sorted$pathway )
    gsea.plot <- ggplot(data=fgsea_sorted, aes(x=pathway, y=NES)) +
      geom_col(aes(fill=padj<0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Top 10 most and least enriched sets") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
      theme_minimal()
      #plot
      pdf(paste(location,"GSEA analysis/","GSEA- ",analysis, " ", goi[1], " vs ", goi[2],".pdf",sep = ""),width=10,height=12,paper='special')
      plot(gsea.plot)
      dev.off()
    }
  # H1
  print("Performing GSEA with hallmark gene sets")
  geneset = geneset_h %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_h1<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_h1<- fgsea_h1[with(fgsea_h1, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_h1,geneset,r_list,analysis="geneset_H")
  # C1
  print("Performing GSEA with C1: positional gene sets")
  geneset = geneset_c1 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c1<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c1<- fgsea_c1[with(fgsea_c1, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c1,geneset,r_list,analysis="geneset_C1")
  # C2
  print("Performing GSEA with C2: curated gene sets")
  geneset = geneset_c2 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c2<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c2<- fgsea_c2[with(fgsea_c2, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c2,geneset,r_list,analysis="geneset_C2")
  # C3
  print("Performing GSEA with C3: positional gene sets")
  geneset = geneset_c3 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c3<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c3<- fgsea_c3[with(fgsea_c3, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c3,geneset,r_list,analysis="geneset_C3")
  # C4
  print("Performing GSEA with C4: computational gene sets")
  geneset = geneset_c4 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c4<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c4<- fgsea_c4[with(fgsea_c4, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c4,geneset,r_list,analysis="geneset_C4")
  # C5
  print("Performing GSEA with C5: GO gene sets")
  geneset = geneset_c5 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c5<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c5<- fgsea_c5[with(fgsea_c5, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c5,geneset,r_list,analysis="geneset_C5")
  # C6
  print("Performing GSEA with C6: oncogenic signatures")
  geneset = geneset_c6 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c6<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c6<- fgsea_c6[with(fgsea_c6, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c6,geneset,r_list,analysis="geneset_C6")
  # C7
  print("Performing GSEA with C7: immunologic signatures")
  geneset = geneset_c7 %>% split(x = .$gene_symbol, f = .$gs_name)
  fgsea_c7<- suppressWarnings(fgsea(pathways = geneset, stats = r_list, nperm=1000))
  fgsea_c7<- fgsea_c7[with(fgsea_c7, order(padj, -NES)), ]
  fgsea_plot (fgsea_result=fgsea_c7,geneset,r_list,analysis="geneset_C7")

  # Save all results
  write.csv(as(fgsea_h1, "data.frame")[,-ncol(as(fgsea_h1, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- H1 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c1, "data.frame")[,-ncol(as(fgsea_c1, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C1 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c2, "data.frame")[,-ncol(as(fgsea_c2, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C2 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c3, "data.frame")[,-ncol(as(fgsea_c3, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C3 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c4, "data.frame")[,-ncol(as(fgsea_c4, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C4 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c5, "data.frame")[,-ncol(as(fgsea_c5, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C5 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c6, "data.frame")[,-ncol(as(fgsea_c6, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C6 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))
  write.csv(as(fgsea_c7, "data.frame")[,-ncol(as(fgsea_c7, "data.frame"))], file = paste(location,"GSEA analysis/","GSEA Analysis- C7 geneset ",goi[1], " vs ", goi[2],".csv",sep = ""))

  print(paste("Well done- your analysis is now complete. Head over to [[", getwd(), "]] to view your results"))
}

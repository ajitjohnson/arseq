#' @title Automated RNASeq Analysis Pipeline
#' @description An easy to use pipeline for analysing RNA Seq data. The package currently supports the following analysis- Differential gene expression analysis using DESeq2, Calculate the most variable genes, PCA analysis, GO enrichment of the differentially expressed genes, KEGG pathway enrichment of the differentially expressed genes, GSEA analysis.
#' @param data Un-normalized counts matrix (please note that you should NOT pass in normalized data).  The counts' table should contain unique gene names as the first column. ENSEMBL ID's are also allowed but no other form of ID's are currently supported. Check example- head(example_data): \code{\link{example_data}}.
#' @param meta A CSV file with information regarding the samples. It is absolutely critical that the columns of the counts' matrix and the rows of the metadata are in the same order. The function will not make guesses as to which column of the count matrix belongs to which row of the metadata, these must be provided in a consistent order. Check example- head(example_meta): \code{\link{example_meta}}.
#' @param design The formula that expresses how the counts for each gene depend on your metadata (used to calculate the necessary data for differential gene expression analysis). Check DESeq2 documentation for designing the formula. In general, you pass in a column name (e.g. treatment) of your metadata file or a combination of column names (e.g. treatment + cell_type).
#' @param contrast Information regarding the groups between which you would like to perform differential gene expression analysis. It could be between two groups or between multiple groups and needs to follow the following format: contrast = list(A = c(" "), B= c(" ")). If you are comparing two groups (e.g. control vs treatment), the constrast argument should look like the following: contrast = list(A = c("control"), B= c("treatment")). In situations where you have multiple groups to compare- (e.g. control vs treatment1 and treatment2), you should do the following- contrast = list(A = c("control"), B= c("treatment1", "treatment2")).
#' @param variable.genes numeric: The number of most variable genes to be identified. By default, the program identifies the top 1000 most variable genes.
#' @param folder.name Custom folder name that you would like to save your results in.
#' @param save.dir Directory to save the results in. Default: Working Directory.
#' @param custom.gsea User defined gene list to perform GSEA. File need to be supplied as a dataframe with each row as a gene list
#' @param qc Logical. When passed in TRUE, the program would run the quality control modules on the entire dataset. If you are planning to perform multiple comparisons using the contrast argument, run  qc = TRUE for the first time and then change it to qc = FALSE for the subsequent comparisons to speed up the analysis.
#' @param dgea Logical. Parameter to define if differential gene expression analysis is to be performed. Default: TRUE
#' @param ensemblmirror String. Values for the mirror argument are: useast, uswest, asia
#' @import DESeq2
#' @import utils
#' @import plyr
#' @importFrom stats as.formula
#' @return Differentially expressed genes.
#' @examples
#' \dontrun{
#' contrast = list(A = c("control"), B= c("drug_A"))
#' arseq2 (data = example_data,meta = example_meta, design = "treatment", contrast = contrast)
#' }
#' @export

arseq <- function(data,meta,design, contrast, dds=NULL, qc= TRUE, dgea=TRUE, variable.genes=1000, folder.name="ARSeq", custom.gsea=NULL, save.dir=getwd(),ensemblmirror="useast"){

  # Set the working directory
  if (save.dir != getwd()){
    setwd(save.dir)
  }

  # Resolve the contrasts and design
  goi <-  c(paste(unlist(contrast[[1]]), collapse="_"),paste(unlist(contrast[[2]]), collapse="_"))
  groups <- c(unlist(contrast[[1]]),contrast[[2]])

  # Figure out the contrast column programatically
  colum.names <- c()
  statements <- c()
  for (i in 1: length(groups)){
    colum.names <- c(colum.names,colnames(meta[, which(meta == groups[i], arr.ind=T)[1,][[2]], drop=FALSE]))
  }
  # IF the groups of interest are found in multiple colums; find if any one column that contains all groups of interest
  for(i in colum.names){statements <- c(statements,all(groups %in% meta[,i]))}
  for(i in colum.names){if (all(groups %in% meta[,i])){con = i}}

  if (all(colum.names == colum.names[1])){
    contrast.by = colum.names[1]
    # Remname the column with contrast information to arseq.group
    colnames(meta)[which(colnames(meta) == contrast.by)] <- "arseq.group"
  }else if (any(statements)){
    contrast.by = con
    colnames(meta)[which(colnames(meta) == contrast.by)] <- "arseq.group"
    print(paste("Meta Data column that has been chosen for contrast determination is [[",con,"]]"))
  }else{stop('All of your contrast groups for Differntial expression analysis need to be listed under the same column within your metadata file')}

  # resolving design
  design <- gsub(contrast.by, "arseq.group", design)
  design <- paste("~",design,sep = "")
  design <- as.formula(design)

  # Creating a folder to save results
  print("Creating a folder to store your analysis results")
  suppressWarnings(dir.create(folder.name))

  # ENSEMBL ID to gene names
  if (all(substr(row.names(data),1,3) == "ENS")){
    data <- arseq.ensembl2genename(data,ensemblmirror)
    }

  # Remove genes that are completely not expressed
  print("Removing genes that are not expressed: if any")
  data <- data[rowSums(data) > 0, ]

  # Rearrange the samples
  meta <- meta[order(meta$arseq.group),]
  data <- data[,row.names(meta)]

  # Generate the DESEq2 data object
  if(is.null(dds)){
    print("Generating the DESeq object")
    dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = design)
    dds <- DESeq(dds)
  }

  # Get the normalized matrix for heatmaps
  print("Normalizing data and saving the normalized data")
  n_data <- log2(counts(dds, normalized=TRUE)+1)
  write.csv(n_data, file = paste(folder.name,"/normalized_data",".csv",sep = ""))

  # Quality Control of the data
  if (isTRUE(qc)){
    # Create a folder to save results
    suppressWarnings(dir.create(paste(folder.name,"/Quality Control",sep = "")))

    # General structure of the dataset
    print("Performing Quality Control on the data")

    # Calculate Euclidean distance between samples
    euclid.dist <- arseq.euclid.dist (dds)
    euclid.plot <- arseq.euclid.dist.plot (euclid.dist) # Euclidean distance heatmap
    arseq.pheatmap(euclid.plot, filename="Euclidean distance.pdf", save.dir=paste(folder.name,"/Quality Control",sep="")) # Save Euclidean distance.pdf

    # Calculate Poisson Distance between samples
    poisd.dist <- arseq.poisd.dist (dds)
    poisd.plot <- arseq.poisd.dist.plot (poisd.dist) # Poisson distance heatmap
    arseq.pheatmap(poisd.plot, filename="Poisson distance.pdf", save.dir=paste(folder.name,"/Quality Control",sep="")) # Save Poisson distance.pdf

    #PCA plot
    pca.plot <- arseq.pca.plot (dds, intgroup="arseq.group")
    arseq.plot.save (pca.plot, filename="PCA plot.pdf", width=7, height=7, save.dir= paste(folder.name,"/Quality Control",sep="")) # Save the PCA plot

    # Eignen vectors of the PC's
    pca.ev <- arseq.pca (dds)
    write.csv(pca.ev, file = paste(folder.name,"/Quality Control/PCA Eigenvectors.csv",sep="")) # save the eigen vectors

    # MSD plot
    mds.plot <- arseq.mds.plot (dds, intgroup="arseq.group")
    arseq.plot.save (mds.plot, filename="MDS plot.pdf", width=7, height=7, save.dir= paste(folder.name,"/Quality Control",sep="")) # Save the MDS plot

    # Most variable genes
    mvg <- arseq.mvg (dds)
    write.csv(mvg, file = paste(folder.name,"/Quality Control/most variable genes.csv",sep=""))

    # Heatmap of the most variable genes
    mvg.plot <- arseq.mvg.plot (mvg,metadata=meta,intgroup="arseq.group",
                                save.plot=TRUE, save.dir=paste(folder.name,"/Quality Control/",sep=""))
  }
  if (isTRUE(dgea)){
    # Create folder to save the differentially expression results
    suppressWarnings(dir.create(paste(folder.name,"/",goi[1], " vs ", goi[2],sep = "")))
    location <- paste(folder.name,"/",goi[1], " vs ", goi[2],"/",sep = "")

    # Manipulate the dds groups to look at the contrasts of interest
    dds$arseq.group <- mapvalues(dds$arseq.group, from=contrast[[1]], to=rep(paste(unlist(contrast[[1]]), collapse="_"),length(contrast[[1]])))
    dds$arseq.group <- mapvalues(dds$arseq.group, from=contrast[[2]], to=rep(paste(unlist(contrast[[2]]), collapse="_"),length(contrast[[2]])))

    # Subset the normalized matrix and dds based on the groups of interest
    colums_to_subset <- row.names(meta[meta$arseq.group %in% groups,,drop=FALSE]) # Changed
    n_data_goi <- data.frame(n_data[,colums_to_subset])
    dds_subset <- dds[ , dds$arseq.group %in% goi ]
    dds_subset$arseq.group <- factor(dds_subset$arseq.group)

    # Differential Gene Expression Analysis
    print("Performing Differential gene expression analysis between the  constrast groups")
    dds <- DESeq(dds)
    deg <- results(dds, contrast=c("arseq.group",goi[1],goi[2]))
    deg <- deg[order(deg$padj),]

    # Save the DEG matrix
    print("Saving the differentially expressed genes")
    suppressWarnings(dir.create(paste(location,"Differential expression/",sep = "")))
    write.csv(deg, file = paste(location,"Differential expression/",goi[1], " vs ", goi[2],".csv",sep = ""))

    # Heatmap of the differentially expressed genes
    deg.plot <- arseq.deg.plot (deg, data=n_data_goi, dds=dds_subset,
                                save.plot=TRUE, save.dir=paste(location,"Differential expression/",sep=""))

    # Calculate Euclidean distance between contrast groups
    suppressWarnings(dir.create(paste(location,"Sample relation plots",sep = "")))
    euclid.dist.subset <- arseq.euclid.dist (dds_subset)
    euclid.plot.subset <- arseq.euclid.dist.plot (euclid.dist.subset) # Euclidean distance heatmap
    arseq.pheatmap(euclid.plot.subset, filename=paste("Euclidean distance- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                   save.dir=paste(location,"Sample relation plots",sep=""))

    # Calculate Poisson Distance between samples
    poisd.dist.subset <- arseq.poisd.dist (dds_subset)
    poisd.plot.subset <- arseq.poisd.dist.plot (poisd.dist.subset) # Poisson distance heatmap
    arseq.pheatmap(poisd.plot.subset, filename=paste("Poisson distance- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                   save.dir=paste(location,"Sample relation plots",sep=""))

    # PCA
    pca.plot.subset <- arseq.pca.plot (dds_subset, intgroup="arseq.group")
    arseq.plot.save (pca.plot.subset, filename=paste("PCA plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                     width=7, height=7,
                     save.dir=paste(location,"Sample relation plots",sep="")) # Save the PCA plot

    # Eignen vectors of the PC's
    pca.ev.subset <- arseq.pca (dds_subset)
    write.csv(pca.ev.subset, file = paste(location,"Sample relation plots/",goi[1], " vs ", goi[2]," PCA Eigenvectors.csv",sep = "")) # save the eigen vectors

    # MSD plot
    mds.plot.subset <- arseq.mds.plot (dds_subset, intgroup="arseq.group")
    arseq.plot.save (mds.plot.subset, filename=paste("MDS plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                     width=7, height=7,
                     save.dir= paste(location,"Sample relation plots",sep="")) # Save the MDS plot

    # Volcano plot
    volcano.plot <- arseq.volcano.plot (deg)
    arseq.plot.save (volcano.plot, filename=paste("Volcano plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                     width=14, height=14,
                     save.dir= paste(location,"Sample relation plots",sep="")) # Save the MDS plot

    # GO Enrichment Analysis
    suppressWarnings(dir.create(paste(location,"GO enrichment analysis/",sep = "")))
    go.enrich <- arseq.go.enrich (deg,Padj=0.05)
    write.csv(go.enrich, file = paste(location,"GO enrichment analysis/","GO Enrichment- ",goi[1], " vs ", goi[2],".csv", sep = ""))

    # GO Enrichment plot
    go.plot <- arseq.go.enrich.plot (go.enrich)
    arseq.plot.save (go.plot, filename=paste("GO Enrichment plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                     width=6, height=20,
                     save.dir= paste(location,"GO enrichment analysis",sep="")) # Save the GO plot

    # Kegg Pathway Analysis for the DEG's
    suppressWarnings(dir.create(paste(location,"KEGG pathway enrichment/",sep = "")))
    kegg.enrich <- arseq.kegg.enrich (deg, kegg.compare="as.group")
    # write the reults
    write.csv(kegg.enrich$up_regulated_pathways, file = paste(location,"KEGG pathway enrichment/","KEGG Up Regulated Pathways- ",goi[1], " vs ", goi[2],".csv",sep = ""))
    write.csv(kegg.enrich$down_regulated_pathways, file = paste(location,"KEGG pathway enrichment/","KEGG Down Regulated Pathways- ",goi[1], " vs ", goi[2],".csv",sep = ""))

    # Kegg Pathway Analysis Plot
    kegg.plot <- arseq.kegg.enrich.plot (kegg.enrich,foldchanges=kegg.enrich$foldchanges,save.dir= paste(location,"KEGG pathway enrichment",sep=""))
    kegg.plot.2 <- arseq.kegg.enrich.plot (kegg.enrich,foldchanges=kegg.enrich$foldchanges,pathway.plots=FALSE)
    arseq.plot.save (kegg.plot.2, filename=paste("KEGG Enrichment plot- ",goi[1], " vs ", goi[2],".pdf",sep = ""),
                     width=8, height=12,
                     save.dir= paste(location,"KEGG pathway enrichment",sep="")) # Save the KEGG plot

    # GSEA Analysis
    ranked.list <- arseq.gsea.preprocess (deg)
    gsea.output <- arseq.gsea.runall (ranked.list, save=TRUE, save.dir=location, custom.gsea=custom.gsea)
  }
  #if (isTRUE(dgea)){return(deg)}
  return (dds)
  # END
  print(paste("Well done- your analysis is now complete. Head over to [[", getwd(), "]] to view your results"))
}

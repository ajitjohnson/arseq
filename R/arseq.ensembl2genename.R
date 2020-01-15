#' @title Gene ID converter
#' @description Convert ensembl id's to hugo gene names. The function also reduces multiple transcripts into a single gene by taking the sum of all transcripts.
#' @param data Expression matrix with ensemble id's as the first column
#' @return Expression matrix with gene names as the first column
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom stats aggregate
#' @examples
#' \dontrun{
#' data <- arseq.ensembl2genename (example_data)
#' }
#' @export

arseq.ensembl2genename <- function(data){
  print("Converting ENSEMBL ID's to gene names")
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)
  data_m <- merge(data, genes, by.x="row.names", by.y= "ensembl_gene_id")[,-1]
  data_m <- data_m[!(data_m$hgnc_symbol==""), ]
  # Reduce multiple transcripts into a single gene by taking the sum of related transcripts
  data_m$hgnc_symbol <- as.factor(data_m$hgnc_symbol)
  data <- aggregate(.~hgnc_symbol, data=data_m, sum)
  rownames(data) <- data[,1]
  data <- data[,-1]
  return(data)
}


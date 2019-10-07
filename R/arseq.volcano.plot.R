#' @title Volcano Plot
#' @description Volcano plot of the differentially expressed genes.
#' @param deg Differentially expressed genes returned by DESeq2 analysis
#' @param pCutoff Cut-off for statistical significance. A horizontal line will be drawn at -log10(pCutoff). DEFAULT = 0.01. OPTIONAL.
#' @return Volcano plot
#' @import EnhancedVolcano
#' @importFrom stats complete.cases
#' @examples
#' \dontrun{
#' volcano.plot <- arseq.volcano.plot (example_deg)
#' }
#' @export

arseq.volcano.plot <- function(deg,pCutoff = 0.01){
  print("Generating a volcano plot between the constrast groups")
  deg.volcano <- data.frame(deg)[complete.cases(data.frame(deg)),]
  EnhancedVolcano(deg.volcano,
                  lab = rownames(deg.volcano),
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  pCutoff = 0.01,
                  FCcutoff = 2.0,
                  gridlines.minor = FALSE,
                  gridlines.major = FALSE,
                  border = 'full',
                  transcriptLabSize = 4.0,
                  transcriptPointSize = 2.0,
                  colAlpha = 1,
                  legend=c('NS','Log2 FC','Adjusted p-value',
                           'Adjusted p-value (<0.01) & Log2 FC'),
                  legendPosition = 'bottom',
                  legendIconSize = 3.0)
}

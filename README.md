# *arseq*: An Automated RNASeq Analysis Pipeline
This is an easy to use R package for *automated basic RNASeq analysis with minimal coding requirement*. This package is designed to be used by biologists with little to no coding experience. <br><br>
For an in depth tutorial checkout the following [blog post](https://ajitjohnson.com/arseq).

### The package currently supports the following analysis

**Over all data structure analysis:**<br>
  - Euclidean distance between samples
  - Poisson distance between samples
  - PCA analysis
  - PCA Eigenvectors
  - Multidimensional scaling (MDS) analysis
  - Most variable genes<br>

**Analysis between groups of interest:**<br>
  - Differential gene expression analysis using DESeq2
  - Volcano Plot of differentially expressed genes
  - Euclidean distance between samples
  - Poisson distance between samples
  - PCA analysis
  - PCA Eigenvectors
  - Multidimensional scaling (MDS) analysis
  - GO enrichment of the differentially expressed genes
  - KEGG pathway enrichment of the differentially expressed genes
  - KEGG pathway diagrams of the top 5 enriched pathways
  - GSEA analysis (H, C1, C2, C3, C4, C5, C6, C7 genesets)

## Requirements
#### Counts Table
A CSV file with **un-normalized unique genes** as rows and samples as columns. Counts table is generally generated after your FASTQ files have been aligned against the reference genome and quantified (not included in this pipeline). Please note that you will have to provide the un-normalized data as input. Using normalized data, will not work with this package. Instead of gene names, you could also feed in the data with **ENSEMBL ID's**. No other form of ID's is supported at the moment. <br><br>
*Example counts table:*<br>
![Example counts table](https://github.com/ajitjohnson/arseq/blob/master/inst/extdata/data.png)<br><br>

#### Meta data
A CSV file with information regarding the samples. The columns of the count matrix and the rows of the meta data (information about samples) must be in the same order. *arseq* will not make guesses as to which column of the count matrix belongs to which row of the metadata, these must be provided to *arseq* already in a consistent order.<br><br>
*Example meta data file:*<br>
![Example counts table](https://github.com/ajitjohnson/arseq/blob/master/inst/extdata/meta.png)

## How to use
Install and load the package.
```R
# For developmental version
if( !require(devtools) ) install.packages("devtools")
devtools::install_github( "ajitjohnson/arseq" )

# For stable version
install.packages("arseq")

# Load the package
library("arseq")
```
Import your counts matrix and meta data file into R environment.
```R
# Set the working directory (path to the folder of where your data is located)
setwd("\path to the folder \of where your data is located\")

# Load your counts table into R
my_data <- read.csv("counts_table.csv", row.names = 1, header = T) # replace counts_table.csv with your file name

# Load your meta data into R
my_meta <- read.csv("meta_data.csv", row.names = 1, header = T) # replace meta_data.csv with your file name
```
Run the analysis
```R
# Run the analysis. The results will be saved in the same folder as your input data.
arseq (data = my_data, meta = my_meta, design = "treatment", contrast = list(A = c("control"), B= c("treatment1")))
```
In the above command,<br><br>
`design` takes in the column name of the metadata file that contains information regarding the groups you would like to perform differential expression on. You could pass more complex designs- Read the documentation of [DESEq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). As an example, in the above image (metadata file), there a column named treatment that contains information regarding which samples are control samples and which samples were treated with different drugs. So if I want to identify the differentially expressed genes between the control samples and treated samples, I would pass  `design = "treatment"`. <br><br>
`contrast` is another argument that you will need to specify. This is simply the groups of samples between which you would like to perform differential expression analysis. It follows the following format `contrast = list(A = c(" "), B= c(" "))`.<br><br>
If you have three groups in your dataset- *Control, treatment-1 and treatment-2*<br><br>
**Comparison- 1:** To identify the differentially expressed genes between *Control vs treatment-1*, you would pass the contrast in the following manner `contrast = list(A = c("Control"), B= c("treatment-1"))`<br><br>
**Comparison- 2:** To identify the differentially expressed genes between *Control vs treatment-1 + treatment-2*, you would pass the contrast in the following manner `contrast = list(A = c("Control"), B= c("treatment-1", "treatment-2"))`

## Example dataset
The package comes with an example dataset. In order to familiarise yourself with the package and its requirements you could play around with the example dataset.
```R
# view the example counts table
head(example_data)

# view the example meta data
head(example_meta)

# Set the working directory. Folder to which you would like to save your results.
setwd("\path to the folder \that you would like to save the results\")

# Run the analysis. Here we are identifying the differences between control samples and treatment1 samples.
arseq (data = example_data, meta = example_meta, design = "treatment", contrast = list(A = c("control"), B= c("treatment1")))
```

## Additional parameters
The `arseq` function can take in a few additional arguments.<br><br>
`general.stats`- Default is TRUE. This will run the general stat module (e.g. PCA, MDS, etc.. for your entire dataset). If you are making multiple comparisons using the `contrast` argument, run  `general.stats = TRUE` for the first time and change it to `general.stats = FALSE` for the subsequent comparisons to speed up the analysis.<br><br>
`variable.genes`- Number of variable genes to be identified. By default the program identifies the top 1000 most variable genes. you could set it to `variable.genes=3000` to calculate the top 3000 most variable genes.

## Cite
If you found this package useful, please do cite this page in your publication. Thank you.

## Issues and Features
If there are any issues please report it at https://github.com/ajitjohnson/arseq/issues

## Additional information
For an in depth tutorial checkout the following [blog post](https://ajitjohnson.com/arseq).<br>
You can also [tweet](https://twitter.com/ajitjohnson_n) me directly for inclusion of new methods into this package.

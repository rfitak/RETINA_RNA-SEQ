# Differential Expression using DESEQ2.
In this section we will look at differential expression using DESEQ2 ([Love at al. 2014, doi: 10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)).  DESEQ2 is able to process each individual gene as a linear model, where we can look at the effect of individual variables in addition to the interaction between them.  We will begin from the count data generated during the mapping with STAR (file [counts.tsv](./Data/counts.tsv) from the "Data folder").  we will also use the samples file [samplesheet.tsv](./Data/samplesheet.tsv) from the data folder.
```
# Make a new working folder
mkdir DESEQ2
cd DESEQ2

# Copy necessary files over
cp ../Data/samplesheet.tsv .
cp ../Data/counts.tsv .
```

Since DESEQ2 is an R package, all the following code is to be used in an R terminal.
```
# Load required packages
library(DESeq2)

# Read in sample information spreadsheet
samplesheet = read.table("samplesheet.tsv", sep = "\t", header = T)

# Make a list of sample names to match those
samplesheet$Library = paste0("Sample_", samplesheet$Library)

```

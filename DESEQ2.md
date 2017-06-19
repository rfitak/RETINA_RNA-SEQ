# Differential Expression using DESEQ2.
In this section we will look at differential expression using DESEQ2 ([Love at al. 2014, doi: 10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)).  DESEQ2 is able to process each individual gene as a linear model, where we can look at the effect of individual variables in addition to the interaction between them.  We will begin from the count data generated during the mapping with STAR (file [counts.tsv](./Data/counts.tsv) from the "Data folder").  we will also use the samples file [samplesheet.tsv](./Data/samplesheet.tsv) from the data folder.
```bash
# Make a new working folder
mkdir DESEQ2
cd DESEQ2

# Copy necessary files over
cp ../Data/samplesheet.tsv .
cp ../Data/counts.tsv .
```

Since DESEQ2 is an R package, all the following code is to be used in an R terminal.
```R
# Load required packages
library(DESeq2)

# Read in sample information spreadsheet
samplesheet = read.table("samplesheet.tsv", sep = "\t", header = T)

# Make a list of sample names to match those in the counts file
samplesheet$Library = paste0("Sample_", samplesheet$Library)

# Load count data as a matrix
data = as.matrix(read.table("counts.tsv", 
   sep = "\t", row.names = 1, header=T))

# Build DESEQ data object
dds = DESeqDataSetFromMatrix(countData = data,
   colData = DataFrame(samplesheet),
   design = ~ Eye + Treatment + Eye:Treatment)

# Relevel treatment to expression is relative to control
dds$Treatment = relevel(dds$Treatment, "CONTROL")

# Remove lowly expressed genes
dds.trim = dds[rowSums(counts(dds))>0,]

# Run DESEQ2
dds.trim = DESeq(dds.trim)

# Filter Results for false discovery rate (FDR) < 0.05
res = results(dds.trim, alpha = 0.05)

# Re-order results by FDR
res.ordered = res[order(res$padj),]
```

## Visualizations

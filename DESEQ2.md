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
dds.trim = dds[rowSums(counts(dds))>10,]

# Run DESEQ2
dds.trim = DESeq(dds.trim, betaPrior = F)

# Get a list of coefficients
resultsNames(dds.trim)
   # [1] "Intercept"                   "Eye_R_vs_L"                 
   # [3] "Treatment_PULSED_vs_CONTROL" "EyeR.TreatmentPULSED" 

# Filter Results for false discovery rate (FDR) < 0.05
res = results(dds.trim, alpha = 0.05)

# Re-order results by FDR
res.ordered = res[order(res$padj),]
   # "EyeR.TreatmentPULSED" : 0 DE genes, 0 FDR < 0.1
   
# Start a list of results.
DESeq2.results=list()

# Append to the list of results
DESeq2.results=c(DESeq2.results, EyeR.TreatmentPULSED = res.ordered)

# Summary of R vs L eye, added to the list of results.
res = results(dds.trim, name = "Eye_R_vs_L", alpha = 0.05)
res.ordered = res[order(res$padj),]
DESeq2.results=c(DESeq2.results, Eye_R_vs_L = res.ordered)
   # "Eye_R_vs_L" : 0 DE genes, 0 FDR < 0.1

# Summary of Pulsed vs Control, added to the list of results.
res = results(dds.trim, name = "Treatment_PULSED_vs_CONTROL", alpha = 0.05)
res.ordered = res[order(res$padj),]
DESeq2.results=c(DESeq2.results, Treatment_PULSED_vs_CONTROL = res.ordered)
   # "Treatment_PULSED_vs_CONTROL" : 1 DE genes, 5 FDR < 0.1
```
The above code looks for differential expression either by comparing overall:
- CONTROL vs PULSED
- LEFT vs RIGHT
- Additive effect of the pulse within a group (EyeR.TreatmentPULSED).


Next, we will repeat the above analysis but this time using the "Group" variable only.  We now can compare different combinations of groups, but lack an interaction effect.
```R
# Starting from the 'data' object above
# Build DESEQ data object
dds = DESeqDataSetFromMatrix(countData = data,
   colData = DataFrame(samplesheet),
   design = ~ Group)
   
# Remove lowly expressed genes
dds.trim = dds[rowSums(counts(dds))>10,]

# Run DESEQ2
dds.trim = DESeq(dds.trim, betaPrior = F)

# Get a list of coefficients
resultsNames(dds.trim)
   # [1] "Intercept"                    "Group_L_PULSED_vs_L_CONTROL" 
   # [3] "Group_R_CONTROL_vs_L_CONTROL" "Group_R_PULSED_vs_L_CONTROL"

# Get results for group of interest with false discovery rate (FDR) < 0.05, then append to results
res = results(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), alpha = 0.05)
res.ordered = res[order(res$padj),]

# Apply shrinkage for very low-expressed genes
resLFC <- lfcShrink(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), res = res)






DESeq2.results=c(DESeq2.results, L_PULSED_vs_L_CONTROL = res.ordered)
   # 0 DE genes, 2 genes FDR < 0.1 
res = results(dds.trim, contrast = c("Group", "R_PULSED", "R_CONTROL"), alpha = 0.05)
res.ordered = res[order(res$padj),]
DESeq2.results=c(DESeq2.results, R_PULSED_vs_R_CONTROL = res.ordered)
   # 1 DE gene, 1 gene FDR < 0.1
res = results(dds.trim, contrast = c("Group", "R_CONTROL", "L_CONTROL"), alpha = 0.05)
res.ordered = res[order(res$padj),]
DESeq2.results=c(DESeq2.results, R_CONTROL_vs_L_CONTROL = res.ordered)
   # 0 DE genes, 0 genes FDR < 0.1 
res = results(dds.trim, contrast = c("Group", "R_PULSED", "L_PULSED"), alpha = 0.05)
res.ordered = res[order(res$padj),]
DESeq2.results=c(DESeq2.results, R_PULSED_vs_L_PULSED = res.ordered)
   # 0 DE genes, 0 genes FDR < 0.1 

# Save all results as an R data file
save(DESeq2.results, file = "DESeq2-results.R")
```
All the differential expression results are now saved as an R data file that can be loaded anytime using `load("DESeq2-results.R")`.
 
## Visualizations
```R
# Convert to log transformed data
rld <- rlog(dds.trim)

# Plot PCA/MDS
plotPCA(rld, intgroup="Group", ntop=1000)

# Plot Dendrogram
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
hc <- hclust(distsRL)
hc$labels <- paste(samplesheet$Group, "_", c(1:24))
plot(hc)

# MA plot
```

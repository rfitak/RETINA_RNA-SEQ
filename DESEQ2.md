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
res = results(dds.trim, name = "EyeR.TreatmentPULSED", alpha = 0.05)

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


Next, we will repeat the above analysis but this time using the "Group" variable only.  We now can compare different combinations of groups, but lack an interaction effect.  By default, a beta Prior is included which shrinks the LFC estimates.
```R
# Starting from the 'data' object above
# Build DESEQ data object
dds = DESeqDataSetFromMatrix(countData = data,
   colData = DataFrame(samplesheet),
   design = ~ Group)
   
# Remove lowly expressed genes
dds.trim = dds[rowSums(counts(dds))>10,]

# Run DESEQ2 (with beta prior, which includes shrinking LFC)
dds.trim = DESeq(dds.trim, betaPrior = T)

# Get a list of coefficients
resultsNames(dds.trim)
   # [1] "Intercept"      "GroupL_CONTROL" "GroupL_PULSED"  "GroupR_CONTROL"
   # [5] "GroupR_PULSED"

# Get results for group of interest with false discovery rate (FDR) < 0.05, then append to results
res1 = results(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), alpha = 0.05)
res2 = results(dds.trim, contrast = c("Group", "R_PULSED", "R_CONTROL"), alpha = 0.05)
res3 = results(dds.trim, contrast = c("Group", "R_CONTROL", "L_CONTROL"), alpha = 0.05)
res4 = results(dds.trim, contrast = c("Group", "R_PULSED", "L_PULSED"), alpha = 0.05)
res1.ordered = res1[order(res1$padj),]   # 0 DE genes, 2 genes FDR < 0.1
res2.ordered = res2[order(res2$padj),]   # 1 DE gene, 1 gene FDR < 0.1
res3.ordered = res3[order(res3$padj),]   # 0 DE genes, 0 genes FDR < 0.1
res4.ordered = res4[order(res4$padj),]   # 0 DE genes, 0 genes FDR < 0.1

# Merge to results (bp = betaprior)
DESeq2.results=c(DESeq2.results, LPvsLC_bpT = res1.ordered, RPvsRC_bpT = res2.ordered, RCvsLC_bpT = res3.ordered, RPvsLP_bpT = res4.ordered)
```
Alternatively, you can run the analysis without the beta prior and perform the shrinkage of LFC estimates separately.  However, the p-values are calculated from the un-shrunken estimates!
```R
# Starting from the 'data' object above
# Build DESEQ data object
dds = DESeqDataSetFromMatrix(countData = data,
   colData = DataFrame(samplesheet),
   design = ~ Group)

# Remove lowly expressed genes
dds.trim = dds[rowSums(counts(dds))>10,]

# Run DESEQ2 (without beta prior, shrinking LFC after)
dds.trim = DESeq(dds.trim, betaPrior = F)

# Get a list of coefficients
resultsNames(dds.trim)
   # [1] "Intercept"                    "Group_L_PULSED_vs_L_CONTROL" 
   # [3] "Group_R_CONTROL_vs_L_CONTROL" "Group_R_PULSED_vs_L_CONTROL"

# Get results for group of interest with false discovery rate (FDR) < 0.05, then append to results
res1 = results(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), alpha = 0.05)
res1.LFC <- lfcShrink(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), res = res1)
res2 = results(dds.trim, contrast = c("Group", "R_PULSED", "R_CONTROL"), alpha = 0.05)
res2.LFC <- lfcShrink(dds.trim, contrast = c("Group", "R_PULSED", "R_CONTROL"), res = res2)
res3 = results(dds.trim, contrast = c("Group", "R_CONTROL", "L_CONTROL"), alpha = 0.05)
res3.LFC <- lfcShrink(dds.trim, contrast = c("Group", "R_CONTROL", "L_CONTROL"), res = res3)
res4 = results(dds.trim, contrast = c("Group", "R_PULSED", "L_PULSED"), alpha = 0.05)
res4.LFC <- lfcShrink(dds.trim, contrast = c("Group", "R_PULSED", "L_PULSED"), res = res4)
res1.ordered = res1.LFC[order(res1.LFC$padj),]
res2.ordered = res2.LFC[order(res2.LFC$padj),]
res3.ordered = res3.LFC[order(res3.LFC$padj),]
res4.ordered = res4.LFC[order(res4.LFC$padj),]

# Merge to results (bp = betaprior)
DESeq2.results=c(DESeq2.results, LPvsLC_bpF = res1.ordered, RPvsRC_bpF = res2.ordered, RCvsLC_bpF = res3.ordered, RPvsLP_bpF = res4.ordered)

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
res1 = results(dds.trim, contrast = c("Group", "L_PULSED", "L_CONTROL"), alpha = 0.05, addMLE = T)
plotMA(res1, MLE = TRUE, alpha = 0.05, ylim = c(-3, 3)) # unshrunk estimates
plotMA(res1, MLE = FALSE, alpha = 0.05, ylim = c(-3, 3)) # shrunken estimates
```

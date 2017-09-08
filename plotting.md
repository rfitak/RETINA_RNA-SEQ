# Additional plotting functions
The various plotting code here is to compare the expression results in the retinae from both the Cufflinks and DESeq2 analysis with the brain expression results from our previous study published in Biology Letters.

```R
# Get list of significant brain genes
a=read.table("/media/rfitak/Seagate\ Expansion\ Drive/TROUT_PROJECT/CUFFLINKS/CUFFDIFF/5-DAYS/gene_exp.diff", sep="\t", header=T)
genes = as.character(subset(a, significant == "yes")$test_id)

# Read in Retina data
load("~/Desktop/TROUT_RETINA_PROJECT/DESEQ2/DESeq2-results.R")

# Get the left and right retina data for these genes
L = DESeq2.results$LPvsLC_bpT[order(rownames(DESeq2.results$LPvsLC_bpT)),]
R = DESeq2.results$RPvsRC_bpT[order(rownames(DESeq2.results$RPvsRC_bpT)),]
retina.R = DESeq2.results$RPvsRC_bpT[which(rownames(DESeq2.results$RPvsRC_bpT) %in% genes),]
retina.L = DESeq2.results$LPvsLC_bpT[which(rownames(DESeq2.results$LPvsLC_bpT) %in% genes),]

# Reorder the retina data by gene name
retina.R = retina.R[order(rownames(retina.R)),]
retina.L = retina.L[order(rownames(retina.L)),]

# Check that gene names match
for (i in 1:nrow(retina.L)){
   if (rownames(retina.L)[i] != rownames(retina.R)[i]) print("Failed")
}

# Because not all the DE brain genes were present in the retina data,
# grab these genes from the brain dataset and reorder
brain.sig = a[which(as.character(a$test_id) %in% rownames(retina.R)),]

# Check that gene names match
for (i in 1:nrow(retina.R)){
   if (rownames(retina.R)[i] != brain.sig$test_id[i]) print("Failed")
}

# Merge datasets
data = cbind(retina.R$log2FoldChange, retina.L$log2FoldChange, brain.sig$log2.fold_change.)
```
The log2FC values for the three comparisons are now stored in the variable `data`. These three comparison are:
- Retina: Right-pulsed vs Right-control
- Retina: Left-pulsed vs Left-control
- Brain: Pulsed vs Control

```R
# Plotting
colors = rep("black", nrow(R))
colors[which(rownames(L) %in% genes)] = "red"
par(mfrow=c(2,2))

# Plot 1
plot.new()
plot.window(xlim = c(-3,3), ylim = c(-3,3))
box()
axis(1)
axis(2, las = 1)
title(xlab = "Right retina", ylab = "Left retina")
#rect(-1,-1,1,1, col = rgb(211/255, 211/255, 211/255, alpha = 0.5), border = NA)
abline(0,1, lty = "dashed")
points(R$log2FoldChange, L$log2FoldChange, pch = 1)
points(data[,1], data[,2], pch = 1, col = "red")

# Plot 2
plot.new()
plot.window(xlim = c(-3,3), ylim = c(-3,3))
box()
axis(1)
axis(2, las = 1)
text(0, 0, "CUFFLINKS DATA HERE")
#title(xlab = "Right retina", ylab = "Left retina")
#rect(-1,-1,1,1, col = rgb(211/255, 211/255, 211/255, alpha = 0.5), border = NA)
#abline(0,1, lty = "dashed")
#points(data[,1], data[,2], pch = 1)

# Plot 3
plot.new()
plot.window(xlim = c(-3,3), ylim = c(-3,3))
box()
axis(1)
axis(2, las = 1)
title(xlab = "Right retina", ylab = "Brain")
#rect(-1,-1,1,1, col = rgb(211/255, 211/255, 211/255, alpha = 0.5), border = NA)
abline(0,1, lty = "dashed")
points(data[,1], data[,3], pch = 1)

# Plot 4
plot.new()
plot.window(xlim = c(-3,3), ylim = c(-3,3))
box()
axis(1)
axis(2, las = 1)
title(xlab = "Left retina", ylab = "Brain")
#rect(-1,-1,1,1, col = rgb(211/255, 211/255, 211/255, alpha = 0.5), border = NA)
abline(0,1, lty = "dashed")
points(data[,2], data[,3], pch = 1)
```

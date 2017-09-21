# Additional plotting functions
The various plotting code here is to compare the expression results in the retinae from both the Cufflinks and DESeq2 analysis with the brain expression results from our previous study published in Biology Letters.

First, make a plot to compare LC in the left and right retinas from CUFFLINKS
```R
# Load libraries
library(cummeRbund)
library(DESeq2)

# Load CUFFLINKS results
cuff = readCufflinks("DIFF")

# Make table of Left and Right LFC
L = diffData(genes(cuff), "LP", "LC")[,c(1,5,8)]
R = diffData(genes(cuff), "RP", "RC")[,c(1,5,8)]
LR = cbind(L, R)
colnames(LR) = c("ID_L","Status_L", "LFC_L","ID_R", "Status_R", "LFC_R" )
LR = subset(LR, Status_L == "OK" & Status_R == "OK")
   # Result: 65110 loci shared as 'OK'
LR = subset(LR, LFC_L < Inf & LFC_L > -Inf & LFC_R < Inf & LFC_R > -Inf)
   # Result: 64805 loci shared as not infinite values

# Build linear model
fit1 <- lm(LFC_L ~ LFC_R, data = LR)
summary(fit1)
   # Result:
      # Residuals:
      #      Min       1Q   Median       3Q      Max 
      # -10.9580  -0.2331  -0.0274   0.2022   6.2554 
      # Coefficients:
      #             Estimate Std. Error t value Pr(>|t|)    
      # (Intercept) 0.020621   0.002302    8.96   <2e-16 ***
      # LFC_R       0.382155   0.003279  116.56   <2e-16 ***
      # ---
      # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      # Residual standard error: 0.5847 on 64803 degrees of freedom
      # Multiple R-squared:  0.1733,	Adjusted R-squared:  0.1733 
      # F-statistic: 1.359e+04 on 1 and 64803 DF,  p-value: < 2.2e-16

# Corellation test
cor.test(LR$LFC_L,LR$LFC_R)
   # Result:
      # t = 116.56, df = 64803, p-value < 2.2e-16
      # alternative hypothesis: true correlation is not equal to 0
      # 95 percent confidence interval:
      #  0.4099161 0.4226460
      # sample estimates:
      #       cor 
      # 0.4163014

# Plot the L vs R LFC comparison
p1 = ggplot(LR, aes(x = LFC_R, y = LFC_L))
p1 = p1 + geom_point(alpha = 0.25, na.rm = T)
p1 = p1 + xlim(c(-10,10))
p1 = p1 + ylim(c(-10,10))
p1 = p1 + geom_smooth(method='lm',formula=y~x, fullrange = T, se = F, na.rm = T)
p1 = p1 + xlab(expression('Right Retina Log'[2]*FC))
p1 = p1 + ylab(expression('Left Retina Log'[2]*FC))
p1 = p1 + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
p1
```


### OLD STUFF BELOW

















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

# Assembly and Quantification using CUFFLINKS v2.2.1
All steps are guided by, or include, the reference annotations we made earlier. The general steps of this set of code are to:
1. Assemble reads from each sample separately
2. Merge the assemblies together
3. Quantify the expression
4. Calculate differential expression
5. Visualize in R

## Step 1:  Assembly
```bash
# Make Cufflinks folder
mkdir CUFFLINKS
cd CUFFLINKS

# Setup a loop for the 24 files
for i in {1..24}; do
   # Make a folder for the individual
   mkdir ASSEMBLY_${i}

   # Run the CUFFLINKS assembly step
   cufflinks \
      -p 8 \
      -o ASSEMBLY_${i} \
      --GTF-guide /work/frr6/RETINA/TROUT_REF/O_mykiss.annot.gtf \
      --frag-bias-correct /work/frr6/RETINA/TROUT_REF/Omykiss.genome.fa \
      --multi-read-correct \
      --library-type fr-firststrand \
      --verbose \
      /work/frr6/RETINA/MAPPING/LIBRARY_${i}/*.mapping-Aligned.sortedByCoord.out.bam

   # Move up a folder
   cd ..
   
# End the loop
done
```
Description of the above parameters:
- -p 8 :: use 8 CPUs (threads)
- -o \<folder\> :: folder to place the output files
- --GTF-guide \<file\>:: GTF file of annotations
- --frag-bias-correct \<file\> :: use bias correction - reference fasta required
- --multi-read-correct :: use 'rescue method' for multi-reads (more accurate)
- --library-type fr-firststrand :: library type, normal for dUTP protocols
- --verbose :: log-friendly verbose processing

## Step 2: Merge Assemblies
```bash
# Make MERGED assembly folder
mkdir MERGED
cd MERGED

# Make list of all merged assembly files
ls /work/frr6/RETINA/CUFFLINKS/ASSEMBLY_*/transcripts.gtf > assembly_GTF_list.txt

# Run cuffmerge
cuffmerge \
   -p 8 \
   -o ./ \
   --ref-gtf /work/frr6/RETINA/TROUT_REF/O_mykiss.annot.gtf \
   --ref-sequence /work/frr6/RETINA/TROUT_REF/Omykiss.genome.fa \
   --min-isoform-fraction 0.05 \
   assembly_GTF_list.txt
```
Description of parameters:
- -p 8 :: use 8 CPUs (threads)
- -o ./ :: folder to place the output files (current folder)
- --ref-gtf \<file\> :: An optional "reference" annotation GTF of original annotations
- --ref-sequence \<file\> :: Genomic DNA sequences for the reference
- --min-isoform-fraction 0.05 :: Discard isoforms with abundance below this.  Range is 0-1

## Step 3: Quantify Expression
```bash
cd /work/frr6/RETINA/CUFFLINKS

# Make new quantification folder
mkdir QUANT
cd QUANT

# Loop through all 24 files to quanitfy
for i in {1..24}; do

   # Make new folder for each sample
   mkdir LIBRARY_${i}
 
   # Run cuffquant
   cuffquant \
      -p 8 \
      -o LIBRARY_${i} \
      --multi-read-correct \
      --library-type fr-firststrand \
      --verbose \
      --frag-bias-correct /work/frr6/RETINA/TROUT_REF/Omykiss.genome.fa \
      /work/frr6/RETINA/CUFFLINKS/MERGED/merged.gtf \
      /work/frr6/RETINA/MAPPING/LIBRARY_${i}/*.mapping-Aligned.sortedByCoord.out.bam

   # Move up a folder
   cd ..
   
# End for loop
done
```
Description of parameters:
- -p 8 :: use 8 CPUs (threads)
- -o ./ :: folder to place the output files (current folder)
- --multi-read-correct :: use 'rescue method' for multi-reads (more accurate)
- --library-type fr-firststrand :: library type, normal for dUTP protocols
- --verbose :: log-friendly verbose processing
- --frag-bias-correct \<file\> :: use bias correction - reference fasta required
- location of merged annotation file
- location of bam file to quanitfy expression

Now clean things up a bit
```bash
# Rename abundance files
for i in {1..24}; do
mv LIBRARY_${i}/abundances.cxb ${i}_abundances.cxb
rm -rf LIBRARY_${i}
done
```
## Step 4:  Differential Expression (DE) Analysis
```bash
cd /work/frr6/RETINA/CUFFLINKS

# Make a folder for DE analysis
mkdir DIFF
cd DIFF

# Note: sample groups
   # LC = left control group
   # RC = right control group
   # LP = left pulsed group
   # RP = right pulsed group

# Setup path to folder of abundance files
pth=/work/frr6/RETINA/CUFFLINKS/QUANT

# Run cuffdiff - all pairwise group comparisons will be made
cuffdiff \
   -p 8 \
   -o ./ \
   -L LC,RC,LP,RP \
   --frag-bias-correct /work/frr6/RETINA/TROUT_REF/Omykiss.genome.fa \
   --multi-read-correct \
   --min-alignment-count 10 \
   --FDR 0.05 \
   --library-type fr-firststrand \
   --verbose \
   /work/frr6/RETINA/CUFFLINKS/MERGED/merged.gtf \
   ${pth}/3_abundances.cxb,${pth}/7_abundances.cxb,${pth}/9_abundances.cxb,${pth}/17_abundances.cxb,${pth}/19_abundances.cxb,${pth}/21_abundances.cxb \
   ${pth}/4_abundances.cxb,${pth}/8_abundances.cxb,${pth}/10_abundances.cxb,${pth}/18_abundances.cxb,${pth}/20_abundances.cxb,${pth}/22_abundances.cxb \
   ${pth}/1_abundances.cxb,${pth}/5_abundances.cxb,${pth}/11_abundances.cxb,${pth}/13_abundances.cxb,${pth}/15_abundances.cxb,${pth}/23_abundances.cxb \
   ${pth}/2_abundances.cxb,${pth}/6_abundances.cxb,${pth}/12_abundances.cxb,${pth}/14_abundances.cxb,${pth}/16_abundances.cxb,${pth}/24_abundances.cxb
```
Description of parameters:
- -p 8 :: use 8 CPUs (threads)
- -o \<folder\> :: folder to place the output files
- -L \<first group,second group\> :: comma-separated list of condition labels 
- --frag-bias-correct \<file\> :: use bias correction - reference fasta required
- --multi-read-correct :: use 'rescue method' for multi-reads (more accurate)
- --min-alignment-count 10 :: minimum number of alignments in a locus for testing, 10 by default
- --FDR 0.05 :: False discovery rate cutoff
- --library-type fr-firststrand :: library type, normal for dUTP protocols
- --verbose :: log-friendly verbose processing
- location of merged annotation file
- comma separated list of files for first group (see --L) (LC)
- comma separated list of files for second group (see --L) (RC)
- comma separated list of files for third group (see --L) (LP)
- comma separated list of files for fourth group (see --L) (RP)

Alternatively, the analysis can also be run using a "contrasts" file to limit the comparisons made to a subset of all pairwise comparisons.  Simply add the parameter:
--contrast-file (file containing list of contrasts (comparisons) to make.  Simple two column design.)

## Step 5: Visualization in R using CUMMERBUND v2.18
```
# R code
library(cummeRbund)
cuff = readCufflinks(".")

# Get list of significant genes
sigGeneIds<-getSig(cuff,alpha=0.05,level="genes")

# Plot Dendrogram
hc = hclust(JSdist(makeprobs(log10(repFpkmMatrix(genes(cuff))+1))))
plot(hc)
```
Results:
CuffSet instance with: 
- 4 samples  
- 126060 genes  
- 273985 isoforms  
- 162818 TSS  
- 46597 CDS  
- 756360 promoters  
- 976908 splicing  
- 252690 relCDS  
No Significant genes, isoforms, TSS, or CDS detected  






# Assembly and Quantification using CUFFLINKS v2.2.1
All steps are guided by, or include, the reference annotations we made earlier. The general steps of this set of code is to:
1. Assemble reads from each sample separately
2. Merge the assemblies together
3. Quantify the expression
4. Calculate differential expression

# Step 1:  Assembly
```
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

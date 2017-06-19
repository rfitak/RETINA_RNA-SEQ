# Read Counting Tools
This set of scripts is to count reads mapping to each gene (expression).
First, the counting of HTSeq v0.7.2 was compared with that of STAR.
I only verified the methods using one BAM file.
```
# Move into the Mapping folder
cd /work/frr6/RETINA/MAPPING

# Count using HTSeq
htseq-count \
   -f bam \
   -r pos \
   -s reverse \
   -a 10 \
   -t cds \
   -m union \
   -i gene_id \
   LIBRARY_1/DC0053L.mapping-Aligned.sortedByCoord.out.bam \
   ../TROUT_REF/O_mykiss.annot.gtf > test.counts
```
I compared the above results with the fourth column of STAR's output file <basename>-ReadsPerGene.out.tab.  The results were identical so no need to count using HTSeq.

Make final table of counts for each file from the STAR output (only the fourth column)
```
# Build table of Raw Counts
cat \
   <(echo "Gene") \
   <(cut -f1 LIBRARY_1/DC0053L.mapping-ReadsPerGene.out.tab | tail -n +5) \
   > counts.tsv
for i in {1..24}
   do
   cut -f4 LIBRARY_${i}/*.mapping-ReadsPerGene.out.tab | \
      tail -n +5 | \
      cat <(echo "Sample_${i}") - > tmp
   paste counts.tsv tmp > tmp2
   mv tmp2 counts.tsv
   rm tmp
   done
```
The file [counts.tsv](./Data/counts.tsv) can be found in the "Data" folder for this Github repository.  It will be used in downstream analyses.

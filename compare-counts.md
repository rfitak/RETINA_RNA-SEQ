# Compare Counting Tools
This set of scripts is to compare the output of STAR's gene counting (expression) with that of HTSeq v0.7.2.
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


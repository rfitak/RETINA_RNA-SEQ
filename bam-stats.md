# Calculate summary statistics from the BAM files
Collect several summary statistics from the BAM files using SAMTOOLS v1.3
```
# Get statistics for bam files
for i in {1..24}
   do
   cd LIBRARY_${i}
   samtools1.3 stats *.mapping-Aligned.sortedByCoord.out.bam > $i.genome.bamstats
   samtools1.3 stats *.mapping-Aligned.toTranscriptome.out.bam > $i.transcriptome.bamstats
   cd ..
done

# Make a header for a summary output file
echo "Library \
   Genome_reads \
   Genome_mapped_reads \
   Transcriptome_reads \
   Transcriptome_mapped_reads" | \
   tr -s " " | \
   sed "s/ /\t/g" > bamstats.tsv
   
# Grab various statistics from each file and append ot the summary file
for i in {1..24}
   do
   a=$(grep $'^SN\tsequences:' LIBRARY_${i}/${i}.genome.bamstats | cut -f3)
   b=$(grep $'^SN\treads mapped:' LIBRARY_${i}/${i}.genome.bamstats | cut -f3)
   c=$(grep $'^SN\tsequences:' LIBRARY_${i}/${i}.transcriptome.bamstats | cut -f3)
   d=$(grep $'^SN\treads mapped:' LIBRARY_${i}/${i}.transcriptome.bamstats | cut -f3)
   echo "Library_${i} $a $b $c $d" | sed "s/ /\t/g" >> bamstats.tsv
done
```

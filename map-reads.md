# Mapping the trimmed and paired reads to the indexed reference genome using STAR v020201
This section performs the mapping, or alignment, of reads to the reference genome.  The following code generates a new folder numbered according to their order in the seqfiles.list file.  The numbers are provided by the SLURM array job queue as $SLURM_ARRAY_TASK_ID
```
# Make a mapping folder and move inside it.
mkdir MAPPING
cd MAPPING

# Make a new folder
mkdir LIBRARY_${SLURM_ARRAY_TASK_ID}

# Get the base name of the forward and reverse reads
file=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../RAW_SEQS/seqfiles.list)

STAR \
   --runMode alignReads \
   --runThreadN 8 \
   --genomeDir /work/frr6/RETINA/TROUT_REF \
   --readFilesIn /work/frr6/RETINA/TRIMMED_SEQS/${file}_trimmed.F.fastq.gz /work/frr6/RETINA/TRIMMED_SEQS/${file}_trimmed.R.fastq.gz \
   --readFilesCommand zcat \
   --outFileNamePrefix LIBRARY_${SLURM_ARRAY_TASK_ID}/${file}.mapping- \
   --quantMode TranscriptomeSAM GeneCounts \
   --outReadsUnmapped Fastx \
   --outSAMtype BAM SortedByCoordinate \
   --outFilterIntronMotifs RemoveNoncanonical \
   --twopassMode Basic \
   --outSAMmode Full \
   --outSAMstrandField none \
   --outSAMattributes All \
   --outSAMunmapped Within \
   --quantTranscriptomeBan IndelSoftclipSingleend
```
Here is a description of the parameters:
- runMode alignReads :: tell STAR to align reads
- runThreadN 8 :: use 8 CPUs (threads)
- genomeDir /work/frr6/RETINA/TROUT_REF :: location of the folder with the indexed reference genome
- readFilesIn /work/frr6/RETINA/TRIMMED_SEQS/${file}_trimmed.F.fastq.gz /work/frr6/RETINA/TRIMMED_SEQS/${file}_trimmed.R.fastq.gz :: location of forward and reverse reads files
- readFilesCommand zcat :: this is used if the read files are compressed
- outFileNamePrefix LIBRARY_${SLURM_ARRAY_TASK_ID}/${file}.mapping- :: prefix for output file names
- quantMode TranscriptomeSAM GeneCounts :: Produce counts in a SAM file as well as that produced by HTseq
- outReadsUnmapped Fastx :: Produce a separate file of unmapped reads
- outSAMtype BAM SortedByCoordinate :: Produce a sorted BAM file output
- outFilterIntronMotifs RemoveNoncanonical :: filter out alignments that contain non-canonical junctions
- twopassMode Basic :: basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
- outSAMmode Full :: full SAM output format
- outSAMstrandField none :: NA
- outSAMattributes All :: provide all attributes "NH HI AS nM NM MD jM jI"
- outSAMunmapped Within :: output unmapped reads within the main SAM file
- quantTranscriptomeBan IndelSoftclipSingleend :: prohibit indels, soft clipping and single-end alignments - compatible with RSEM

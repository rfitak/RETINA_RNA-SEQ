# This code is for generating the index using STAR v020201
This set of code is for generating an index of the trout reference genome using STAR.  In the previous section, the reference genome, annotations, and mitochondrial genome with annotations were all merged together.
```
# Move into the correct folder
cd /work/frr6/RETINA/TROUT_REF

# Index using STAR
STAR \
   --runThreadN 8 \
   --runMode genomeGenerate \
   --genomeDir . \
   --genomeFastaFiles Omykiss.genome.fa \
   --sjdbGTFfeatureExon CDS \
   --sjdbGTFtagExonParentTranscript transcript_id \
   --sjdbGTFtagExonParentGene gene_id \
   --sjdbGTFfile O_mykiss.annot.gtf \
   --sjdbOverhang 74
```
Here is a description of the parameters:
- runThreadN 8 :: use 8 CPUs (threads)
- runMode genomeGenerate :: tells STAR to generate an index
- genomeDir . :: write output files to the current folder
- genomeFastaFiles Omykiss.genome.fa :: reference genome file
- sjdbGTFfeatureExon CDS :: feature type in GTF file to be used as exons for building transcripts
- sjdbGTFtagExonParentTranscript transcript_id :: tag name to be used as exons' transcript-parents (default "transcript_id" works for GTF files)
- sjdbGTFtagExonParentGene gene_id :: tag name to be used as exons' gene-parents (default "gene_id" works for GTF files)
- sjdbGTFfile O_mykiss.annot.gtf :: path to the GTF file with annotations
- sjdbOverhang 74 :: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)

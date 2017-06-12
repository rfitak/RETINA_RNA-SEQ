# RETINA_RNA-SEQ
The code for analysis of the trout retina RNA-seq data.  The outline and associated code are found below:
1. [Trim the raw sequencing data](./trim-seqs.md)
2. [Download and prepare the reference genome and annotations](./prep-reference.md)
3. [Index the reference genome](./index-reference.md)
4. [Map trimmed reads to reference](./map-reads.md)
 1. [Calculate summary mapping statistics](./bam-stats.md)
 2. [Verify quantification with STAR and HTSeq](./compare-counts.md)
5. [Assembly and quantification in CUFFLINKS](./cufflinks-workflow.md)

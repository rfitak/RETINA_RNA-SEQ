# Differential expression in the retina of rainbow trout exposed to a magnetic pulse
The code for analysis of the trout retina RNA-seq data.  The results from the following analyses are published in:  

Fitak RR, Schweikert LE, Wheeler BR, Ernst DA, Lohmann KJ, and Johnsen S. (2018) Near absence of differential gene expression in the retina of rainbow trout after exposure to a magnetic pulse: Implications for magnetoreception. *Biology Letters*. **14**(6): 20180209. doi: [10.1098/rsbl.2018.0209](https://doi.org/10.1098/rsbl.2018.0209)  

The raw sequences have been deposited into the NCBI short read archive (SRA): [BioProject PRJNA391176](https://www.ncbi.nlm.nih.gov/bioproject/391176)  

Additional data for this study can be found in:
Fitak R, Schweikert L, Wheeler B, Ernst D, Lohmann K, Johnsen S (2018) Data from: Near absence of differential gene expression in the retina of rainbow trout after exposure to a magnetic pulse: implications for magnetoreception. Dryad Digital Repository. [https://doi.org/10.5061/dryad.24732](https://doi.org/10.5061/dryad.24732)  

The outline and associated code are found below:
1. [Trim the raw sequencing data](./trim-seqs.md)
2. [Download and prepare the reference genome and annotations](./prep-reference.md)
3. [Index the reference genome](./index-reference.md)
4. [Map trimmed reads to reference](./map-reads.md)
   * [Calculate summary mapping statistics](./bam-stats.md)
5. [Verify quantification with STAR and HTSeq](./counts.md)
6. [Differential Expression using DESEQ2](./DESEQ2.md)
7. [Assembly and quantification in CUFFLINKS](./cufflinks-workflow.md)
8. [Analysis of the CUFFLINKS reference annotation-based transcript (RABT) assembly](./Cufflinks-RABT.md)
9. [Various plotting code and functions](./plotting.md)

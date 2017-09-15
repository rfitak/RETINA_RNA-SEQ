# Analysis of the CUFFLINKS reference annotation-based transcript (RABT) assembly
This section will focus on processing and analyzing the novel transfrags (assembled isoforms) produced during the RABT assembly using cufflinks.

```bash
# Get a general summary of the merged assembly
cuffcompare \
   -V \
   -r O_mykiss.annot.gtf \
   -C \
   -G \
   -o merged 
   merged.gtf
```
Summary of the parameters:
- -V: verbose output
- -r: reference annotation in GTF format
- -C: include the "contained" transcripts in the .combined.gtf file
- -G: generic GFF input file(s): do not assume Cufflinks GTF, do not discard any intron-redundant transfrags)
- -o: output file prefix
- the merged gtf file from CUFFMERGE

Summary of the output file 'merged.stats'
```
#= Summary for dataset: ../MERGED/merged.gtf :
#     Query mRNAs :  274039 in  126080 loci  (203849 multi-exon transcripts)
#            (34430 multi-transcript loci, ~2.2 transcripts per locus)
# Reference mRNAs :   46622 in   46602 loci  (43310 multi-exon)
# Super-loci w/ reference transcripts:    33682
#--------------------|   Sn   |  Sp   |  fSn |  fSp  
        Base level: 	100.0	 39.1	  - 	  - 
        Exon level: 	104.9	 47.1	100.0	 58.0
      Intron level: 	100.0	 61.6	100.0	 71.3
Intron chain level: 	100.0	 23.8	100.0	 55.1
  Transcript level: 	100.1	 17.0	100.0	 18.7
       Locus level: 	100.0	 33.4	100.0	 33.4

     Matching intron chains:   48520
              Matching loci:   46602

          Missed exons:       0/323551	(  0.0%)
           Novel exons:  204913/720957	( 28.4%)
        Missed introns:       0/276929	(  0.0%)
         Novel introns:   79556/449303	( 17.7%)
           Missed loci:       0/46602	(  0.0%)
            Novel loci:   80624/126080	( 63.9%)

 Total union super-loci across all input datasets: 126080
```

Now we will make a list of all the isoforms along with its gene id (XLOC\_) and class code
```bash
perl -ne \
   'm/gene_id "([^"]*)";.*transcript_id "([^"]*)".*class_code "([^"]*)"/; print "$1\t$2\t$3\n"' merged.gtf | \
   sort | \
   uniq > merged.class_codes.tsv
```
These class codes are as follows:
- = - Complete match of intron chain
- j - Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript
- o - Generic exonic overlap with a reference transcript
- s - An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)
- u - Unknown, intergenic transcript
- x - Exonic overlap with reference on the opposite strand

Next, we will use R to make a barplot of the class code distribution
```R
library(ggplot2)
a = read.table("merged.class_codes.tsv", sep = "\t", header = F)
names = c("Potentially novel isoform",
   "Potentially novel gene",
   "Exact match to reference",
   "Partial exonic match to reference",
   "Exonic overlap with reference on the opposite strand",
   "Intronic match to reference (possibly from read mapping errors)")
Percent = summary(a$V3)[c(order(summary(a$V3), decreasing = T))] / nrow(a)
data = data.frame(names = names, Percent = Percent)
p = ggplot(data, aes(x = names, y = Percent))
p = p + geom_bar(stat = 'identity')
p = p + coord_flip()
p = p + scale_x_discrete(limits = rev(names))
p = p + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), axis.title.y=element_blank())
p
```

```R
a=scan("transcripts-per-gene.list")
b = ifelse( a > 20, 20, a)
qplot(b, geom="histogram", binwidth = 1) + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + labs(x = "Isoforms per gene", y = "Count")
```

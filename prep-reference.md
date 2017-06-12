# Download reference genome and annotations and merge together with mitochondrial genome and annotations
In this section the general steps are:
- Install NCBI E-utilities
- Download the reference nuclear genome and annotation
- Download the reference mitochondrial genome and annotation
- Merge genomes and annotations together

Installing NCBI E-utilities to my local /bin folder
```
cd ~
   perl -MNet::FTP -e \
      '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
       $ftp->login; $ftp->binary;
       $ftp->get("/entrez/entrezdirect/edirect.zip");'
   unzip -u -q edirect.zip
   rm -rf edirect.zip
   mv edirect/* ~/bin
   rm rf edirect
   ./bin/setup.sh
```

Download the nuclear and mitochondrial genomes and merge together
```
# Make a new folder for the reference files
mkdir TROUT_REF
cd TROUT_REF

# Download the trout reference nuclear genome
wget http://www.genoscope.cns.fr/trout/data/Oncorhynchus_mykiss_chr.fa.gz

# Downloaded trout mitogenome fasta file from NCBI (refseq): NC_001717.1
efetch \
   -db nuccore \
   -format fasta \
   -id NC_001717.1 > mito.fa

# Merge together fasta files
zcat Oncorhynchus_mykiss_chr.fa.gz | \
   cat - mito.fa > Omykiss.genome.fa
```

In this part the nuclear annotations are downloaded and reformatted for use with STAR
```
# Download annotation
wget http://www.genoscope.cns.fr/trout/data/Oncorhynchus_mykiss_chr_annot.gff.gz

# Convert genome GFF3 to GTF format.
gunzip -c Oncorhynchus_mykiss_chr_annot.gff.gz | \
   cut -f9 | \
   sed 's/.* GS/Name=GS/' | \
   paste <(gunzip -c Oncorhynchus_mykiss_chr_annot.gff.gz | \
      cut -f1-8) \
   - > out.gff
gffread out.gff -T -o out.gtf
rm out.gff
mv out.gtf O_mykiss.annot.gtf

# Add "gene_id" to GTF file
cut -f9 O_mykiss.annot.gtf | \
   perl -ne \
   '$_=~ m/\"Name=(.+)\";/; \
   print "gene_id \"$1\"; transcript_id \"$1\";\n"' \
   O_mykiss.annot.gtf | \
   paste <(cut -f1-8 O_mykiss.annot.gtf) - > tmp
mv tmp O_mykiss.annot.gtf

# Remove files not needed
rm Oncorhynchus_mykiss_chr_annot.gff.gz
```
To get the mitochondrial annotations in gff3 format:
- Go to the NCBI page for [NC_001717.1](https://www.ncbi.nlm.nih.gov/nuccore/5835261/)
- Select Send -> Complete Record -> File -> GFF3 -> Create File
- Rename this file to mito.gff3

Finally, merge together with the nuclear annotations
```
# Convert to gtf format
gffread mito.gff3 -T -E --force-exons -o mito.gtf

# Then made a copy of each 'CDS' line and changed to 'exon' using VIM editor

# Change the RNA features to the appropriate format
sed -i "s/transcript_id\( \"rna[0-9]*\"\)/gene_id\1; transcript_id\1/g" mito.gtf

# Merge annotations together
cat O_mykiss.annot.gtf mito.gtf > tmp
mv tmp O_mykiss.annot.gtf
```
The file O_mykiss.annot.gtf is now ready for indexing, and includes both nuclear and mitochondrial annotations.

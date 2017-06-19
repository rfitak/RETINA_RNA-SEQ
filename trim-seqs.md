# This code is for trimming the raw Illumina sequencing data
All that is required is a file listing the base names of the paired sequencing data.  This file is called seqfiles.list.  The following code was actually run as a SLURM script to the cluster, therefore we had to assign each job in the SLURM array to the variable "job".  The raw sequencing files and seqfiles.lsit are in the RAW_SEQS folder.

```bash
# Move into the appropriate folder
cd /work/frr6/RETINA

# Make an output folder for trimmed sequence files
mkdir TRIMMED_SEQS

# Assign SLURM array number to $job
job=$SLURM_ARRAY_TASK_ID

# Get the base name for each job
name=$(sed -n "${job}"p RAW_SEQS/seqfiles.list)

# Assign forward and reverse variable to each sequencing file
forward=RAW_SEQS/${name}.F.fastq.gz
reverse=RAW_SEQS/${name}.R.fastq.gz
```

Now two variables, $forward and $reverse, contain the name of matching forward and reverse (paired) sequencing files.  Now we begin the trimming procedure using the program TRIMMOMATIC v0.36.  Before trimming, we have to assign "trimmomatic" to the actual command to run the program.  This is a simple bash script called "trimmomatic" and is in the ~/bin folder where all programs are stored.
```bash
#!/bin/bash

cmd="java -jar /dscrhome/frr6/PROGRAMS/Trimmomatic-0.35/trimmomatic-0.35.jar"
for i in "$@"
do
cmd="$cmd $i"
done
$cmd
```

Now to trim the sequence files.
```bash
trimmomatic \
   PE \
   -threads 4 \
   -phred33 \
   $forward \
   $reverse \
   /work/frr6/RETINA/TRIMMED_SEQS/${name}_trimmed.F.fastq.gz \
   /work/frr6/RETINA/TRIMMED_SEQS/${name}_trimmed.F.SE.fastq.gz \
   /work/frr6/RETINA/TRIMMED_SEQS/${name}_trimmed.R.fastq.gz \
   /work/frr6/RETINA/TRIMMED_SEQS/${name}_trimmed.R.SE.fastq.gz \
   ILLUMINACLIP:/dscrhome/frr6/PROGRAMS/Trimmomatic-0.35/adapters/all.fa:2:30:7 \
   LEADING:20 \
   TRAILING:20 \
   SLIDINGWINDOW:4:20 \
   AVGQUAL:20 \
   MINLEN:50
```
Here is a desciption of the parameters:
- PE :: the input data are paired-end reads
- -threads 4 :: use four CPUs (threads)
- $forward :: the name of the forward reads file
- $reverse :: the name of the reverse reads file
- name of the trimmed and paired forward reads file
- name of the trimmed and unpaired forward reads file
- name of the trimmed and paired reverse reads file
- name of the trimmed and unpaired reverse reads file
- ILLUMINACLIP:file:2:30:7 :: fasta file of adapter sequences to trim
  * 2 :: seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
  * 30 :: palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
  * 7 :: simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
 - LEADING:20 :: trim bases from start of the read with a Q < 20
 - TRAILING:20 :: trim bases from end of the read with a Q < 20
 - SLIDINGWINDOW:4:20 :: using a 4-base window, remove the last base if average Q < 20
 - AVGQUAL:20 :: remove the entire read if the average quality is < 20
 - MINLEN:50 :: remove reads less than 50 bases after quality trimming
 

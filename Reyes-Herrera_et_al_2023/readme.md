# *Fusarium oxysporum* f. sp. *cubense* TR4 genomics
Here you can find our code used for the analyses performed by David Studholme and Eliana Torres Bedoya included in the following paper:

**Reyes-Herrera, P. H., Torres-Bedoya, E., Lopez-Alvarez, D., Burbano-David, D., Carmona, S. L., Bebber, D. P., Studholme, D. J., Betancourt, M., & Soto-Suarez, M.
(2023)**.
Genome sequence data reveal at least two distinct incursions of the Tropical Race 4 variant of *Fusarium* wilt into South America.
*Phytopathology* **113**: 90–97.
https://doi.org/10.1094/PHYTO-01-22-0034-R.

## Download and index the reference genome sequence
```
wget --no-clobber https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/994/515/GCA_007994515.1_ASM799451v1/GCA_007994515.1_ASM799451v1_genomic.fna.gz

gunzip GCA_007994515.1_ASM799451v1_genomic.fna.gz

bwa index GCA_007994515.1_ASM799451v1_genomic.fna
```

## Trim to remove adaptors and low-quality sequences
### Illumina reads
```
trim_galore -q 20 --paired SRR10054446_1.fq SRR10054446_2.fq SRR10054447_1.fq SRR10054447_2.fq SRR10054448_1.fq SRR10054448_2.fq SRR10054449_1.fq SRR10054449_2.fq SRR10054450_1.fq SRR10054450_2.fq SRR10103605_1.fq SRR10103605_2.fq SRR10125423_1.fq SRR10125423_2.fq SRR10747097_1.fq SRR10747097_2.fq SRR15514269_1.fq SRR15514269_2.fq SRR15514270_1.fq SRR15514270_2.fq SRR15514271_1.fq SRR15514271_2.fq SRR15514272_1.fq SRR15514272_2.fq SRR7226877_1.fq SRR7226877_2.fq SRR7226878_1.fq SRR7226878_2.fq SRR7226879_1.fq SRR7226879_2.fq SRR7226880_1.fq SRR7226880_2.fq SRR7226881_1.fq SRR7226881_2.fq SRR7226882_1.fq SRR7226882_2.fq SRR7226883_1.fq SRR7226883_2.fq SRR9733598_1.fq SRR9733598_2.fq
```
### Nanopore reads
```
./bin/canu -correct -p FAO10428 genomeSize=50m -pacbio FAO10428.fastq

./bin/canu -correct -p FAO15599 genomeSize=50m -pacbio FAO15599.fastq

./bin/canu -correct -p FAO28907 genomeSize=50m -pacbio FAO28907.fastq
```

## Align sequences against against the reference genome
### Illumina reads
```
for i in SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598; do echo $i; bwa mem -t 8 GCA_007994515.1_ASM799451v1_genomic.fna $i_trimmed_1.fq $i_trimmed_2.fq > $i.sam; done
```
### Nanopore reads
```
minimap2 -ax map-ont GCA_007994515.1_ASM799451v1_genomic.fasta FAO10428_trimmed.fastq > FAO10428.sam

minimap2 -ax map-ont GCA_007994515.1_ASM799451v1_genomic.fasta FAO10428_trimmed.fastq > FAO15599.sam

minimap2 -ax map-ont GCA_007994515.1_ASM799451v1_genomic.fastaFAO10428_trimmed.fastq > FAO28907.sam
```

## Convert SAM files to BAM files
```
for i in FAO10428 FAO15599 FAO28907 SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598; do echo $i; samtools view -b $i.sam > $i.bam; done
```

## Index the BAM files
```
for alignmentFile in FAO10428.bam FAO15599.bam FAO28907.bam SRR10054446.bam SRR10054447.bam SRR10054448.bam SRR10054449.bam SRR10054450.bam SRR10103605.bam SRR10125423.bam SRR10747097.bam SRR15514269.bam SRR15514270.bam SRR15514271.bam SRR15514272.bam SRR7226877.bam SRR7226878.bam SRR7226879.bam SRR7226880.bam SRR7226881.bam SRR7226882.bam SRR7226883.bam SRR9733598.bam; do echo $i samtools index $alignmentFile; done
```

## Generate the pileup files from BAM files
```
for alignmentFile in FAO10428 FAO15599 FAO28907 SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598;do echo $i samtools mpileup -f GCA_007994515.1_ASM799451v1_genomic.fna $alignmentFile > $i.pileup; done
```



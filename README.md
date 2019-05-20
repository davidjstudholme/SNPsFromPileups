# SNPsFromPileups
Script for inferring SNPs from pileup files


First, .bam file was generated using Burrows-Wheeler Aligner (BWA) mem version 0.7.15-r1140.

Candidate SNVs were identified using Sequence Alignment/Map tools (SAMtools)/binary call format tools (BCFtools) package, version 1.6, using the following command-lines: 

```
~/samtools-1.6/samtools/samtools mpileup -u -f genome.fasta alignment.bam > alignment.bcf

~/bcftools-1.6/bcftools call -m -v -Ov alignment.bcf  > alignment.vcf
```

The candidate variants were then filtered using the following command line:
```
~/bcftools-1.6/bcftools filter --SnpGap 100 --include '(REF="A" | REF="C" | REF="G" | REF="T") & %QUAL>=35 & MIN(IDV)>=2 & MIN(DP)>=5 & INDEL=0' alignment.vcf > alignment.filtered.vcf
```

This filtering step eliminates indels with low-confidence single-nucleotide variant calls.

It also eliminates candidate SNVs within 10 base pairs of an indel, since alignment artefacts are relatively common in the close vicinity of indels.

Allele frequencies at each SNP site were estimated from frequencies of each base
(adenine (A), cytosine (C), guanine (G) or thymine (T)) among the aligned reads.

Thus, we would expect an allele frequency of close to zero or one for homozygous sites and approximately 0.5 for heterozygous sites in a diploid genome.

The binary alignment/map (BAM)-formatted BWA-mem alignments were converted to pileup format using the samtools mpileup command in SAMtools version 1.6. From the resulting pileup files, we used a custom Perl script (included in Supplementary material) to detect SNPs.

For SNP detection, we considered only sites where depth of coverage by aligned reads was at least 5× for all 17 datasets. 

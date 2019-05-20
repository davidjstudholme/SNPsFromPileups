# SNPsFromPileups
Script for inferring SNPs from pileup files


Burrows-Wheeler Aligner (BWA) mem version 0.7.15-r1140 with default options and parameter values.

Candidate SNVs were identified using Sequence Alignment/Map tools (SAMtools)/binary call format tools (BCFtools) package, version 1.6, using the following command-lines: 

```
samtools mpileup -u -f genome.fasta alignment.bam 4 alignment.bcf and.
bcftools call -m -v –Ov alignment.bcf 4 alignment.vcf
```

The candidate variants were then filtered using the following command line:
```
bcftools filter –SnpGap 100 –include ’(REF¼"A" | REF¼"C" | REF¼"G" | REF¼"T") & %QUAL4¼35 & MIN(IDV)4¼2 & MIN(DP)4¼5 & INDEL¼0’ alignment.vcf 4 align-ment.filtered.vcf
```

This filtering step eliminates indels with low-confidence single-nucleotide variant calls. It also
eliminates candidate SNVs within 10 base pairs of an indel, since alignment artefacts are relatively common in the close vicinity of indels. Allele frequencies at each SNP site were estimated from frequencies of each base (adenine (A),
cytosine (C), guanine (G) or thymine (T)) among the aligned reads. Thus, we would expect an allele frequency of close to zero or one for homozygous sites and approximately 0.5 for heterozygous sites in a diploid genome. The binary alignment/map (BAM)-formatted BWA-mem alignments were con- verted to pileup format using the samtools mpileup command in SAMtools [33] version 1.6 with default options and parameter values. From the resulting pileup files, we used a custom Perl script (included in Supplementary material) to detect SNPs. For SNP detection, we considered only sites where depth of coverage by aligned reads was at least 5× for all 17 datasets. 

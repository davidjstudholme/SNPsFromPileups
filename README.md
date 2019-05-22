# SNPsFromPileups

The purpose of this script is as part of a simple pipeline for calling SNPs and tabulating the allele-frequencies
across a set of individuals.

The scenario is that we have performed genomic re-sequencing on a set of several biological samples.
The genomic sequence reads have been aligned against a reference genome sequence using a tool such as BWA to generate a ```.bam``` file.
We want to identify single-nucleotide positions in the genome that show variation between samples (i.e. SNPs) and we want to estimate the allele frequencies in each sample at each of those SNP sites.
* The samples could be individuals, for example each could be an individual plant or animal. If these individuals are diploid, then we expect that the allele frequencies will always be 0, 0.5 or 1. If the individuals are triploid then we would expect 0, 0.33, 0.67, or 1. And so on for other levels of ploidy. There is a discrete set of possible values for allele frequency.
* Alternatively, the samples could be populations rather than individuals; for example, they could be microbial cultures containing many millions of individuals of various different genotypes. In this case, the allele frequency depends on the relative abundance of each genotype within the population. There is a continuous set of possible values for allele frequency.

So, let's assume that we have one ```.bam``` file for each sample (genomic sequence reads aligned against reference genome ```genome.fasta```).
The first step is to identify candidate SNPs using SAMtools version 1.6 and BCFtools version 1.6: 

```
foreach my alignmentFile *.bam; do samtools mpileup -u -f genome.fasta $alignmentFile > $alignmentFile.bcf; done

foreach my alignmentFile *.bam; do bcftools call -m -v -Ov $alignmentFile.bcf > $alignmentFile.vcf

```

So, now we have a set of .vcf files containing candidate SNPs. Let's filter these to keep only high-confidence ones:
```

~/bcftools-1.6/bcftools filter --SnpGap 100 --include '(REF="A" | REF="C" | REF="G" | REF="T") & %QUAL>=35 & MIN(IDV)>=2 & MIN(DP)>=5 & INDEL=0' alignment.vcf > alignment.filtered.vcf
```

This filtering step eliminates indels with low-confidence single-nucleotide variant calls.

It also eliminates candidate SNVs within 10 base pairs of an indel, since alignment artefacts are relatively common in the close vicinity of indels.

Allele frequencies at each SNP site were estimated from frequencies of each base
(adenine (A), cytosine (C), guanine (G) or thymine (T)) among the aligned reads.

Thus, we would expect an allele frequency of close to zero or one for homozygous sites and approximately 0.5 for heterozygous sites in a diploid genome.

The BAM alignments were converted to pileup format using the samtools mpileup command in SAMtools version 1.6. From the resulting pileup files, we used a custom Perl script (get_snps_from_pileups.pl) to detect SNPs.

```
perl get_snps_from_pileups.pl 10 alignment.filtered.vcf *.pileup > snps.csv
```

For SNP detection, we considered only sites where depth of coverage by aligned reads was at least 10× for all datasets. 

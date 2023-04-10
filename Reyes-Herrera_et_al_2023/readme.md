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

### Download the Illumina FASTQ files
Use NCBI' SRA tools: https://www.metagenomics.wiki/tools/short-read/ncbi-sra-file-format
```
fasterq-dump --split-3 SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598 -p
```
### Download Oxford Nanopore FASTQ files
```
fasterq-dump SRR16568737 SRR16568738 SRR16568739 -p
```

## Trim to remove adaptors and low-quality sequences
### Illumina reads
```
for i in SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598
do echo $i
trim_galore -q 30 --paired $i"_1.fastq" $i"_2.fastq" && gzip $i"_1_val_1.fq" SRR10054446_2_val_2.fq && rm $i"_1.fastq" $i"_2.fastq"
done   
```
### Nanopore reads
Using: https://github.com/wdecoster/chopper
```
for i in  SRR16568737 SRR16568738 SRR16568739; do echo $i 
chopper -q 10 -l 500 < $i.fastq > $i.chopper.fq
if [ -s $i.fastq ]; then
  rm $i.fastq
fi
done
```

## Align sequences against against the reference genome to generate indexed BAM files
Using SAMtools: 
**Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup** (2009). 
The Sequence Alignment/Map format and SAMtools. 
*Bioinformatics*, **25** 2078–2079. https://doi.org/10.1093/bioinformatics/btp352.
### Illumina reads
Using BWA:
**Li, H., & Durbin, R.** (2009).
Fast and accurate short read alignment with Burrows-Wheeler transform. 
*Bioinformatics*, **25** 1754–1760. https://doi.org/10.1093/bioinformatics/btp324.
```
for i in SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598
do echo $i
if ! [ -s $i.sorted.bam ]; then
  bwa mem -t 4 GCA_007994515.1_ASM799451v1_genomic.fna SRR10054446_1_val_1.fq.gz SRR10054446_2_val_2.fq.gz > $i.sam
  samtools view -h -b -q 1 $i.sam > $i.bam && rm $i.sam
  samtools sort $i.bam -o $i.sorted.bam && rm $i.bam
  samtools index $i.sorted.bam 
fi
done
```
### Nanopore reads
Using Minimap2: 
**Li H.** (2018). 
Minimap2: pairwise alignment for nucleotide sequences. 
*Bioinformatics*, **34**, 3094–3100. 
https://doi.org/10.1093/bioinformatics/bty191.
```
for i in SRR16568737 SRR16568738 SRR16568739
do echo $i
if ! [ -s $i.sorted.bam ]; then
  minimap2 -ax map-ont GCA_007994515.1_ASM799451v1_genomic.fna $i.chopper.fq > $i.sam
  samtools view -h -b -q 1 $i.sam > $i.bam && rm $i.sam
  samtools sort $i.bam -o $i.sorted.bam && rm $i.bam
  samtools index $i.sorted.bam 
fi
done
```

## Generate the pileup files from BAM files
### Illumina data
```
for i in SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598
do echo $i 
if ! [ -s $i.pileup ]; then
  samtools mpileup -q 1 -f GCA_007994515.1_ASM799451v1_genomic.fna $i.sorted.bam > $i.pileup
fi
done
```
### Nanopore data
```
for i in SRR16568737 SRR16568738 SRR16568739
do echo $i 
if ! [ -s $i.pileup ]; then
  samtools mpileup -B -q 1 -f GCA_007994515.1_ASM799451v1_genomic.fna $i.sorted.bam > $i.pileup
fi
done
```

## Generate VCF using Pilon
**Walker, B. J., Abeel, T., Shea, T., Priest, M., Abouelliel, A., Sakthikumar, S., Cuomo, C. A., Zeng, Q., Wortman, J., Young, S. K., & Earl, A. M.**
(2014).
Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement. 
*PLoS One*, **9**: e112963.
https://doi.org/10.1371/journal.pone.0112963.
```
for i in SRR10054446 SRR10054447 SRR10054448 SRR10054449 SRR10054450 SRR10103605 SRR10125423 SRR10747097 SRR15514269 SRR15514270 SRR15514271 SRR15514272 SRR7226877 SRR7226878 SRR7226879 SRR7226880 SRR7226881 SRR7226882 SRR7226883 SRR9733598 SRR16568737 SRR16568738 SRR16568739
do echo $i
if ! [ -e pilon_$i ]; then
  pilon --genome GCA_007994515.1_ASM799451v1_genomic.fna --bam $i.sorted.bam --output pilon_$i --vcf
fi
done
```

## Filter the VCF files with bcftools
**Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H.**
(2021).
Twelve years of SAMtools and BCFtools.
*GigaScience* **10**: giab008.
https://doi.org/10.1093/gigascience/giab008.

```
for alignmentFile in FAO10428.bam FAO15599.bam FAO28907.bam SRR10054446.bam SRR10054447.bam SRR10054448.bam SRR10054449.bam SRR10054450.bam SRR10103605.bam SRR10125423.bam SRR10747097.bam SRR15514269.bam SRR15514270.bam SRR15514271.bam SRR15514272.bam SRR7226877.bam SRR7226878.bam SRR7226879.bam SRR7226880.bam SRR7226881.bam SRR7226882.bam SRR7226883.bam SRR9733598.bam; do echo $alignmentFile; bcftools filter --include '(REF="A" | REF="C" | REF="G" | REF="T") & (ALT="A" | ALT="C" | ALT="G" | ALT="T")' pilon_$alignmentFile.vcf > $alignmentFile.filtered.vcf; done
```

## Get the SNP-calling scripts from GitHub
```
git clone https://github.com/davidjstudholme/SNPsFromPileups.git
```

## Perform SNP-calling from pileup files.
```
perl SNPsFromPileups/get_snps_from_pileups_small_genome.pl 10 FAO10428.bam.filtered.vcf FAO15599.bam.filtered.vcf FAO28907.bam.filtered.vcf SRR10054446.bam.filtered.vcf SRR10054447.bam.filtered.vcf SRR10054448.bam.filtered.vcf SRR10054449.bam.filtered.vcf SRR10054450.bam.filtered.vcf SRR10103605.bam.filtered.vcf SRR10125423.bam.filtered.vcf SRR10747097.bam.filtered.vcf SRR15514269.bam.filtered.vcf SRR15514270.bam.filtered.vcf SRR15514271.bam.filtered.vcf SRR15514272.bam.filtered.vcf SRR7226877.bam.filtered.vcf SRR7226878.bam.filtered.vcf SRR7226879.bam.filtered.vcf SRR7226880.bam.filtered.vcf SRR7226881.bam.filtered.vcf SRR7226882.bam.filtered.vcf SRR7226883.bam.filtered.vcf SRR9733598.bam.filtered.vcf FAO10428.bam.pileup FAO15599.bam.pileup FAO28907.bam.pileup SRR10054446.bam.pileup SRR10054447.bam.pileup SRR10054448.bam.pileup SRR10054449.bam.pileup SRR10054450.bam.pileup SRR10103605.bam.pileup SRR10125423.bam.pileup SRR10747097.bam.pileup SRR15514269.bam.pileup SRR15514270.bam.pileup SRR15514271.bam.pileup SRR15514272.bam.pileup SRR7226877.bam.pileup SRR7226878.bam.pileup SRR7226879.bam.pileup SRR7226880.bam.pileup SRR7226881.bam.pileup SRR7226882.bam.pileup SRR7226883.bam.pileup SRR9733598.bam.pileup > snps_final.csv
```

## Convert the SNPs into Nexus format for input into iqtree
```
perl SNPsFromPileups/get_haplotypes_and_aligned_fasta_from_csv.pl snps_final.csv
```

## Perform phylogenetic analysis using iqtree
**Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q.** (2015).
IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. 
*Molecular Biology and Evolution* **32**: 268–274.
https://doi.org/10.1093/molbev/msu300.
```
iqtree -s snps_final.csv.haplotype.nex -m GTR+ASC
```

## Perform bootstrapping
```
iqtree -s snps_final.csv.haplotype.nex.uniqueseq.phy -m TIM2+I+G -bb 1000
```



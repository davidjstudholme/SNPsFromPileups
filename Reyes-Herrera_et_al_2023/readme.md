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

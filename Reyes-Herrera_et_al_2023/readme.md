# *Fusarium oxysporum* f. sp. *cubense* TR4 genomics
Here you can find our code used for the analyses performed by David Studholme and Eliana Torres Bedoya included in the following paper:

**Reyes-Herrera, P. H., Torres-Bedoya, E., Lopez-Alvarez, D., Burbano-David, D., Carmona, S. L., Bebber, D. P., Studholme, D. J., Betancourt, M., & Soto-Suarez, M.
(2023)**.
Genome sequence data reveal at least two distinct incursions of the Tropical Race 4 variant of *Fusarium* wilt into South America.
*Phytopathology* **113**: 90–97.
https://doi.org/10.1094/PHYTO-01-22-0034-R.

### Download and index the reference genome sequence
```
wget --no-clobber https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/994/515/GCA_007994515.1_ASM799451v1/GCA_007994515.1_ASM799451v1_genomic.fna.gz

gunzip GCA_007994515.1_ASM799451v1_genomic.fna.gz

bwa index GCA_007994515.1_ASM799451v1_genomic.fna
```

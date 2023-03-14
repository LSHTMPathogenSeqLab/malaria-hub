# Env preparation

```
Tools: bcftools
R packages: moimix, SeqArray, optparse, readr, dplyr
```
# Input data

## Filter genotyped merged VCF
Fws metric needs to be calculated based on the relevant SNPs (coding region, bi-allelic, with high MAF). Specification of thresholds should be adjusted for any specific set.

### Include only coding regions

Extract exon information from GFF file and create tab-separated bed file with chr, start, end columns. Use `bcftools view -R <coding.regions.bed>` option.

### Bi-allelic SNPs

If passed through `filter_merged_vcf.py` this step has been already done (_bi_ suffix).  If not use `bcftools view -m2 -M2` option.

### MAF

Use `bcftools view -q` option foe filtering. Adjust threshold to not remove to many SNPs.

-----------------------------------------
## Split metadata to populations
Calculations of Fws are done per population, therefore it is needed to create subsets.
Metadata needs to contain only samples that passed filtering (filter_merged_vcf.py_) and are in final VCF with filtered SNPs. Aim is to create sub file with sample names for each population.

```{r}
# Example how to create sub files in R
library(readr)
library(dplyr)

met <- readr::read_tsv('pf_metadata.tsv') # Load metadata
suffix <- "pf_meta" # Add additional suffix after country name

met %>%
  group_by(country) %>%
  dplyr::select(sra_run, country) %>%
  group_walk(~ write.table(.x, paste0(.y$country, suffix, ".txt"), quote = F, row.names = F, col.names = F))
```
-----------------------------

## Split genotyped merged VCF
Based on sub files create VCF subsets per population
```
ls *<suffix>.txt > population_files_list.txt

cat population_files_list.txt | xargs -I {} sh -c 'bcftools view -c1 -S {} -Oz --threads 20 -o {}.genotyped.vcf.gz <input_vcf>'
```
---------------------------------------------

# Run moimix for each population
Requires installation of `SeqArray`, `moimix` and `optparse` packages. Script transforms input VCF file to GDS file (`SeqArray::seqVCF2GDS`). Based on created GDS file calculates Fws metric (`moimix::getFws()`). Output in form of two column tab-separated file for each population `<population>_moi_fws.tsv`.

## Single population
```{bash}
Rscript ~/software/malaria-hub/moi/calculate_fws.R -d <workdir> -f <population_vcf>'
```
## Multiple populations
```
ls *.genotyped.vcf.gz | sed 's/.genotyped.vcf.gz//g' > population_names_list.txt

cat population_names_list.txt | xargs -I {} sh -c 'Rscript ~/software/malaria-hub/moi/calculate_fws.R -d <workdir> -f {}.genotyped.vcf.gz'
```

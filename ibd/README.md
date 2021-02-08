# Analysis in few words

IBD analysis is performed in two-step manner: input files preparation and running [hmmIBD](https://malariajournal.biomedcentral.com/articles/10.1186/s12936-018-2349-7) then parsing and annotating results. Scripts for populations are designed to be run in command line by specifying arguments. If calling script from R level is prefered, replace variables defined in top-level and omit optparse chunk.

Firstly, input data is parsed, filtered and recoded to match requirements.
Script processes one population at a time (take only populations that have sufficient sample size). Preparing data can be done on few resolution levels region/country/site/ancestry. Population is subsetted based on metadata description.
Following, hmmIBD should produce two files `_out.hmm.txt` and `_out.hmm_fract.txt`. Description of output files can be studied [here](https://github.com/glipsnort/hmmIBD). Log from hmmIBD run is stored in `hmmIBD_run_<pop_name>.log`.  

Filtering steps:
* select snps with MAF > th (default 0.01)
* select samples with Fws > th (default 0.95)
* recode data to 0 (REF), 1 (ALT), -1 (MISSING)
* !! IMPORTANT !! all missing calls are replaced to reference and all mixed calls to alternative

Secondly, files produced by hmmIBD are combined across selected categories and annotated. In order to calculate cumulative fraction across genome (per population) sliding window analysis was used. Specigy window size with `--window_size` argument. Example plots and visualization can be explored in `vis_hmmIBD_results.R`.

# Env preparation

### Create conda env
```{bash}
conda create --name ibd r-base r-dplyr r-data.table r-optparse r-stringr r-ggplot2 r-ggrepel r-tidyr r-readr --channel conda-forge
conda activate ibd
```

### Install hmmIBD in ~/software dir
```{bash}
mkdir ~/software
cd software
git clone https://github.com/glipsnort/hmmIBD.git
cc -o hmmIBD -O3 -Wall hmmIBD.c -lm
```

### Clone repository to ~/software dir
```{bash}
cd ~/software
git clone https://github.com/LSHTMPathogenSeqLab/malaria-hub.git
```

# Input data

* binary matrix - matrix needs to be after standard filtering (__bi-allelic__ snps only). Selection analysis will be done only for nucleus chromosomes. Support matrix without apicoplast and mitochondrion chromosomes or use `--remove_chr` argument to exlude them in analysis (REQUIRED). In order to transform chromosome numbering, as default `<clone>_<number>_<version>` (Pf3D7_01_v3) pattern is used. If your species organism is not matching standard format use `--regex_chr` and `--regex_groupid` to specify your own pattern.

* metadata - tab-separated file with description of every sample. Requires field with Fws scores, country, region or site category. Column names describing sample id (`--label-id`), country and region `--label_category`, fws score (`--label_fws`) can be defined on arguments level.

* annotation - tab-separated file with description of every snp. Required columns __Chr, Pos, Ref, Alt_1, Gene_name_1__. Gene naming need to be in accordance. Can be prepared with SnpEff/CSQ.

* gene/product file - tab-separated file with information about gene/product. Required columns __chr, pos_start, pos_end, gene_id, product, gene_name__ Used to annotate candidate regions. Can be extracted per species from PlasmoDB. Go to _Genes_ > _Annotation,curation and identifiers_ > _Updated annotation at GeneDB_ > _Species filter_ >  _Apply_ > _Download table_.

# Arguments

## run_hmmIBD_per_category.R ##

* `-d --workdir` - working directory to save results
* `-b --binary_matrix` - absolute location of binary matrix
* `-m --metadata` - absolute location of metadata. Requires to have all columns specified in label_id, label_category, label_fws
* `-c --category` - name of country or region ex. Peru or Oceania that analysus will be done for
* `--label_category` - column name in metadata file with category name ex. country/region
* `--label_fws` - columna name in metadata file with Fws score
* `--fws_th` - threshold for Fws score (default 0.95)
* `--label_id` - column name in metadata with sample id
* `--maf` - threshold for MAF (default 0.01)
* `--na_char` - specify missing calls character (default NA)
* `--threads` - optional argument to specife threads usage
* `--remove_chr` - field to specify non nuclear chromosomes that need to be removed
* `--regex_chr` - regex pattern to detect chromosome numbering
* `--regex_groupid` - group id for regex with numbering

## parse_annotate_hmmIBD_results.R ##

* `-d --workdir` - working directory to save results
* `--list_category` - file with category names per line
* `-l --legend` - __ibd_matrix_hap_leg.tsv__ legend file (produced with previous script)
* `--gene_product` - absolute location of annotation file. Requires to have chr, pos_start, pos_end, gene_id, product, gene_name fields. Used to annotate candidate regions.
* `-r --ref_index` - Reference index (.fai)
* `--window_size` - Window size for genome wide fraction calculations (default 10000)
* `--quantile_cutoff` - 
* `--suffix` - suffix for output files

* `-t --threads` - optional argument to specife threads usage
* `--remove_chr` - field to specify non nuclear chromosomes that need to be removed
* `--regex_chr` - regex pattern to detect chromosome numbering
* `--regex_groupid` - group id for regex with numbering

# Run scripts

## Part I
### Single population
```{bash}
Rscript ~/software/malaria-hub/ibd/run_hmmIBD_per_category.R \
-d <workdir> \
-b <binary_matrix> \
-m <metadata> \
-c DRC \
--label_category region \
--label_fws fws \
--fws_th 0.95 \
--label_id sra_run \
--maf 0.01 \
--na_char NA \
--remove_chr Pf3D7_API_v3,Pf_M76611

```
### Multiple population can be run with xrags or parallel
```{bash}
cat country_list.txt | xargs -I {} -P 2 sh -c 'Rscript ~/software/malaria-hub/ibd/run_hmmIBD_per_category.R \
-d <workdir> \
-b <binary_matrix> \
-m <metadata> \
-c {} \
--label_category region \
--label_fws fws \
--fws_th 0.95 \
--label_id sra_run \
--maf 0.01 \
--na_char NA \
--remove_chr Pf3D7_API_v3,Pf_M76611'

```
## Part II
```{bash}
Rscript ~/software/malaria-hub/ibd/parse_annotate_hmmIBD_results.R\
-d <workdir>\
--list_category country_list.txt \
--legend ibd_matrix_hap_leg.tsv \
--gene_product <gene_product>
--ref_index Pfalciparum.genome.fasta.fai \
--maf 0.01 \
--windwo_size 10000 \
--quantile_cutoff 0.95 \
--remove_chr Pf3D7_API_v3,Pf_M76611
```

# Output data

## run_hmmIBD_per_category.R ##

_Results per category_

* `ibd_matrix_hap_leg.tsv` - matrix legend (same for each iteration)
* `ibd_matrix_hap_<category>.tsv` - recoded matrix per category
* `hmmIBD_<category>_maf<th_maf>.txt` - input file for hmmIBD (after filtering)
* `hmmIBD_run_<category>.log` - logfile to observe output from hmmIBD
* `hmmIBD_<category>_maf<th_maf>_out.hmm.txt` - 
* `hmmIBD_<category>_maf<th_maf>_out.hmm.fract.txt` - 

## parse_annotate_hmmIBD_results.R ##

* `<suffix>_hmmIBD_fraction_results_combined.tsv` - table with pairwise fraction collected for all specified categories
* `<suffix>_hmmIBD_ibd_results_combined.tsv` - table with re-calculated genome-wide IBD fractions collected for all specified categories
* `<suffix>_hmmIBD_ibd_annotated_q0.95.tsv` - table with annotated candidate regions collected for all specified categories

## Visualization
Example visualization script in `vis_hmmIBD_results.R`

![Example plots](https://github.com/LSHTMPathogenSeqLab/malaria-hub/blob/hmmibd/ibd/example_plots.png?raw=true)
# Analysis in few words

Selection analysis is performed in two-step manner with usage of `rehh` package. Scripts for populationsmain are designed to be run in command line by specifying arguments. If calling script from R level is prefered, replace variables defined in top-level and omit optparse chunk.
To familiarize with analysis recommended to read vignette [here.](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html)

--------

Firstly, parsing, filtering and recoding input to match requirements of rehh package. Process populations that have sufficient size. Script processes one population at a time. Preparing data can be done on few resolution levels region/country/site. Population is collected based on metadata description. Results in form of tab-separated file with calculated IHH, IES, INES metrics. 

__Filtering step:__
* selecting snps with MAF > th (_default_ 0.01)
* selecting samples with Fws > th (_default_ 0.95)
* filtering missing data in `data2haplohh()` function in samples (min_perc_geno.hap) and snps (min_perc_geno.mrk).
 Can be modified with arg `rehh_min_perc_hap` & `rehh_min_perc_mrk`) (_default_ min_perc_geno.hap = 80 & min_perc_geno.mrk = 70)
* __!! IMPORTANT !!__ As default missing or mixed calls are not modified. If data is of high quality replacement can be considered. To proceed with forced data recoding add `--forced_recode` argument (missing calls -> reference and mixed calls -> alternative). To only modify mixed calls -> alternative add `--forced_mixed` argument.

-----

Secondly, using outcomes from first part calculating standard metrics __iHS__ (per population), __Rsb__ and __XP-EHH__ (pairwise comparisons between populations). Metric results are represented in form of manhattan plots, tables with highly significant selection snps and candidate regions. All results are annotated. Collection of results are done for all populations at once.

__Filtering steps:__
* ihh2ihs(min_maf = 0.0, freqbin = 0.05)
* iHS (pvalue >= `ihs_th`), Rsb (pvalue >= `--rsb_th`), XP-EHH (pvalue >= `--xpehh_th`) 
* calculate_candidate_regions(threshold = `--ihs_th`,
                              pval = TRUE,
                              window_size = 2E4,
                              overlap = 1E4,
                              min_n_extr_mrk = 2)

-----

__Functions used from rehh package:__

```{r}
rehh::data2haplohh()
rehh::scan_hh()
rehh::ihh2ihs()
rehh::ines2rsb()
rehh::ies2xpehh()
rehh::calc_candidate_regions()
```

# Env preparation

### Create conda env
```{bash}
conda create -n selection r-base r-dplyr r-data.table r-optparse r-stringr r-ggplot2 r-ggrepel r-tidyr --channel conda-forge
conda activate selection
```

### Install rehh package in R
```{r}
install.packages('rehh')
```

### Clone repository to ~/software dir
```{bash}
cd ~/software
git clone https://github.com/LSHTMPathogenSeqLab/malaria-hub.git
```

# Input data

* binary matrix - matrix needs to be after standard filtering (__bi-allelic__ snps only). Selection analysis will be done only for nucleus chromosomes. Support matrix without apicoplast and mitochondrion chromosomes or use `--remove_chr` argument to exlude them in analysis (REQUIRED)

* metadata - tab-separated file with description of every sample. Requires field with Fws scores, country, region or site category. Column names describing sample id `--label-id`, country and region `--label_category`, fws score `--label_fws` can be defined on arguments level.

* annotation - tab-separated file with description of every snp. Required columns __chr, pos, ref, alt, Gene_name__. Gene naming need to be in accordance. Can be prepared with SnpEff or CSQ.

* gene/product file - tab-separated file with information about gene/product. Required columns __chr, pos_start, pos_end, gene_id, gene_product, gene_name__ Used to annotate candidate regions. Can be extracted per species from PlasmoDB. Go to _Genes_ > _Annotation,curation and identifiers_ > _Updated annotation at GeneDB_ -> Species filter >  Apply > Download table.

# Arguments

## prepare_input_rehh_per_category.R ##

* `-d --workdir` - working directory to save results
* `-b --binary_matrix` - absolute location of binary matrix
* `-m --metadata` - absolute location of metadata. Requires to have all columns specified in label_id, label_category, label_fws
* `-a --annotation` - absolute location of annotation file. Requires to have chr, pos, ref, alt, Gene_name fields.
* `-c --category` - name of country or region ex. Peru or Oceania that analysis will be done for
* `--label_category` - column name in metadata file with category name ex. country/region
* `--label_fws` - columna name in metadata file with Fws score
* `--fws_th` - threshold for Fws score (_default_ 0.95)
* `--label_id` - column name in metadata with sample id
* `--maf` - threshold for MAF (_default_ 0.01)
* `--rehh_min_perc_hap` - Threshold data2haplohh % of missing data for sample  (_default_ 80)
* `--rehh_min_perc_mrk` - Threshold data2haplohh % of missing data for snp (_default_ 70)
* `--na_char` - Specify character for missing calls (_default_ NA)
* `--forced_recode` - Recode missing calls to reference and mixed calls to alternative (_optional_)
* `--forced_mixed` - Recode mixed calls to alternative (_optional_)
* `--remove_chr` - field to specify non nuclear chromosomes that need to be removed
* `--regex_chr` - regex pattern to detect chromosome numbering
* `--regex_groupid` - group id for regex with numbering
* `--threads` - argument to specife threads usage (_optional_)

## calculate_rehh_metrics.R ##

* `-d --workdir` - working directory to save results
* `-p --prefix` - prefix name with scan_hh objects (default scanned_haplotypes)
* `--list_category` - file with category names per line
* `--annotation` - absolute location of annotation file. Requires to have chr, pos, ref, alt, Gene_name fields.
* `--gene_product` - absolute location of gene product file. Requires to have chr, pos_start, pos_end, gene_id, product, gene_name fields. Used to annotate candidate regions.
* `--remove_chr` - field to specify non nuclear chromosomes that need to be removed
* `--regex_chr` - regex pattern to detect chromosome numbering
* `--regex_groupid` - group id for regex with numbering
* `--ihs_th` - log p-value threshold for iHS metric (_default_ 4)
* `--rsb_th` - log p-value threshold for Rsb metric (_default_ 5)
* `--xpehh_th` - log p-value threshold for XP-EHH metric (_default_ 5)
* `--threads` - argument to specife threads usage (_optional_)

# Run scripts

## Part I
### Single population
```{bash}
Rscript ~/software/malaria-hub/selection/prepare_input_rehh_per_category.R \
-d <workdir> \
-b <binary_matrix> \
--remove_chr Pf3D7_API_v3,Pf_M76611 \
-m <metadata> \
--annotation <annotation> \
-c West_Africa \
--label_category region \
--label_fws fws \
--fws_th 0.95 \
--label_id sra_run \
--forced_mixed
```
### Multiple population can be run with xrags or parallel
```{bash}
cat country_list.txt | xargs -I {} -P 2 sh -c 'Rscript ~/software/malaria-hub/selection/prepare_input_rehh_per_category.R \
-d <workdir> \
-b <binary_matrix> \
--remove_chr Pf3D7_API_v3,Pf_M76611 \
-m <metadata> \
--annotation <annotation> \
-c {} \
--label_category country \
--label_fws fws \
--fws_th 0.95 \
--label_id sra_run \
--forced_mixed'

```
## Part II
```{bash}
Rscript ~/software/malaria-hub/selection/calculate_rehh_metrics.R\
-d <workdir>\
--prefix scanned_haplotypes \
--remove_chr Pf3D7_API_v3,Pf_M76611 \
--list_category region_list.txt \
--annotation <annotation> \
--gene_product <gene_product>
```

# Output data

## prepare_input_rehh_per_category.R ##

_Results per category_

* `hap_chr_*_<category>` - haplotypes per chromosome
* `snp.info.inp.<category>` - list of input snps
* `<category>.log` - logfile to observe output from haplotype processing data2haplohh() & scan_hh()
* `scanned_haplotypes_<category>.tsv` - tab-separated file for snps passing filtering process with calculation of IHH_A, IHH_D, IES, INES. Input for second script.

## calculate_rehh_metrics.R ##

* `plots_<category>.pdf` - collection of iHS, Rsb, XP-EHH plots per catgeory. SNPs that pass threshold (PVALUE > 4 or 5) are highlighted and genes that have at least two hits are labelled
__!!NOTE!!__ only chromosomes where snps available are displayed
* `plots_<category>.rds` - separate plots for each category and metric in `.rds` format. Can be loaded and edited in R with `readRDS` function.
* `high_<metric>_all_categories.tsv` - table with signigicant snps for all categories
* `cr_<metric>_all_categories_annot.tsv` - table with annotated candidate regions for all categories
* `<category>_metrics.log` - output from metrics calculations








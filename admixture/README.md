# ADMIXTURE analysis in few steps

## As starting point have genotyped VCF for all populations after standard filtering steps

## 1. Additional filtering of SNPs in VCF
To select applicable filters check methods section from published work with similar species/dataset.
It is also advised to check how the filters influence admixture results.
Filters that can be considered:
* LD - can be done with `--indep-pairwise` option in plink
* MAF - can be done with `-q` option in bcftools

## 2. Convert VCF to BED

```
plink --vcf <vcf_file> --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out <prefix>
```

## 3. Convert chromosome names to numbers in .bim file
### Example command for P. vivax genome
NOTE: This command will automatically make changes in the `<prefix>.bim` file.
If you prefer original version to remain, remove flag `i` from command and `>` outcome to new file

```
sed -ie 's/PvP01_//g; s/_v1//g; s/^0//g; s/API/15/g; s/MIT/16/g' <prefix>.bim
```

## 4. Run admixture
Admixture can be installed with conda. Documentation can be found [here](https://vcru.wisc.edu/simonlab/bioinformatics/programs/admixture/admixture-manual.pdf)

NOTE:
* For Plasmodium we use `--haploid` option, it might not be required for other species
* It is advised to run analysis with different seeds (flag `-s`) to check reproducibility of the results

### Example command for Plasmodium single K=1 run
```
admixture --cv=10 -j8 --haploid="*" -s 12345 <prefix>.bed 1 | tee log1.cv10.haploid.seed12345.out
```
### Example command for Plasmodium for multiple K runs
NOTE:
* Specify all K in file `K_runs.txt` from 1 to 10 and execute command with xargs
```
cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j8 --haploid="*" -s 12345 <prefix>.bed {} | tee log{}.cv10.haploid.seed12345.out'
```

## 5. Inspect CV error it logs (documentation page 4)
Plot CV errors and inspect infleciton point
```
grep -h CV *out
```

## 6. Visualize in admixture-like plot (documentation page 6)
Script to plot admixture-like barplots is in `malaria-hub/admixture/generate_admix_barplot.R`. It produces graph in `tiff` format for __region__ and __country__ and optionally for __site__. Graphs can single K or multiple Ks on one plot.

### Prepare ENV
```bash
conda create -n radmix r-essentials r-base
R
```
```r
install.packages(c("unikn", "countrycode", "optparse"))
```


Inputs:
* tab-separated metadata file with sample identifier, region, coutnry and optionaly site (`--label_site`)
* output files from ADMIXTURE runs `<prefix>.<K>.Q`
* output file from PLINK with sample order `<prefix>.nosex`
* single K or comma-separated Ks (multiple give just instant comparison but they are sorted by the first K)
* If you want specific order of regions to display use `region_order`
* filter population based on N `<filter_N>`
* select specific country `<select_country>`


```
Options:
	-d CHARACTER, --workdir=CHARACTER
		Specify main directory

	--prefix=CHARACTER
		Prefix of admixture output files <prefix>.Q

	--filename=CHARACTER
		Filename for plots

	-k CHARACTER, --kval=CHARACTER
		Comma-separated K values

	-m CHARACTER, --metadata=CHARACTER
		Metadata filepath with sample, region, country (optionally) site information

	--region_order=CHARACTER
		Region order

	-s CHARACTER, --select_country=CHARACTER
		Filter per country

	-f INTEGER, --filter_N=INTEGER
		Select population with N > [default: 20]

	--label_region=CHARACTER
		Region label name in metadata file

	--label_country=CHARACTER
		Country label name in metadata file

	--label_site=CHARACTER
		Region label name in metadata file

	--label_id=CHARACTER
		Category label name in metadata file

	--country_code=LOGICAL
		Country codes

	--axisx_angle=NUMERIC
		Rotation angle for X axis labels

	-h, --help
		Show this help message and exit
```

## Example command for one K:

```bash
Rscript ~/software/malaria_hub/admixture/generate_admix_barplot.R \
-d <workdir> \
--prefix <prefix> \
--kval 7
-m <metadata> \
--filter_N 20 \
--label_id sample_id \
--label_region region \
--label_country country \
--label_id sample_id \
--country_code TRUE
```

## Example command for multi Ks (all in one plot but sorted based on the first):

```bash
Rscript ~/software/malaria_hub/admixture/generate_admix_barplot.R \
-d <workdir> \
--prefix <prefix> \
--kval 5,6,7 
-m <metadata> \
--filter_N 20 \
--label_id sample_id \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE
```

## Example to select few countries:
```
Rscript ~/software/malaria_hub/admixture/generate_admix_barplot.R \
-d <workdir> \
--prefix <prefix> \
--kval7 
-m <metadata> \
--filter_N 20 \
--label_id sample_id \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE \
--select_country Thailand,Vietnam
```

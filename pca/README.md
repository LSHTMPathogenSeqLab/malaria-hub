# Example PCA analysis


## Input files / Dependencis

R (amap, ggplot2), Plink

* VCF file bi-allelic after all filtering steps
* Metadata with ID, country and region indication

## Run PLINK
Use plink to calculate distances. Output files
* <prefix>.dist - distance matrix NxN dimenstion (N is number of samples)
* <prefix>.dist.id - order of samples

```
plink --vcf <vcf-file> --distance square --double-id --allow-extra-chr --out <prefix>
```

## Work in R based on example script
Work through script `example_pca_plink.R`, specify workdir, prefix and metadata matching your input data.

Figures are saved to `.png` files.


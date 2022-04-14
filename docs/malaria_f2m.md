 # fastq2vcf

Process raw samples with fastq2matrix pipeline with gVCF files as output and BQSR correction.
 
 # Create or extend GenomicsDB
 Create DB and import new isolates to store together information from multiple VCFs. Pass the DB name (`--prefix`) and list of samples you want to include (`--sample-file`) along with a directory of your gVCFs (`--vcf-dir`) and `--vcf-extension` (default `g.vcf.gz`). Parameter `--num-genome-chunks` can be the same as default - information will be stored in `<database_name>.dbconf.json` file at the end of the run. A successful run produces a complete set of `genomics_db` folders that should exist for each genomic subregion (number of chromosomes x number of genome chunks).

 Creating and extending the database on top-level uses the same command. To extend the database make sure you use the correct database name (`--prefix`) and run command from the same directory as the existing `json` file and `genomics_db` folders. Assure that to-be-imported isolates have no overlap with previously uploaded samples that exist already in DB.
[Read more here](https://gatk.broadinstitute.org/hc/en-us/articles/360035891051-GenomicsDB)

Tip1: Create DB where you have spare storage, especially for growing data.
Tip2: If you are extending DB (multiple imports) make a local copy of your GenomicsDB (`<database_name>_<genomic_loc>_genomics_db` folders and `<database_name>.dbconf.json`) in case of the DB gets corrupted during an erroneous import
Tip3: Pipe logs to the file. It might be helpful to trace any ERRORs or WARNINGs

```
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file samples_test.txt \
    --prefix <database_name> \
    --ref <reference.fasta> \
    --vcf-dir .
    --threads 4
```

# Genotype and merge VCFs
After import information is genotyped from GenomicsDB and merged into a single VCF file. Only VCFs imported to GenomicsDB will be included in the file. Output `<prefix>.genotyped.vcf.gz`

```
python ~/software/fastq2matrix/scripts/merge_vcfs.py genotype \
    --prefix <database_name> \
    --ref <reference.fasta>
    --threads 4
```

# Filter genotyped VCF
File created with previous script is a main input (`--merged-file`) along with VCF with known, true snps (`​--bqsr-vcf`), bed file with core genome (`--include-region`) and annotation file (`--gff-file`). Last stage is annotation with `bcftools csq`. After each filter you get 'checkpoint' file produced (can be rocognized by suffix), that can be later used in any analysis.

Here, threshold like missingness and GT need to be adjusted to the overall quality of your data (less/more stringent).

__Example: `<prefix>.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz`__

Filtering stage | args | default | suffix | samples | snps |
--- | --- | --- | --- | --- | --- |
(0) Select core genome | --include-region | | | | x|
(1) Select SNPs | | | snps |  | x |
(2) Select indels | | | indel |  | x |
(3) Missing calls | --missing-sample-cutoff | 0.2 | miss | x | |
(4) GATK VQSR filtering | --VQSLOD | 0 | vqslod.filt |  | x |
(5) Mixed calls | --cutoff-mix-GT | 0.8 | GT |  | x |
(6) Bi-allelic calls |  | | bi |  | x |


```
python ~/software/fastq2matrix/scripts/filter_merged_vcf.py \
    --merged-file <genotyped_file.vcf.gz> \
    --prefix <prefix> \
    --ref <reference.fasta> \
    --bqsr-vcf <vcf_with_known_snps.vcf.gz> \
    --include-region <core_genome.bed> \
    --vqslod 0 \
    --missing-sample-cutoff 0.2 \
    --cutoff-mix-GT 0.8 \
    --gff-file <annotation.gff3>
```

# Convert to matrix
Any VCF file can be converted to matrix format including both nucleotide and binary version. Specify format of missing data with `--na` argument. For nucleotide version keep arg `--no-iupacgt`.

Output:
* `mat.bin` - binary matrix (0 - reference, 1 - alternative, 0.5 - mixed, specifies as in `--na` - missing)
* `noiupacgt.mat` - nucleotide matrix (specifies as in `--na` - missing)

```
python ~/software/fastq2matrix/scripts/vcf2matrix.py \
    --vcf <filtered_file.vcf.gz> \
    --no-iupacgt \
    --na NA \
    --threads 4
```

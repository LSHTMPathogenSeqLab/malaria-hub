options(scipen = 999)

# Load packages
require(data.table)
require(dplyr)
require(rehh)
library(optparse)
library(tidyr)
library(stringr)

source("~/software/malaria-hub/utils/helpers.R")

option_list = list(
  make_option(c("-d", "--workdir"), type = "character", default = '.',
              help = "Specify main directory",
              metavar = "character"),
  make_option(c("-b", "--matrix_binary"), type = "character", default = NULL,
              help = "Input filename of filtered binary matrix",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Full dir to metadata file",
              metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Input annotation file",
              metavar = "character"),
  make_option(c("-c", "--category"), type = "character", default = NULL,
              help = "Input filename of filtered binary matrix",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Label name in metadata for category column - [default %default]",
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = "fws",
              help = "Label name in metadata for fws column - [default %default]",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold [default %default]",
              metavar = "number"),
  make_option(c("--label_id"), type = "character", default = "sra_run",
              help = "Label name in metadata for id column - [default %default]",
              metavar = "character"),
  make_option(c("--maf"), type = "numeric", default = 0.01,
              help = "MAF threshold [default %default]",
              metavar = "number"),
  make_option(c("--na_char"), type = "character", default = "NA",
              help = "Specify NA characters",
              metavar = "character"),
  make_option("--forced_recode", type = "logical", default = FALSE,
              action = "store_true",
              help = "Recode missing to REF and mixed to ALT"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to remove ex. Pf3D7_API_v3,Pf_M76611",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection. Default matches Pf3D7_01_v3",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Regex pattern group",
              metavar = "numeric"),
  make_option(c("--threads"), type = "integer", default = 4,
              help = "Specify threads [default %default]",
              metavar = "number")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

setDTthreads(opt$threads)

## Storing binary matrix file
bin_mat_file <- opt$matrix_binary
## Storing metadata file name
metadata_file <- opt$metadata
## Annotation file
annotation_file <- opt$annotation
## Category
category <- opt$category
## Storing metadata_file field wth category groupings
label_category <- opt$label_category
## Storing metadata_file field with sample ids
label_id <- opt$label_id
## Storing metadata_file id with fws information
label_fws <- opt$label_fws
## Storing threshold for fws
threshold_fws <- opt$fws_th
## MAF threshold
th_maf <- opt$maf
# Working directory
workdir <- opt$workdir
# Missing calls character
na_char <- opt$na_char
# Recode data
forced_recode <- opt$forced_recode
# Remove chromosomes
rm_chr <- opt$remove_chr
# Pattern for chromosome detection
pattern <- opt$regex_chr
# Pattern group
groupid <- opt$regex_groupid

# Load annotation file
annotation <- read.table(annotation_file, sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = TRUE)

# Load metadata and check available categories
metadata <- read.csv(metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!all(c(label_id, label_category, label_fws) %in% colnames(metadata))) {
  stop("Wrongly specified metadata labels. Check again. Stopping...")
}

if (label_category %in% colnames(metadata)) {
   av_category <- metadata %>% select(!!sym(label_category)) %>% distinct() %>% pull()
    message("Available category options: ", paste0(av_category, collapse = ", "))
    cat("\n")
} else {
    stop("Wrongly specified category. Column name does not exist in metadata file.\n")
}

# Load subset of SNP binary matrix for selected category
if (!is.null(category)) {
  if (length(category) == 1 & category %in% av_category) {
    message(category, " found, Processing...\n")
    metadata <- metadata %>% filter(!!sym(label_category) == category)
    samples <- c("chr", "pos", "ref", (metadata %>% select(!!sym(label_id)) %>% pull() %>% as.vector()))
    snp <- fread(bin_mat_file, sep = "\t", select = samples, header = TRUE, data.table = FALSE)
  } else {        
    stop("Can not find ", category, " category. Exitting...\n")
  }
} else {
  stop("Category not specified. Exitting...\n")
}

# Trim category name
category_str <- as.character(gsub(" ", "_", category))

# Filter chromosome from matrix and annotation
if (!is.null(rm_chr)) {
  rm_chr <- strsplit(rm_chr, ",")[[1]]
  if (all(rm_chr %in% unique(snp$chr))) {
    snp <- snp %>% filter(!chr %in% rm_chr)
    annotation <- annotation %>% filter(!Chr %in% rm_chr)
  } else {
    stop("Wrong name for chromosomes to remove. Stopping...")
  }
} else {
  message("None chromosomes removed. Api and mito should be removed!")
}

# Transform chromosome from string to numeric
snp$chr <- as.numeric(stringr::str_match(snp$chr, pattern)[, groupid])
annotation$Chr <- as.numeric(stringr::str_match(annotation$Chr, pattern)[, groupid])

# Check if all samples match between binary matrix and metadata file
metadata <- metadata %>% filter(!!sym(label_id) %in% colnames(snp[, -c(1:3)]))
if (all(metadata[[label_id]] == colnames(snp[, -c(1:3)]))) {
  message("Matrix matches metadata. Proceeding...")
}

# Recode missing data
if (!is.na(na_char)) {
    snp[snp == na_char] <- NA
    snp[snp == "."] <- NA
}

# Separate SNP calls only matrix
snp_c <- as.data.frame(snp[, -(1:3)])
snp_c[] <- lapply(snp_c, as.numeric)

# Seperate SNP description chr,pos,ref
snp_d <- as.data.frame(snp[, (1:3)])

rm(snp)

maj3 <- snp_c
if (forced_recode) {
# Set NA to ref i.e. 0
  maj3[is.na(snp_c)] <- 0
  maj3[snp_c == 0.5] <- 1
} else {
  maj3[snp_c == 0.5] <- NA
}

rm(snp_c)

# Calculate MAF
maf_sti <- calculate_maf(maj3)
to_keep_sti <- which(maf_sti >= th_maf)

rm(maf_sti)

# Create MAP file with MAF filtered SNPs
snp_d <- snp_d[to_keep_sti, ]
snp_annot <- snp_d %>% left_join(annotation, by = c("chr" = "Chr", "pos" = "Pos", "ref" = "Ref")) %>%
              tidyr::unite("info", c(chr, pos), sep = "_", remove = FALSE)

map <- snp_annot %>% select(info, chr, pos, ref, Alt_1)
write.table(map, file.path(workdir, sprintf("snp.info.inp.%s", category_str)),
 quote = FALSE, row.names = FALSE, col.names = FALSE)

# Apply MAF filter on snps and Fws for samples
maj4 <- maj3[to_keep_sti, ]
maj4 <- maj4 %>% select(metadata %>%
    filter(!!sym(label_fws) >= threshold_fws) %>%
    select(!!sym(label_id)) %>%
    pull() %>%
    as.vector())

rm(maj3)

# Create hap file
hap <- t(maj4)
hap_c <- hap
hap_c[hap == 1] <- 2
hap_c[hap == 0] <- 1
hap_c[is.na(hap)] <- 0
colnames(hap_c) <- NULL
rownames(hap_c) <- NULL

rm(maj4, hap)

hap_c <- as.data.frame(hap_c, stringsAsFactors = FALSE)

i <- sapply(hap_c, is.character)
hap_c[i] <- lapply(hap_c[i], as.numeric)

# Create haplotypes per each chromosome with rehh::data2haplohh() and rehh::scan_hh()
u_chr <- unique(map$chr)

sink(file.path(workdir, paste0(category_str, ".log")), append = FALSE)
cat("## Create haplotypes: data2haplohh() & scan_hh() ## \n")
results_hh <- c()
for (uchr in u_chr) {
    hap_chr_s <- hap_c[, which(map$chr == uchr)]
    write.table(hap_chr_s, file.path(workdir, sprintf("hap_chr%d_%s", uchr, category_str)),
     sep = "\t", col.names = FALSE, quote = FALSE)
    hap_chr_pop <- data2haplohh(hap_file = file.path(workdir, sprintf("hap_chr%d_%s", uchr, category_str)),
                                map_file = file.path(workdir, sprintf("snp.info.inp.%s", category_str)),
                                recode.allele = FALSE,
                                chr.name = uchr,
                                min_perc_geno.hap = 80,
                                min_perc_geno.mrk = 70,
                                min_maf = 0)
    res_chr_s <- scan_hh(hap_chr_pop)
    results_hh <- rbind(results_hh, res_chr_s)
    cat("\n")
}
sink()

# Store scan_hh results per category
results_hh$category <- category_str
write.table(results_hh,
            file.path(workdir, sprintf("scanned_haplotypes_%s.tsv", category_str)),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

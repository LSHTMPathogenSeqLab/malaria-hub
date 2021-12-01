options(scipen = 999)
require(optparse)
require(scales)
require(data.table)
require(dplyr)
require(crayon)

source("~/software/malaria-hub/utils/helpers.R")

option_list = list(
  make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Specify main directory",
              metavar = "character"),
  make_option(c("-b", "--binary_matrix"), type = "character", default = NULL,
              help = "Input filename of filtered binary matrix",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Full dir to metadata file",
              metavar = "character"),
  make_option(c("-c", "--category"), type = "character", default = NULL,
              help = "Name of country/region",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Label name in metadata for category column",
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = "fws",
              help = "Label name in metadata for fws column",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold",
              metavar = "number"),
  make_option(c("--label_id"), type = "character", default = "sra_run",
              help = "Label name in metadata for id column",
              metavar = "character"),
  make_option(c("--maf"), type = "numeric", default = 0.01,
              help = "MAF threshold [default %default]",
              metavar = "number"),
  make_option(c("--na_char"), type = "character", default = "NA",
              help = "Specify NA characters",
              metavar = "character"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Specify threads number",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to remove ex. Pf3D7_API_v3,Pf_M76611",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection. Default matches Pf3D7_01_v3",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Regex pattern group",
              metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Working directory
workdir <- opt$workdir
# Binary matrix file name
bin_mat_file <- opt$binary_matrix
# Metadata file name
met_file <- opt$metadata
# Category
category <- opt$category
# Metadata field wth category groupings
label_category <- opt$label_category
# Metadata field with sample ids
label_id <- opt$label_id
# Metadata field with fws information
label_fws <- opt$label_fws
# Fws hreshold
threshold_fws <- opt$fws_th
# MAF threshold
th_maf <- opt$maf
# Threads number
threads <- opt$threads
# Missing calls character
na_char <- opt$na_char
# Remove chromosomes
rm_chr <- opt$remove_chr
# Pattern for chromosome detection
pattern <- opt$regex_chr
# Pattern group
groupid <- opt$regex_groupid

# Specify threads number
setDTthreads(threads)

# Load metadata
metadata <- read.csv(met_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (label_category %in% colnames(metadata)) {
    av_category <- metadata %>%
      select(!!sym(label_category)) %>%
      distinct() %>% pull()
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
        snp <- fread(bin_mat_file, sep = "\t", select = samples,
            header = TRUE, data.table = FALSE)
    } else {
        stop("Can not find " %+% bgBlue(category) %+% " category. Exitting...\n")
    }
} else {
    message("Category not specified. Processing for all...\n")
    snp <- fread(bin_mat_file, sep = "\t", header = TRUE, data.table = FALSE)
}

# Trim category name
category_str <- as.character(gsub(" ", "_", category))

# Filter chromosome from matrix
if (!is.null(rm_chr)) {
    rm_chr <- strsplit(rm_chr, ",")[[1]]
    if (any(rm_chr %in% unique(snp$chr))) {
        snp <- snp %>% filter(!chr %in% rm_chr)
    } else {
        stop("Wrong name for chromosomes to remove. Stopping...")
    }
} else {
    message("None chromosomes removed. Api and mito should be removed!")
}

# Transform chromosome from string to numeric
snp$chr <- as.numeric(stringr::str_match(snp$chr, pattern)[, groupid])

# Check if all samples match between binary matrix and metadata file
metadata <- metadata %>% filter(!!sym(label_id) %in% colnames(snp[, -c(1:3)]))
if (all(metadata[[label_id]] == colnames(snp[, -c(1:3)]))) {
  message("Matrix matches metadata. Proceeding...")
}

# Create input file for hmmIBD
# Recode missing data
if (!is.na(na_char)) {
    snp[snp == na_char] <- NA
    snp[snp == "."] <- NA
}

# Reformat matrix to fit hmmIBD format
snp_hmmibd <- snp
snp_hmmibd <- snp_hmmibd[, -3]
snp_hmmibd[snp_hmmibd == 0.5] <- 1
snp_hmmibd[is.na(snp_hmmibd)] <- -1
colnames(snp_hmmibd)[1:2] <- c("chrom", "pos")
snp_hmmibd_1 <- snp_hmmibd[, 1:2]
snp_hmmibd_2 <- snp_hmmibd[, -c(1:2)]
snp_hmmibd_2[] <- lapply(snp_hmmibd_2, as.numeric)

rm(snp, snp_hmmibd)

# Subset reformatted matrices
snp_hmmibd_02_1 <- snp_hmmibd_1
write.table(snp_hmmibd_02_1,
    file.path(workdir, "ibd_matrix_hap_leg.tsv"),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

snp_hmmibd_02_2 <- snp_hmmibd_2
write.table(snp_hmmibd_02_2,
    file.path(workdir, sprintf("ibd_matrix_hap_%s.tsv", category_str)),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Run hmmIBD for population
message(category)
  
# Calculate MAF
maf_sti <- calculate_maf(snp_hmmibd_02_2)
to_keep_sti <- which(maf_sti >= th_maf)

# Select samples of category with Fws >= threshold_fws
fws_samples <- metadata %>%
    filter(!!sym(label_fws) >= threshold_fws) %>%
    select(!!sym(label_id)) %>%
    pull() %>%
    as.vector()
  
# Select SNPs that passed MAF threshold
snp_hmmibd_02_2 <- snp_hmmibd_02_2[to_keep_sti, ]
snp_hmmibd_leg <- snp_hmmibd_02_1[to_keep_sti, ]

# Select samples that passed Fws threshold
snp_hmmibd_02_2 <- snp_hmmibd_02_2 %>% select(all_of(fws_samples))

snp_hmmibd_merged <- cbind(snp_hmmibd_leg, snp_hmmibd_02_2)

# Write country matrix to file
write.table(format(snp_hmmibd_merged, digits = 0),
    file.path(workdir, sprintf("hmmIBD_%s_maf%s.txt", category_str, as.character(th_maf))),
    sep = "\t", quote = FALSE, row.names = FALSE)
  
# Run hmmIBD
string_i <- file.path(workdir,
    sprintf("hmmIBD_%s_maf%s.txt", category_str, as.character(th_maf)))
string_o <- file.path(workdir,
    sprintf("hmmIBD_%s_maf%s_out", category_str, as.character(th_maf)))
output <- system(command = sprintf("~/software/hmmIBD/hmmIBD -i %s -o %s", string_i, string_o), intern = TRUE)
write.table(output, file.path(workdir,
    sprintf("hmmIBD_run_%s.log", category_str)), quote = FALSE, col.names = FALSE, row.names = FALSE)
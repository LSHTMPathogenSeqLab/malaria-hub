options(scipen = 999)
library(dplyr)
library(readr)
library(tidyr)
library(data.table)
library(ggplot2)
library(optparse)

source("~/software/malaria-hub/utils/helpers.R")

option_list = list(
  make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Specify main directory",
              metavar = "character"),
  make_option(c("--list_category"), type = "character", default = NULL,
              help = "Specify category list",
              metavar = "character"),
  make_option(c("-l", "--legend"), type = "character",
              default = "ibd_matrix_hap_leg.tsv",
              help = "Specify SNPs legend",
              metavar = "character"),
  make_option("--gene_product", type = "character", default = NULL,
              help = "Gene product file",
              metavar = "character"),
  make_option(c("-r", "--ref_index"), type = "character", default = NULL,
              help = "File name for reference index",
              metavar = "character"),
  make_option("--maf", type = "numeric", default = 0.01,
              help = "MAF",
              metavar = "numeric"),
  make_option(c("-w", "--window_size"), type = "numeric", default = 10000,
              help = "Specify window size",
              metavar = "numeric"),
  make_option("--quantile_cutoff", type = "numeric", default = 0.95,
              help = "Quantile cut-off for annotated IBD segemnts",
              metavar = "numeric"),
  make_option("--suffix", type = "character",
              default = format(Sys.time(), "%d_%m_%Y")),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Specify threads number",
              metavar = "numeric"),
  make_option("--remove_chr", type = "character", default = NULL,
              help = "Chromosomes to remove ex. Pf3D7_API_v3,Pf_M76611",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection. Default matches Pf3D7_01_v3",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Regex pattern group",
              metavar = "numeric"),
  make_option("--NSNP", type = "numeric", default = 0,
              help = "Minimal numbrt of SNPs in segment",
              metavar = "numeric"),
  make_option("--LSEGMENT", type = "numeric", default = 0,
              help = "Minimal length of segment",
              metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Workdir
workdir <- opt$workdir
# Country list
category_list <- opt$list_category
# Legend file
legend <- opt$legend
# Annotation file
gene_product_file <- opt$gene_product
# Ref index
ref_index <- opt$ref_index
# MAF threshold
th_maf <- opt$maf
# Window size
window_size <- opt$window_size
# Quantile cut-off
th_quantile <- opt$quantile_cutoff
# Suffix
suffix <- opt$suffix
# Threads
threads <- opt$threads
# Remove chromosomes
rm_chr <- opt$remove_chr
# Pattern for chromosome detection
pattern <- opt$regex_chr
# Pattern group
groupid <- opt$regex_groupid
# Number of SNPs in segment
min_n_snp <- opt$NSNP
# Lenght of segment
min_l_seg <- opt$LSEGMENT

# Specify threads number
setDTthreads(threads)

# Load category list
if (file.exists(file.path(workdir, category_list))) {
    category_list <- read.table(file.path(workdir, category_list),
    header = FALSE, stringsAsFactor = FALSE)$V1
} else {
    stop("Can not locate file with category list. Stopping...")
}
# Load legend
if (file.exists(file.path(workdir, legend))) {
  snp_hmmibd_02_1 <- read_tsv(file.path(workdir, legend)) %>%
      as.data.frame()
} else {
  stop("Can not locate legend file. Stopping...")
}
# Load reference index
if (file.exists(ref_index)) {
  fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
    rename(chr = V1, end_chr = V2) %>%
    select(chr, end_chr) %>%
    mutate(start_chr = 1)

  # Remove chromosomes
  if (!is.null(rm_chr)) {
    rm_chr_ <- strsplit(rm_chr, ",")[[1]]
    if (any(rm_chr_ %in% unique(fai$chr))) {
      rm_chr_ <- rm_chr_[rm_chr_ %in% unique(fai$chr)]
      fai <- fai %>% filter(!chr %in% rm_chr_)
    } else {
      stop("Wrong name for chromosomes to remove. Stopping...")
    }
  } else {
    message("None chromosomes removed. API and MITO should be removed!")
  }
  fai$chr <- as.numeric(stringr::str_match(fai$chr, pattern)[, groupid])
} else {
  stop("Can not locate reference index. Stopping...")
}
    
combined_ibd_r <- c()
combined_fraction_r <- c()

for (category_n in category_list) {
    message(category_n)
    # Define data for category
    fraction <- file.path(workdir, sprintf("hmmIBD_%s_maf%s_out.hmm_fract.txt", category_n, as.character(th_maf)))
    message(fraction)
    hmm_ibd <- file.path(workdir, sprintf("hmmIBD_%s_maf%s_out.hmm.txt", category_n, as.character(th_maf)))
    message(hmm_ibd)

    # Read data
    if (file.exists(fraction) & file.exists(hmm_ibd)) {
        ibd_frac <- read_tsv(fraction, col_types = cols()) %>% as.data.frame()
        ibd_data <- read_tsv(hmm_ibd, col_types = cols()) %>% as.data.frame()
        if (nrow(ibd_frac) != 0 & nrow(ibd_data) != 0) {
          # Parse IBD information
          ibd_data <- ibd_data %>%
              unite(id, sample1, sample2, sep = "_", remove = FALSE) %>%
              mutate(total = n_distinct(id),
                    length = end - start + 1)
          ibd_conf <- ibd_data %>%
              filter(different == 0 & Nsnp > min_n_snp & length > min_l_seg) %>%
              arrange(chr, start, end) %>%
              rename(pos_start = start, pos_end = end)

          # Parse everything
          chr_vec <- unique(snp_hmmibd_02_1[, 1])
          results <- c()
          total <- unique(ibd_conf$total)
          #TODO collapse to functions
          for (k in seq_along(chr_vec)) {              
              data_chr <- ibd_conf %>% filter(chr == chr_vec[k])
              length_chr <- fai %>% filter(chr == chr_vec[k]) %>%
                  select(end_chr) %>% pull()
              
              # Specify all windows for chromosome (1:length_chr)
              windows <- data.frame(chr = chr_vec[k],
                                    win_start = seq(1, length_chr, by = window_size),
                                    win_end = seqlast(window_size, length_chr, by = window_size))
              
              x <- data.table(start = as.numeric(as.character(data_chr$pos_start)),
                end = as.numeric(as.character(data_chr$pos_end)))
              y <- data.table(start = as.numeric(as.character(windows$win_start)),
                end = as.numeric(as.character(windows$win_end)))

              setkey(y, start, end)
              overlaps <- foverlaps(x, y, type = "any", which = TRUE)
              matched <- cbind(data_chr[overlaps$xid, ], as.data.frame(y[overlaps$yid])) %>% #TODO make nicer
                      mutate(rel_start = ifelse(pos_start > start, pos_start, start),
                            rel_end = ifelse(pos_end < end, pos_end, end),
                            rel_length = rel_end - rel_start + 1)

              matched_w <- windows %>% full_join(matched, by = c("chr", "win_start" = "start", "win_end" = "end"))
              matched_wm <- matched_w %>% group_by(win_start, win_end) %>%
                  mutate(m_sum = sum(rel_length, na.rm = T),
                         m_sum_av = ifelse(!is.na(m_sum), m_sum / window_size, 0),
                         fraction = m_sum_av / total) %>%
                  select(chr, win_start, win_end, m_sum, m_sum_av, fraction) %>%
                  distinct()
              results <- rbind(results, matched_wm)
          }
          results_it <- results %>% as.data.frame() %>%
              mutate(category = category_n)
          combined_ibd_r <- rbind(combined_ibd_r, results_it)

          # Parse fraction information
          frac_it <- ibd_frac %>% unite(id, sample1, sample2, sep = "_") %>%
                                  mutate(category = category_n) %>%
                                  select(c(id, category, fract_sites_IBD))
          combined_fraction_r <- rbind(combined_fraction_r, frac_it)
      } else {
        message("hmmIBD results empty. Skipping...")
      }
    } else {
        message("hmmIBD results not found. Check path once again. Skipping...")
    }
}

message("Saving...")
if (length(combined_ibd_r) != 0 & length(combined_fraction_r) != 0) {
  #colnames(combined_ibd_r) <- c("chr", "start", "end", "fraction", "category")
  write.table(combined_ibd_r,
    file.path(workdir, sprintf("%s_hmmIBD_ibd.tsv", suffix)), #TODO save window, Nsnp, Length
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  colnames(combined_fraction_r) <- c("id", "category", "fraction")
  write.table(combined_fraction_r,
    file.path(workdir, sprintf("%s_hmmIBD_fraction.tsv", suffix)),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
} else {
  stop("None results to save. Stopping...")
}

# Annotate results
message("Annotating...")
if (length(combined_ibd_r) != 0 & length(combined_fraction_r) != 0) {
  # Load annotation
  if (!is.null(gene_product_file)) {
    if (file.exists(gene_product_file)) {
    annotation <- readr::read_tsv(gene_product_file, col_types = cols())

    # Remove chromosomes
    if (!is.null(rm_chr)) {
      rm_chr_ <- strsplit(rm_chr, ",")[[1]]
      if (any(rm_chr_ %in% unique(annotation$chr))) {
        rm_chr_ <- rm_chr_[rm_chr_ %in% unique(annotation$chr)]
        annotation <- annotation %>% filter(!chr %in% rm_chr_)
      } else {
        stop("Wrong name for chromosomes to remove. Stopping...")
      }
    } else {
      message("None chromosomes removed. API and MITO should be removed!")
    }

    # Transform chromosome from string to numeric
    annotation$chr <- as.numeric(stringr::str_match(annotation$chr, pattern)[, groupid])
    ibd_regions <- combined_ibd_r %>% rename("start" = "win_start", "end" = "win_end") %>% select(c(chr, start, end)) %>% distinct()

    # Find overlap
    res_annot <- annotate_candidate_regions(ibd_regions, annotation)

    # Calculate quantiles
    quantile <- combined_ibd_r %>% group_by(category) %>%
      mutate(qfrac = quantile(fraction, th_quantile, na.rm = T)) %>%
      filter(fraction >= qfrac) %>% ungroup() %>% distinct()

    quantile_annot <- quantile %>%
      rename("start" = "win_start", "end" = "win_end") %>%
      inner_join(res_annot) %>% ungroup()
    write.table(quantile_annot, file.path(workdir, sprintf("%s_hmmIBD_ibd_annotated_q%s.tsv", suffix, as.character(th_quantile))),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    # Wide format
    quantile_regs <- quantile %>% select(chr, win_start, win_end) %>% distinct()
    quantile_annot_wide <- combined_ibd_r %>% right_join(quantile_regs) %>%
      rename("start" = "win_start", "end" = "win_end") %>%
      left_join(res_annot) %>%
      select(-c(pos_start, pos_end)) %>%
      group_by(chr, start, end, category, fraction) %>%
      dplyr::summarise(genes = paste0(gene_name, "(", X6, ")", collapse = "; "),
                       products = paste0(gene_product, collapse = "; ")) %>%
      mutate(genes = gsub("\\(NA\\)", "", genes)) %>%
      mutate(fraction = round(fraction, 3)) %>%
      tidyr::pivot_wider(names_from = "category", values_from = "fraction") %>% ungroup()

      write.table(quantile_annot_wide, file.path(workdir, sprintf("%s_hmmIBD_ibd_annotated_q%s_wide.tsv", suffix, as.character(th_quantile))),
          sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
   } else {
          stop("Annotation file not found. Skiping....\n")
    }
  } else {
    stop("Specify annotation file. Stopping...")
  }
} else {
  stop("None results to annotate. Stopping...")
}
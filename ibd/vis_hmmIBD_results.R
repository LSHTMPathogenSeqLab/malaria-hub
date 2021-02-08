library(ggplot2)
library(readr)
library(dplyr)

source("~/software/malaria-hub/utils/helpers.R")

workdir <- "~/hmmIBD_tests"
metadata_file <- "~/hmmIBD_tests/pf_metadata_collapsed_w_region_np_June_2020.tsv"
country_label <- "country"
region_label <- "region"
suffix <- "04_02_2021"
ref_index <- "~/hmmIBD_tests/Pfalciparum.genome.fasta.fai"
rm_chr <- c("Pf_M76611", "Pf3D7_API_v3")
category_order <- c("South_America", "South_Central_Africa")

pattern <- "(.*?)_(.+)_(.*)"
groupid <- 3

# Load IBD
combined_ibd_r <- read_tsv(file.path(workdir, sprintf("%s_hmmIBD_ibd_results_combined.tsv", suffix)), col_types = cols())
# Load IBD fractions
fraction_ibd_r <- read_tsv(file.path(workdir, sprintf("%s_hmmIBD_fraction_results_combined.tsv", suffix)), col_types = cols())
# Load metadata
metadata <- read_tsv(metadata_file, col_types = cols()) %>%
  select(c(country_label, region_label))
# Reference index
fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  rename(chr = V1, end_chr = V2) %>%
  select(chr, end_chr) %>%
  mutate(start_chr = 1) %>%
  filter(!chr %in% rm_chr)

fai$chr <- as.numeric(stringr::str_match(fai$chr, pattern)[, groupid])
chrom_ends <- c(0, fai$end_chr[-max(NROW(fai$end_chr))])
transpose_chr <- data.frame(chr = fai$chr, tr_chr = chrom_ends) %>%
                     mutate(ind = seq(1, nrow(.)))

# Combine results wih region
combined_ibd_r <- combined_ibd_r %>% left_join(metadata, by = c("category" = country_label))
fraction_ibd_r <- fraction_ibd_r %>% left_join(metadata, by = c("category" = country_label))

# Boxplot fractions
g <- ggplot(data = fraction_ibd_r, aes(x = category, y = fraction, fill = !!sym(region_label))) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 9,
               size = 1, color = "yellow") +
  theme_classic() +
  guides(fill = guide_legend(title = "Region:", nrow = 1)) +
  theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(size = 12, color = "black",
                                     angle = 90, vjust = -0.5),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 15, color = "black", hjust = 0.5),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          strip.background = element_blank(),
          legend.position = "bottom") +
  labs(x = "Country", y = "Pairwise fraction IBD")

# Transform to genetic distribution
ibd_frac_tr <- combined_ibd_r %>% group_by(chr) %>%
    mutate(trans = get_chrom_transposition(transpose_chr, chr)) %>%
    mutate(pos_bp_ed = as.numeric(as.numeric(start) + trans),
            Fraction = as.numeric(fraction)) %>%
    ungroup()

# Establish order in plot
# Arrange according to order
ibd_frac_tr <- ibd_frac_tr %>% mutate(region = factor(!!sym(region_label), levels = category_order)) %>% arrange(region)
ibd_frac_tr <- ibd_frac_tr %>% mutate(country = factor(category, levels = unique(category)))

# IBD pairwise fraction in 10kb windows
p <- ggplot(data = ibd_frac_tr) +
    geom_line(aes(x = pos_bp_ed, y = fraction, color = region)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 1), labels = c("0.0", "1.0")) +
    facet_grid(category ~ ., space = "free_x") +
    labs(x = "Chromsome", y = "IBD Fraction") +
    guides(color = guide_legend(title = "Region:", nrow = 1, byrow = FALSE)) +
    theme_classic() +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(size = 12, color = "black", angle = 0, vjust = -0.5),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 15, color = "black", hjust = 0.5),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          strip.background = element_blank(),
          legend.position = "bottom") +
    geom_vline(data = ibd_frac_tr, aes(xintercept = trans), color = "black",
               alpha = 0.5, linetype = "longdash", size = 0.2) +
    scale_x_continuous(breaks = unique(ibd_frac_tr$trans),
                       labels = unique(ibd_frac_tr$chr))

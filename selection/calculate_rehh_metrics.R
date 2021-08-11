library(rehh)
library(dplyr)
library(data.table)
library(optparse)
library(ggplot2)
library(ggrepel)
library(stringr)
library(showtext)
showtext_auto()


source("~/software/malaria-hub/utils/helpers.R")

option_list = list(
    make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Specify main directory", metavar = "character"),
    make_option(c("-p", "--prefix"), type = "character", default = "scanned_haplotypes",
              help = "Prefix", metavar = "character"),
    make_option(c("--list_category"), type = "character", default = NULL,
              help = "Category list", metavar = "character"),
    make_option(c("--annotation"), type = "character", default = NULL,
              help = "Annotation", metavar = "character"),
    make_option(c("--gene_product"), type = "character", default = NULL,
              help = "Gene product", metavar = "character"),
    make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to remove ex. Pf3D7_API_v3,Pf_M76611",
              metavar = "character"),
    make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
                help = "Regex pattern for chromosome detection. Default matches Pf3D7_01_v3",
                metavar = "character"),
    make_option("--regex_groupid", type = "numeric", default = 3,
                help = "Regex pattern group",
                metavar = "numeric"),
    make_option(c("--ihs_th"), type = "integer", default = 4,
              help = "iHS p-value threshold",
              metavar = "number"),
    make_option(c("--rsb_th"), type = "integer", default = 5,
              help = "Rsb p-value threshold",
              metavar = "number"),
    make_option(c("--xpehh_th"), type = "integer", default = 5,
              help = "XP-EHH p-value threshold",
              metavar = "number"),
    make_option(c("--threads"), type = "integer", default = 4,
              help = "Specify threads [default %default]",
              metavar = "number")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# thresholds for iHS, Rsb, XP-EHH - p-value
ihs_th <- opt$ihs_th
rsb_th <- opt$rsb_th
xpehh_th <- opt$xpehh_th

# Set threads for DT
setDTthreads(opt$threads)

# workdir
workdir <- opt$workdir
# prefix
prefix <-  opt$prefix
# list_category
category_list <- opt$list_category
# annotation
annotation_file <- opt$annotation
# gene_product
gene_product_file <- opt$gene_product
# Chromosomes to remove
rm_chr <- opt$remove_chr
# Pattern for chromosome detection
pattern <- opt$regex_chr
# Pattern group
groupid <- opt$regex_groupid

# Y axis labels
ihs_expr <- expression("-" * log[10] * "[1" ~ "-" ~ "2" ~ "|" ~ Phi[scriptstyle(italic(iHS))] ~ "-" ~ 0.5 * "|]")
rsb_expr <- expression(-log[10] ~ "(" * italic(p) * "-value)")
xpehh_expr <- expression(-log[10] ~ "(" * italic(p) * "-value)")

# Load categories file
categories <- read.table(category_list, sep = "\n")$V1 %>% as.vector()

# Load annotation file Chr, Pos, Ref, Alt_1, Gene_name_1
annotation <- read.table(annotation_file, sep = "\t", fill = TRUE,
                         header = TRUE, stringsAsFactors = FALSE)
annotation <- annotation %>% select(c(chr, pos, ref, alt, Gene_name)) %>% distinct()

# Load gene/product file
gff_table <- read.csv(gene_product_file, sep = "\t",
                      header = TRUE, stringsAsFactors = FALSE) %>%
    as.data.frame()

# Filter chromosome from annotation and gene product table
if (!is.null(rm_chr)) {
  rm_chr <- strsplit(rm_chr, ",")[[1]]
  if (any(rm_chr %in% unique(annotation$chr))) {
    annotation <- annotation %>% filter(!chr %in% rm_chr)
    gff_table <- gff_table %>% filter(!chr %in% rm_chr)
  } else {
    stop("Wrong name for chromosomes to remove.")
  }
} else {
  message("None chromosomes removed")
}

# Transform chromosome names to numeric
annotation$chr <- as.numeric(stringr::str_match(annotation$chr, pattern)[, groupid])
gff_table$chr <- as.numeric(stringr::str_match(gff_table$chr, pattern)[, groupid])


high_ihs_all <- c()
cr_ihs_all <- c()

high_rsb_all <- c()
cr_rsb_all <- c()

high_xpehh_all <- c()
cr_xpehh_all <- c()

# Read IHH, IES, INES metrics and calculate iHS, Rsb, XPEHH metric
# * Plot manhattan plot
# * Filter sites with high significance
# * Detect candidate regions
for (category in categories) {
  # Per category
  message(category, "\n")
  pdf(file.path(workdir, sprintf("plots_%s.pdf", category)))
  sink(file.path(workdir, paste0(category, "_metrics.log")))

  filec <- file.path(workdir, paste0(prefix, "_", category, ".tsv"))
  if (file.exists(filec)) {
    #### IHS #####
    cat("\n## iHS ##\n")
    ihh <- fread(filec, sep = "\t", header = TRUE, data.table = FALSE)
    ihs <- rehh::ihh2ihs(ihh, min_maf = 0.0, freqbin = 0.05)

    if (nrow(ihs$ihs) > 1) {
      ihsA <- ihs$ihs %>% left_join(annotation, by = c("CHR" = "chr", "POSITION" = "pos"))
      gg_data <- gg_to_plot <- modify_df_ggplot(ihsA, th = ihs_th)

      gg <- generate_manhattan_ggplot(gg_data$df_vis, gg_data$df_axis,
                                th = ihs_th,
                                name = bquote(italic("iHS") ~ .(gsub("_", " " , category))),
                                yname = ihs_expr,
                                hcolor = "lightblue")
      saveRDS(gg, file.path(workdir, paste0("plot_", category, ".iHS.rds")))

      # Annotation
      high_ihs <- ihsA %>% filter(LOGPVALUE >= ihs_th)
      if (nrow(high_ihs) > 1) {
        high_ihs$category_name <- category
        high_ihs_all <- rbind(high_ihs_all, high_ihs)
      }
      # Candidate regions
      cr_ihs <- rehh::calc_candidate_regions(ihs$ihs,
                                             threshold = ihs_th,
                                             pval = TRUE,
                                             window_size = 2E4,
                                             overlap = 1E4,
                                             min_n_extr_mrk = 2)
      if (nrow(cr_ihs) > 1) {
        cr_ihs$category_name <- category
        cr_ihs_all <- rbind(cr_ihs_all, cr_ihs)
      }
    }

    # Pairwise comparisons
    if (length(categories >= 2)) {
      other_categories <- categories[-which(categories == category)]
      for (contr_category in other_categories) {
        cat(paste0("\n", contr_category, "\n"))

        fileoc <- file.path(workdir, paste0(prefix, "_", contr_category, ".tsv"))
        if (file.exists(fileoc)) {
          ihh_oc <- fread(fileoc, sep = "\t", header = TRUE, data.table = FALSE)
          #### Rsb ####
          cat("\n## Rsb ##\n")
          rsb <- rehh::ines2rsb(ihh, ihh_oc)

          if (nrow(rsb) > 1) {
            rsbA <- rsb %>% left_join(annotation, by = c("CHR" = "chr", "POSITION" = "pos"))
            gg_data <- gg_to_plot <- modify_df_ggplot(rsbA, th = rsb_th)

            gg <- generate_manhattan_ggplot(gg_data$df_vis, gg_data$df_axis,
                                      th = rsb_th,
                                      name = bquote(italic("Rsb") ~ .(gsub("_", " " , category)) ~ "vs." ~ .(gsub("_", " ", contr_category))),
                                      yname = rsb_expr,
                                      hcolor = "red")
            saveRDS(gg, file.path(workdir, paste0("plot_", category, ".", contr_category, ".Rsb.rds")))


            # High significance
            high_rsb <- rsbA %>% filter(LOGPVALUE >= rsb_th)
            if (nrow(high_rsb) > 1) {
              high_rsb$category_name <- paste0(sort(c(category, contr_category)), collapse = "|")
              high_rsb_all <- rbind(high_rsb_all, high_rsb)
            }

            # Candidate regions
            cr_rsb <- rehh::calc_candidate_regions(rsb,
                                                   threshold = rsb_th,
                                                   pval = TRUE,
                                                   window_size = 2E4,
                                                   overlap = 1E4,
                                                   min_n_extr_mrk = 2)
            if (nrow(cr_rsb) > 1) {
              cr_rsb$category_name <- paste0(sort(c(category, contr_category)), collapse = "|")
              cr_rsb_all <- rbind(cr_rsb_all, cr_rsb)
            }
          }

          #### XPEHH #####
          cat("\n## XPEHH ##\n")
          xpehh <- rehh::ies2xpehh(ihh, ihh_oc)

          if (nrow(xpehh) > 1) {
            xpehhA <- xpehh %>% left_join(annotation, by = c("CHR" = "chr", "POSITION" = "pos"))
            gg_data <- gg_to_plot <- modify_df_ggplot(xpehhA, th = xpehh_th)

            gg <- generate_manhattan_ggplot(gg_data$df_vis, gg_data$df_axis,
                            th = xpehh_th,
                            name = bquote(italic("XP-EHH") ~ .(gsub("_", " " , category)) ~ "vs." ~ .(gsub("_", " ", contr_category))),
                            yname = xpehh_expr,
                            hcolor = "purple")
            saveRDS(gg, file.path(workdir, paste0("plot_", category, ".", contr_category, ".XPEHH.rds")))

            high_xpehh <- xpehhA %>% filter(LOGPVALUE >= xpehh_th)
            if (nrow(high_xpehh) > 1) {
              high_xpehh$category_name <- paste0(sort(c(category, contr_category)), collapse = "|")
              high_xpehh_all <- rbind(high_xpehh_all, high_xpehh)
            }

            cr_xpehh <- rehh::calc_candidate_regions(xpehh,
                                                     threshold = xpehh_th ,
                                                     pval = TRUE,
                                                     window_size = 2E4,
                                                     overlap = 1E4,
                                                     min_n_extr_mrk = 2)
            if (nrow(cr_xpehh) > 1) {
              cr_xpehh$category_name <- paste0(sort(c(category, contr_category)), collapse = "|")
              cr_xpehh_all <- rbind(cr_xpehh_all, cr_xpehh)
            }
          }
        }
      }
    } else {
      message('Only iHS results calculated. Not enough populations for Rsb, XPEHH.')
    }
  }
  sink()
  dev.off()
}

# Save iHs, Rsb, XP-EHH results for all categories
# iHS
if (length(high_ihs_all) != 0) {
   write.table(high_ihs_all, file.path(workdir, "high_ihs_all_categories.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
 
if (length(cr_ihs_all) != 0) {
  cr_ihs_all <- cr_ihs_all %>%
    rename("chr" = "CHR", "start" = "START", "end" = "END") %>% distinct()
  cr_ihs_ann <- annotate_candidate_regions(cr_ihs_all, gff_table) %>% select(-c(pos_start, pos_end)) %>%
      group_by(chr, start, end, category_name, N_MRK, MEAN_MRK, MAX_MRK, N_EXTR_MRK, PERC_EXTR_MRK, MEAN_EXTR_MRK) %>%
      dplyr::summarise(genes = paste0(gene_id, "(", gene_name, ")", collapse = "; "),
                       products = paste0(gene_product, collapse = "; ")) %>%
      mutate(genes = gsub("\\(\\)", "", genes),
             products = gsub("\\t", "", products)) %>% ungroup() %>%
      mutate_if(is.numeric, round, 3) %>% distinct()
  write.table(cr_ihs_ann, file.path(workdir, "cr_ihs_all_categories_annot.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

# Rsb
if (length(high_rsb_all) != 0) {
  write.table(high_rsb_all, file.path(workdir, "high_rsb_all_categories.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
if (length(cr_rsb_all) != 0) {
  cr_rsb_all <- cr_rsb_all %>%
    rename("chr" = "CHR", "start" = "START", "end" = "END") %>% distinct()
    write.table(cr_rsb_all, file.path(workdir, "test.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  cr_rsb_ann <- annotate_candidate_regions(cr_rsb_all, gff_table) %>% select(-c(pos_start, pos_end)) %>%
      group_by(chr, start, end, category_name, N_MRK, MEAN_MRK, MAX_MRK, N_EXTR_MRK, PERC_EXTR_MRK, MEAN_EXTR_MRK) %>%
      dplyr::summarise(genes = paste0(gene_id, "(", gene_name, ")", collapse = "; "),
                       products = paste0(gene_product, collapse = "; ")) %>%
      mutate(genes = gsub("\\(\\)", "", genes),
             products = gsub("\\t", "", products)) %>% ungroup() %>%
      mutate_if(is.numeric, round, 3) %>% distinct()
  write.table(cr_rsb_ann, file.path(workdir, "cr_rsb_all_categories_annot.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

# XPEHH
if (length(high_xpehh_all) != 0) {
  write.table(high_xpehh_all, file.path(workdir, "high_xpehh_all_categories.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

if (length(cr_xpehh_all) != 0) {
  cr_xpehh_all <- cr_xpehh_all %>%
    rename("chr" = "CHR", "start" = "START", "end" = "END") %>% distinct()
  cr_xpehh_ann <- annotate_candidate_regions(cr_xpehh_all, gff_table) %>%
    select(-c(pos_start, pos_end)) %>%
    group_by(chr, start, end, category_name, N_MRK, MEAN_MRK, MAX_MRK, N_EXTR_MRK, PERC_EXTR_MRK, MEAN_EXTR_MRK) %>%
    dplyr::summarise(genes = paste0(gene_id, "(", gene_name, ")", collapse = "; "),
                     products = paste0(gene_product, collapse = "; ")) %>%
    mutate(genes = gsub("\\(\\)", "", genes),
           products = gsub("\\t", "", products)) %>% ungroup() %>%
    mutate_if(is.numeric, round, 3) %>% distinct()
  write.table(cr_xpehh_ann, file.path(workdir, "cr_xpehh_all_categories_annot.tsv"),
  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
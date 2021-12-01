## Script to hold helper functions for selection scripts
## prepare_input_rehh_per_category.R
## calculate_rehh_metrics.R

# Count reference snp
ref <- function(x) {sum(1 * (as.numeric(x) == 0), na.rm = TRUE)}
# Count alternative snp
alt <- function(x) {sum(1 * (as.numeric(x) == 1), na.rm = TRUE)}

# Calculate MAF frequency for bi-allelic snps
calculate_maf <- function(x) {
  cref <- apply(x, 1, ref)
  calt <- apply(x, 1, alt)
  af <- calt / (calt + cref)
  maf <- ifelse(af > 0.5, 1 - af, af)
}

# Modify rehh table to fit qqman:manhattan() input requirements
modify_df_qqman <- function(df) {
  dfsel <- df %>% mutate(CHR = as.numeric(CHR),
                        BP = as.numeric(POSITION),
                        P = as.numeric(LOGPVALUE))
}


# Generate simple manhattan plot with qqman package
manhattan_plot <- function(df, title, yname, colors = NULL) {
    if (is.null(colors)) {
      colors <- c("black", "grey")
    }
    ymax <- ceiling(max(df$P, na.rm = TRUE))
    manhattan(df,
              logp = FALSE,
              main = title,
              genomewideline = 5,
              suggestiveline = 20,
              ylab = yname,
              col = colors,
              ylim = c(0, ymax))
}

# Modify rehh table for ggplot manhattan visualization
modify_df_ggplot <- function(df, th=4) {
  # Add highlight for snps with P > th
  # Add label for genes with at least 2 significant snps
  df <- df %>% 
    select(CHR, POSITION, LOGPVALUE, alt, Gene_name) %>%
    mutate(CHR = as.numeric(CHR),
          POSITION = as.numeric(POSITION),
          LOGPVALUE = as.numeric(LOGPVALUE)) %>%
    filter(!is.na(LOGPVALUE)) %>%
    group_by(Gene_name) %>%
    mutate(pc = sum(LOGPVALUE >= th),
           pcmax = max(LOGPVALUE, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(is_annotate = ifelse(Gene_name != "" & pc >= 2 & LOGPVALUE == pcmax, "yes", "no"),
           is_highlight = ifelse(LOGPVALUE >= th, "yes", "no"))
          
  # Compute chromosome sizes
  mod_df <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POSITION), .groups = "drop") %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df, ., by = c("CHR" = "CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POSITION) %>%
  mutate(BPcum = POSITION + tot) %>%
  ungroup()

  axis_df <- mod_df %>%
   group_by(CHR) %>%
   summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = 'drop')

  list("df_vis" = mod_df, "df_axis" = axis_df)
}

# Generate manhattan plot for iHS, Rsb, XPEHH results
generate_manhattan_ggplot <- function(df, axis, th, name, yname, hcolor = "orange") {

  p <- ggplot(data = df, aes(x = BPcum, y = LOGPVALUE)) +
      geom_point(aes(color = as.factor(CHR)), alpha = 1, size = 1.3) +
      scale_color_manual(values = rep(c("black", "grey"), 22)) +
      scale_x_continuous(label = axis$CHR, breaks = axis$center) +
      scale_y_continuous(breaks = seq(0, ceiling(max(df$LOGPVALUE, na.rm=TRUE)) + 1, 2),
                         limits=c(first(seq(0, ceiling(max(df$LOGPVALUE, na.rm=TRUE)) + 1, 2)),
                                  last(seq(0, ceiling(max(df$LOGPVALUE, na.rm=TRUE)) + 1, 2)))) +
      geom_point(data = subset(df, is_highlight == "yes"),
                 color = hcolor, size = 1.5) +
      geom_hline(yintercept = th, color = "red", alpha = 1) +
      geom_label_repel(data = (df %>% filter(is_annotate == "yes") %>%
                       group_by(Gene_name) %>% top_n(1, row_number())),
                       aes(label = Gene_name),
                       size = 3,
                       box.padding = 0.25,
                       label.padding = 0.35,
                       position = "identity") +
      labs(title = name,
           x = "Chromosome",
           y = yname) +
      theme_classic() +
      theme(
          legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color = "black"),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_text(size=12, color = "black"),
          plot.title = element_text(size=15, color="black", hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
      )
    print(p)
}

# Annotate genomic regions
annotate_candidate_regions <- function(cr_res, annot) {
  # Creating overlap to identify genes in candidate regions
  x <- data.table(chr = as.numeric(annot$chr),
                  start = as.numeric(as.character(annot$pos_start)),
                  end = as.numeric(as.character(annot$pos_end)))
  y <- data.table(chr = as.numeric(as.character(cr_res$chr)),
                  start = as.numeric(as.character(cr_res$start)),
                  end = as.numeric(as.character(cr_res$end)))
  
  data.table::setkey(y, chr, start, end)
  overlaps <- data.table::foverlaps(x, y, type = "any", which = TRUE)

  cr_res <- cr_res %>% 
      arrange(chr, start, end) %>%
      mutate(idR = row_number())
  annot <- annot %>% dplyr::mutate(idA = row_number())
  df_overlaps <- as.data.frame(overlaps)

  # Merging
  res_annot <- cr_res %>%
    left_join(df_overlaps, by = c("idR" = "yid")) %>%
    left_join(annot, by = c("xid" = "idA", "chr"))
  res_annot <- res_annot %>% dplyr::select(-c("idR", "xid"))

  return(res_annot)
}

# Transpose genomic positions
get_chrom_transposition <- function(chrom_map, str_chr) {
  if (length(str_chr) > 1) {
    chromosome <- as.character(unique(str_chr))
  } else {
    chromosome <- as.character(str_chr)
  }
  if (length(chromosome) == 1) {
    chrom_map <- chrom_map %>% mutate(chr = as.character(chr))
    i <- chrom_map[which(chrom_map$chr == chromosome), ]$ind
    result <- chrom_map %>% filter(ind <= i) %>% select(tr_chr) %>% sum()
    if (length(str_chr) > 1) {
      result <- rep(result, length(str_chr))
    }
  } else {
    stop("Cannot use function.")
  }
  return(result)
}

# Seq vector with last element
# https://stackoverflow.com/questions/28419281/missing-last-sequence-in-seq-in-r
seqlast <- function (from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if (tail(vec, 1) != to) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(forcats)
library(optparse)
library(unikn)
#library(showtext)
library(countrycode)
library(stringi)
#showtext_auto()

options(scipen = 999, digits = 3)

option_list = list(
 make_option(c("-d", "--workdir"), type = "character", default = getwd(),
             help = "Specify main directory", metavar = "character"),
 make_option(c("--prefix"), type = "character", default = NULL,
             help = "Prefix of admixture output files <prefix>.Q",
             metavar = "character"),
 make_option(c("--filename"), type = "character",
             default = paste("barplot", format(Sys.Date(), "%d_%m_%y"),
                    stringi::stri_rand_strings(1, 10), sep = "_"),
             help = "Filename for plots", metavar = "character"),
 make_option(c("-k", "--kval"), type = "character", default = NULL,
             help = "Comma-separated K values", metavar = "character"),
 make_option(c("-m", "--metadata"), type = "character", default = NULL,
             help = "Metadata filepath with sample, region, country
                   (optionally) site information", metavar = "character"),
 make_option(c("--region_order"), type = "character", default = NULL,
             help = "Region order", metavar = "character"),
 make_option(c("-s", "--select_country"), type = "character", default = NULL,
             help = "Filter per country", metavar = "character"),
 make_option(c("-f", "--filter_N"), type = "numeric", default = 20,
             help = "Select population with N > [default: %default]",
             metavar = "integer"),
 make_option(c("--label_region"), type = "character", default = "region",
            help = "Region label name in metadata file", metavar = "character"),
 make_option(c("--label_country"), type = "character", default = "country",
            help = "Country label name in metadata file",
            metavar = "character"),
 make_option(c("--label_site"), type = "character", default = NULL,
            help = "Region label name in metadata file", metavar = "character"),
 make_option(c("--label_id"), type = "character", default = "sra_run",
            help = "Category label name in metadata file",
            metavar = "character"),
 make_option(c("--country_code"), type = "character", default = FALSE,
            help = "Country codes", metavar = "logical"),
 make_option(c("--axisx_angle"), type = "numeric", default = 90,
            help = "Rotation angle for X axis labels",
            metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

workdir <- opt$workdir
prefix <- opt$prefix
filename <- opt$filename
Kval <- stringr::str_split(opt$kval, ",")[[1]] %>% as.numeric()
metadata_file <- opt$metadata
order <- opt$region_order
selected <- opt$select_country
n_filter <- as.numeric(opt$filter_N)
label_id <- opt$label_id
label_country <- opt$label_country
label_region <- opt$label_region
label_site <- opt$label_site
met_cols <- c(label_id, label_region, label_country, label_site)

pals <- unikn::usecol(c(Seeblau, Petrol, Bordeaux, Pinky, Seegruen, "gold"),
        n = max(Kval))

## ggplot fun to plot admix-like plots
plot_gg_admix <- function(df, label_id, label_cat, type, palette, l_angle) {
  p <- ggplot(data = df, aes(x = !!sym(label_id), y = Q, fill = K)) +
    geom_col(width = 1) +
    facet_grid(as.formula(paste0("def", "~", label_cat)),
               switch = "x", scales = "free", space = "free") +
    theme_minimal(base_size = 20) +
    labs(x = "Individuals",
         title = sprintf("Admixture - %s", type),
         y = "Ancestry",
         fill = "K") +
    scale_y_continuous(expand = c(0, 0.05), breaks = c(0, 0.5, 1.0)) +
    scale_fill_manual(values = palette, guide = guide_legend(nrow = 1)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    theme(
      legend.position = "top",
      strip.text.x = element_text(angle = l_angle,
                                  hjust = 0.5, vjust = 0.5,
                                  size = 20),
      panel.spacing.x = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
print(p)
}

# Read sample list and metadata
samples <- read.table(file.path(workdir, paste0(prefix, ".nosex")),
                      header = FALSE, stringsAsFactors = FALSE)$V2
metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE)

# Check if column labels in metadata
if (all(met_cols %in% colnames(metadata))) {
  metadata <- metadata %>%
    select(all_of(met_cols))
} else {
  stop("Labels not existing in metadata. Support all information")
}

# Select only population size > filter_N
metadata <- metadata %>%
  filter(!!sym(label_country) %in% (metadata %>%
  select(!!sym(label_country)) %>%
  plyr::count() %>%
  filter(freq >= n_filter) %>%
  select(!!sym(label_country)) %>%
  pull()))

# Select country
if (!is.null(selected)) {
  selected <- strsplit(selected, ",")[[1]]
  if (all(selected %in% metadata[[label_country]])) {
    metadata <- metadata %>% filter(!!sym(label_country) %in% selected)
  }
}

# Adjust order
if (!is.null(order)) {
  if (all(order %in% metadata[[label_region]])) {
    order <- gsub("_", " ", strsplit(order, ",")[[1]])
  } else {
    order <- gsub("_", "", sort(unique(metadata[[label_region]])))
  }
} else {
  order <- gsub("_", " ", sort(unique(metadata[[label_region]])))
}

admix_all <- c()
for (kk in Kval) {
  qval <- read.table(file.path(workdir, paste0(prefix, ".", kk, ".", "Q")),
                     header = FALSE, stringsAsFactors = FALSE)
  colnames(qval) <- gsub("V", "K", colnames(qval))

  qval[[label_id]] <- samples

  qval <- qval %>%
    mutate(ind = row_number()) %>%
    pivot_longer(cols = matches("K"), names_to = "K", values_to = "Q") %>%
    group_by(ind) %>%
    mutate(max.K = which.max(Q)) %>%
    ungroup()

  df <- qval %>%
    inner_join(metadata, by = label_id) %>%
    mutate(K = as.numeric(gsub("K", "", K)))

  df$def <- paste(kk)

  admix_all <- rbind(admix_all, df)
}

# Select max.K
admix_all_f <- admix_all %>%
  mutate(!!sym(label_region) := gsub("_", " ", !!sym(label_region)),
         !!sym(label_country) := gsub("_", " ", !!sym(label_country))) %>%
  mutate(Q = as.numeric(Q)) %>%
  mutate(K = as.factor(K)) %>%
  mutate(def = as.numeric(def))


message("\nPlots saved here:")
# Generate plots
# Region
admix_all_r <- admix_all_f %>%
  mutate(!!sym(label_region) := factor(!!sym(label_region), levels = order)) %>%
  arrange(def, !!sym(label_region), max.K, desc(Q), K) %>%
  mutate(!!sym(label_id) := factor(!!sym(label_id),
         levels = unique(!!sym(label_id))))

pdf(file.path(workdir, paste0(filename, ".region.pdf")))
plot_gg_admix(admix_all_r, label_id, label_region, "region", pals, opt$axisx_angle)
invisible(dev.off())

message(file.path(workdir, paste0(filename, ".region.pdf")))

# Country
if (opt$country_code) {
  admix_all_c <- admix_all_f %>%
    mutate(!!sym(label_country) := countrycode(sourcevar = !!sym(label_country),
                                               origin = "country.name",
                                               destination = "iso3c",
                                               warn = FALSE,
                                               nomatch = NA))
} else {
  admix_all_c <- admix_all_f
}

admix_all_c <- admix_all_c %>%
    mutate(!!sym(label_region) := factor(!!sym(label_region),
           levels = order)) %>%
    arrange(!!sym(label_region)) %>%
    mutate(!!sym(label_country) := factor(!!sym(label_country),
           levels = unique(!!sym(label_country)))) %>%
    arrange(def, !!sym(label_region), !!sym(label_country), max.K, desc(Q), K) %>%
    mutate(!!sym(label_id) := factor(!!sym(label_id),
            levels = unique(!!sym(label_id))))

pdf(file.path(workdir, paste0(filename, ".country.pdf")))
plot_gg_admix(admix_all_c, label_id, label_country, "country", pals, opt$axisx_angle)
invisible(dev.off())
message(file.path(workdir, paste0(filename, ".country.pdf")))


if (!is.null(label_site)) {
  # Site
  admix_all_s <- admix_all_c %>%
  arrange(!!sym(label_region), !!sym(label_country), !!sym(label_site), max.K, desc(Q), K) %>%
  mutate(!!sym(label_id) := factor(!!sym(label_id),
         levels=unique(!!sym(label_id))))

  pdf(file.path(workdir, paste0(filename, ".site.pdf")))
  plot_gg_admix(admix_all_c, label_id, label_site, "site", pals, opt$axisx_angle)
  invisible(dev.off())
  message(file.path(workdir, paste0(filename, ".site.pdf")))
}

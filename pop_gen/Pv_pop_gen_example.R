require(scales)
require(data.table)
library(unikn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(amap)
library(ape)
library(forcats)
library(phyloch)
library(showtext)
showtext_auto()

source("~/software/malaria-hub/utils/plots.R")
source("~/software/malaria-hub/utils/helpers.R")

# Read matrix
snp <- fread("/mnt/storage9/emilia/Pvivax/DB/unfiltered_2021_02_19/Pv_gatk_05_maf_corrected_db.2021_02_19_filt.coding.ns.bi.GT.miss0.4.vqslod.filt.snps.pub.mat.bin",
              sep = "\t", header = TRUE, data.table = FALSE) %>%
              select(-c("PV14_com", "PV09_com", "SRR572650"))

# Metadata
met <- read.csv("~/Pvivax/metadata/pv_metadata_fws_final_set_public_samples_Mar2021.tsv",
                sep = "\t", header = T)

met <- as.data.frame(met)
met <- met %>% filter(sample_id %in% colnames(snp[, -c(1:3)])) %>%
               mutate(region = gsub("_", " ", region),
                      country = gsub("_", " ", country))

all(met$sample_id == colnames(snp[, -c(1:3)]))
snp_c <- snp[, -c(1:3)]

# Colors palettes
# REGION
reg_cols <- c("South Asia" = "palevioletred3",
              "East Africa" = "royalblue4",
              "South East Asia" = "royalblue",
              "Southern SEA" = "tomato2",
              "South America" = "seagreen")
order <- c("East Africa", "South America", "South Asia", "South East Asia", "Southern SEA")

# SOUTH AMERICA
SA_cols <- c("Brazil" = "firebrick2",
            "Colombia" = "yellow3",
            "Mexico" = "violetred3",
            "Guyana" = "turquoise4",
            "Panama" = "deeppink4",
            "Peru" = "orange1")
SA_order <- c("Brazil", "Colombia", "Mexico", "Guyana", "Panama", "Peru")

# ADMIXTURE COLORS
admix_cols <- unikn::usecol(c(Seeblau, Karpfenblau, Bordeaux, Seegruen, "darkorange", Pinky, "darkred", Petrol), n = 10)
names(admix_cols) <- paste0("K", c(1:10))


# Calculate distance matrix and MDS
#dist_pv <- Dist(t(snp_c), method = "manhattan", nbproc = 12)
#saveRDS(dist_pv, file = file.path(workdir, "pv_844_NO_public_manhattan_dist_Mar2021.rds"))
dist_pv <- readRDS("/mnt/storage9/emilia/Pvivax/analysis_2021_02_19/pca_tree/pv_844_NO_public_manhattan_dist_Mar2021.rds")
dist_pv <- as.matrix(dist_pv)
cmd_pv <- cmdscale(dist_pv, k = 10, eig = TRUE, x.ret = TRUE)
vars <- calc_variance_explained(cmd_pv)

# Load ADMIX K=10
workdirA <- "/mnt/storage9/emilia/Pvivax/analysis_2021_02_19/admixture"
suffix <- "Pv_gatk_05_maf_corrected_db.2021_02_19_filt.ld0.1.ns.bi.GT.miss0.4.vqslod.filt.snps.admix"
admix_Q <- read.table(file.path(workdirA, sprintf('ld0.1.bi.GT/seed12345/%s.10.Q', suffix)), as.is = T)
colnames(admix_Q) <- gsub("V", "K", colnames(admix_Q))
admix_Q <- admix_Q %>%
    mutate(ind = row_number()) %>%
    pivot_longer(cols = matches("K"), names_to = "K", values_to = "Q") %>%
    group_by(ind) %>%
    mutate(max.K = paste0("K", which.max(Q)))
fam <- read.table(file.path(workdirA, "ld0.1.bi.GT/seed12345", paste0(suffix, ".nosex")), as.is = T) %>%
    rename("sample" = "V2") %>%
    select(sample) %>%
    mutate(ind = row_number())
admix_Q <- admix_Q %>% left_join(fam) %>% ungroup()

############## PCA #############
ndfr <- combine_PCA_w_metadata(cmd_pv, met)
dfr <- ndfr %>%
  select(c(PC1, PC2, PC3, sample_id, country, region)) %>%
  left_join(admix_Q, by = c("sample_id" = "sample"))


# PCA plots with region colors
pc12 <- plot_PCA_gg(dfr, vars, "PC1", "PC2", "region", "Region: ", reg_cols, order) +
  theme_PCA()

pc13 <- plot_PCA_gg(dfr, vars, "PC1", "PC3", "region", "Region: ", reg_cols, order) +
  theme_PCA()

# PCA with ADMIXTURE
pc12A <- plot_PCA_gg(dfr, vars, "PC1", "PC2", "max.K", "", admix_cols, names(admix_cols), Lnr=4, Lnc=5) +
  theme_PCA()

pc13A <- plot_PCA_gg(dfr, vars, "PC1", "PC3", "max.K", "", admix_cols, names(admix_cols), Lnr=4, Lnc=5) +
  theme_PCA()


###### NJ Tree ######
tree_pv <- ape::nj(dist_pv)

# Specify clades
# Region clades
clades <- split(dfr$sample_id, dfr$region)
pals <- reg_cols[names(clades)]
ecol <- edge.color(tree_pv, clades, col = pals, what = "stem")
# Plot NJ tree withRegion colors
plot(tree_pv, type = "unrooted", cex = 1, lab4ut = "azial",
              edge.width = 1.5, edge.color = ecol, show.tip.label = FALSE)

# Admix clades
clades <- split(dfr$sample_id, dfr$max.K)
pals <- admix_cols[names(clades)]
ecol <- edge.color(tree_pv, clades, col = pals, what = "stem")
# Plot NJ tree with ADMIXTURE colors
plot(tree_pv, type = "unrooted", cex = 1, lab4ut = "azial",
              edge.width = 1.5, edge.color = ecol, show.tip.label = FALSE)


### ADMIXTURE ###
# Select populations with at least 5N
pop5N <- met %>% select(country) %>% plyr::count() %>% filter(freq >= 5) %>%
  select(country) %>% distinct() %>% pull()
  
admix_gg <- dfr %>%
  filter(country %in% pop5N) %>%
  mutate(K = as.numeric(gsub("K", "", K))) %>%
  select(sample_id, region, country, max.K, K, Q) %>%
  mutate(Q = as.numeric(Q)) %>%
  mutate(K = as.factor(K))

# Order
admix_gg <- admix_gg %>% mutate(region = factor(region, levels=order)) %>% arrange(region)
admix_gg <- admix_gg %>% mutate(country = factor(country, levels=unique(admix_gg$country)))

## COUNTRY ##
# Order samples for COUNTRY level
admix_gg_order <- admix_gg %>%
  filter(Q != 0) %>%
  arrange(region, country, max.K, desc(Q), K) %>%
  mutate(sample = factor(sample_id, levels = unique(sample_id)))

# Remove K from color palette names
aadmix_cols <- admix_cols
names(aadmix_cols) <- gsub("K", "", names(aadmix_cols))

# Plot
ggplot(data = admix_gg_order, aes(x = sample, y = Q, fill = K)) +
  geom_col(width = 1) +
  facet_grid(~country, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry", fill = "K") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(breaks = names(aadmix_cols),
                    values = aadmix_cols) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme_ADMIX(Fangle = 90)

## REGION ##
# Order samples for REGION level
admix_gg_order <- admix_gg %>%
  filter(Q != 0) %>%
  arrange(region, max.K, desc(Q), K) %>%
  mutate(sample = factor(sample_id, levels = unique(sample_id)))

# Plot
ggplot(data = admix_gg_order, aes(x = sample, y = Q, fill = K)) +
  geom_col(width = 1) +
  facet_grid(~region, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry", fill = "K") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(breaks = names(aadmix_cols),
                    values = aadmix_cols) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme_ADMIX(Fangle = 0, Fhjust=0.5)

## Cumulative boxplot
ggplot(admix_gg) +
    geom_bar(aes(x = country, y = ..count.., fill = max.K), position = "fill") +
    scale_fill_manual(breaks = names(admix_cols),
                      values = admix_cols) +
    labs(fill = "") +
    theme_PCA(Fsize = 8, Fangle = 90, Fhjust = 1)

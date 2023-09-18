library(tess3r)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(showtext)
showtext_auto()

## Load funcitons with themes for plots
# For PCA
theme_PCA <- function(Fsize = 10, Fcolor = "black", Fangle = 0, Fhjust = 0.5) {
    theme_classic() +
    theme(axis.text.x = element_text(size = Fsize, color = Fcolor, angle = Fangle,
          hjust= Fhjust),
          axis.text.y = element_text(size = Fsize, color = Fcolor),
          axis.title.x = element_text(size = Fsize, color = Fcolor),
          axis.title.y = element_text(size = Fsize, color = Fcolor),
          legend.text = element_text(color = Fcolor, size = Fsize - 1,
                                     margin = margin(l = 0, r = 10)),
          legend.title = element_text(color = Fcolor, size = Fsize,
                                      hjust = Fhjust),
          legend.key.size = unit(0.2, "mm"),
          legend.position = "bottom")
}

# For ADMIXTURE
theme_ADMIX <- function(Fsize = 10, Fcolor = "black", Fangle = 0, Lposition = "none", Fhjust = 1) {
    theme_minimal() +
        theme(strip.text.x = element_text(angle = Fangle, hjust = Fhjust,
                                         size = Fsize - 3, color = "black",
                                         margin = margin(t = 0.25, r = 0.5, b = 0.25, l = 0.5, unit = "lines")),
            panel.spacing.x = unit(0.2, "lines"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = Fsize, color = "black"),
            axis.text.y = element_text(size = Fsize, color = "black",
                                       angle = 90),
            legend.text = element_text(color = Fcolor, size = Fsize - 1,
                                     margin = margin(l = 0, r = 20)),
            panel.grid = element_blank(),
            legend.position = Lposition)
}
# For MAP
theme_tess3r_map <- function(Fsize = 10, Fcolor = "black", Lposition="none") {
    theme_classic() +
    theme(axis.line.x = element_line(color = Fcolor),
          axis.line.y = element_line(color = Fcolor),
          axis.text.x = element_text(size = Fsize, color = Fcolor),
          axis.text.y = element_text(size = Fsize - 2, color = Fcolor),
          axis.title.x = element_text(size = Fsize, color = Fcolor),
          axis.title.y = element_text(size = Fsize, color = Fcolor),
          plot.title = element_text(size = Fsize, color = Fcolor, hjust = 0.5),          
          legend.position = Lposition)
}

#################################################################################
workdir <- "/mnt/storage9/emilia/Leen/tess3r"

metadata <- readr::read_csv(file.path(workdir, "Metadata_877samples_withGPS_corrected.csv")) %>%
    mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude)) %>%
    mutate(country = gsub("GM_Coastal", "Gambia", country))
matbinT <- readRDS(file.path(workdir, "Pf_matbinT.mis.rds"))
coordsM <- readRDS(file.path(workdir, "Pf_coordsM.rds"))
samples <- readRDS(file.path(workdir, "Pf_877_samples.rds"))

# Load tess3 object
tess3_all <- list()
Ks <- seq(1,10)
for (K in Ks) {
    tess3_all[[K]] <- readRDS(file.path(workdir, sprintf("Pf_877_rep50_K%d.rds", 1)))
}
res5 <- readRDS(file.path(workdir, sprintf("Pf_877_rep50_K%d.rds", 5)))

# Plot error
Gettess3res(res1, K = 1)$crossentropy
Gettess3res(res1, K = 1)$rmse
Gettess3res(res5, K = 5)$rmse

# Colors
cols <- unikn::usecol(c("goldenrod2", "coral2", "darkmagenta", "darkgreen", "darkorchid3", "dodgerblue"), n = 8)
names(cols) <- paste0("K", c(1:8))

# Cumulative plots
kk <- 5
res <- res5
q.matrix <- qmatrix(res, K = kk)
p.values <- pvalue(res, K = kk)

df_matrix <- as.data.frame(q.matrix)
colnames(df_matrix) <- paste0("K", c(1:ncol(df_matrix)))
df_matrix$sample_id <- samples

df_matrix_G <- df_matrix %>% mutate(ind = row_number()) %>%
    tidyr::pivot_longer(cols = matches("K"), names_to = "K", values_to = "Q") %>%
    group_by(ind) %>%
    mutate(max.K = paste0("K", which.max(Q))) %>%
    left_join(metadata %>% select(sample_id, country, region)) %>% ungroup()
 
# Cumulative barplot max.K
# Region
cmplot <- ggplot(df_matrix_G) +
    geom_bar(aes(x = region, y =..count.., fill = max.K), position="fill") +
    scale_fill_manual(breaks = cols,
                      values = cols) +
    labs(fill="") +
    theme_PCA(Fsize = 8, Fangle = 90)


cumR <- df_matrix_G %>% group_by(region, K) %>%
     summarise(avg = sum(Q)/n())

# Cumulative barplot all K
ggplot(cumR, aes(x = region, y = avg, fill = K, group=K)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(breaks = names(cols),
                      values = cols,
                      labels = names(cols)) +
    geom_text(aes(x = region, y = avg, label = round(avg, 3)*100),
                  position = position_dodge(width = 1),
                  hjust = -0.2, vjust = 0.4, size = 3.5, angle = 90) +
    labs(x = "Region", y = "K% per region", fill="") +
    ylim(c(0, 1)) +
    theme_PCA(Fsize = 8, Fangle = 90)


### Admixture data
admix_gg <- df_matrix_G %>%
  mutate(K = as.numeric(gsub("K", "", K))) %>%
  select(sample_id, region, max.K, K, Q) %>%
  mutate(Q = as.numeric(Q)) %>%
  mutate(K = as.factor(K))

admix_gg_order <- admix_gg %>%
     mutate(region = factor(region, levels=unique(admix_gg$region))) %>%
     arrange(region) %>%
     filter(Q != 0)  %>%
     arrange(region, max.K, desc(Q), K) %>%
    mutate(sample = factor(sample_id, levels = unique(sample_id))) %>%
    mutate(K = paste0("K", K))

## Admixture plot
adm <- ggplot(data = admix_gg_order, aes(x = sample, y = Q, fill = K)) +
  geom_col(width = 1, alpha = 0.9) +
  facet_grid(~region, switch = "x", scales = "free", space = "free") +
  labs(x = "Individuals", y = "Ancestry", fill = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(add = 3)) +
  theme_ADMIX(Fangle = 90, Lposition = "bottom", Fhjust = 0.5, Fsize = 8)


### Map
my_palette <- CreatePalette(cols[1:kk])
coordsMM <- as.data.frame(coordsM) %>% left_join(metadata %>%
 dplyr::select(country, longitude, latitude) %>% distinct(), c("V1"="longitude", "V2"="latitude")) %>%
 distinct() %>% group_by(country) %>% top_n(1, V1)
map.polygon <- rworldmap::getMap(resolution = "coarse")

pl <- ggtess3Q(q.matrix, coordsM, map.polygon = map.polygon,
  col.palette = my_palette)

pl <- pl +
  geom_path(data = map.polygon, aes(x = lon, y = long, group = group), color = NA) +
    xlim(-30, 50) +
    ylim(-25, 20) +
    coord_equal() +
    geom_point(data = as.data.frame(coordsM), aes(x = V1, y = V2), color = "black", size = 0.5) + 
    geom_label_repel(data = as.data.frame(coordsMM), aes(x = V1, y = V2, label = country),
        size = 3, max.overlaps = 20) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_tess3r_map()

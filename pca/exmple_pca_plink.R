library(showtext)
library(dplyr)
library(ggplot2)
showtext_auto()

workdir <- "" # Working directory with plink files
prefix <- "" # Prefix for plink files
metadata <- "" # File path to metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "ID"))

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA ##
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling - might take a while
# saveRDS(cmd, paste0(prefix, ".dist.rds") # save to RDS format
#cmd <- readRDS(file.path(workdir, paste0(prefix, ".dist.rds"))
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Overlay region, country info
df <- as.data.frame(cmd$points, stringsAsFactors = F)
df$country <- gsub("_", " ", desc$country)
df$region <- gsub("_", " ", desc$region)
df$sample_id <- rownames(matrix)
colnames(df) <- gsub("V", "PC", colnames(df))

color_by <- "region" # specify if coloured by region or country

# Graph with PC1 an PC2
png("figure_PC12.png") # Save to PNG file
ggplot(data = df, aes(x = PC1, y = PC2,
       color = !!sym(color_by))) +
    geom_point() +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC2", " (", vars["PC2"], "%)")) +
    theme_classic() +
    theme(legend.position = "bottom")
dev.off()


# Graph with PC1 and PC3
png("figure_PC13.png") # Save to PNG file
ggplot(data = df, aes(x = PC1, y = PC3,
       color = !!sym(color_by))) +
    geom_point() +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"),
            y = paste0("PC3", " (", vars["PC3"], "%)")) +
    theme_classic() +
    theme(legend.position = "bottom")
dev.off()

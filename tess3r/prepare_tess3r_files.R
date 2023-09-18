library(tess3r)
library(readr)
library(dplyr)

wd <- "/mnt/storage9/emilia/Leen/tess3r"

# METADATA
# Make GPS coordinates uniform across site/location
# Keep columns numeric
metadata <- readr::read_csv(file.path(wd, "Metadata_925samples_withGPScoordinates_withHAFandMad.csv"), col_names = T) %>%
    select(sample_id, latitude, longitude, Site, country, region) %>%
     mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# MATRIX
# Binary matrix
# Extract header and save sample order
matbin <- data.table::fread(file.path(wd, "Pf_NECTAR_MalariaGEN_Oct2022.filt.bi.GT.miss0.4.vqslod.filt.snps.mat.bin"))
matbin_s <- matbin[ , -c(1:3)]
samples <- colnames(matbin_s)

# Check overlap
# Sample save
metadata <- metadata %>% filter(sample_id %in% samples)
matbin_s <- matbin_s %>% select(metadata$sample_id)
saveRDS(colnames(matbin_s), file.path(wd, "Pf_925_samples.rds"))

# Matrix save
# Transposed recoded mixed calls to alternative - 0.5 to 1
matbin_s[matbin_s == 0.5] <- 1
all(colnames(matbin_s) == metadata$sample_id)
matbinT <- t(matbin_s)
saveRDS(matbinT, file.path(wd, "Pf_925_matbinT.mis.rds"))

# Coordinates
# Save as matrix
coords <- metadata %>% select(longitude, latitude) %>%
      mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude))
colnames(coords) <- c("V1", "V2")
coordsM <- data.matrix(coords, rownames.force=FALSE)
saveRDS(coordsM, file.path(wd, "Pf_925_coordsM.rds"))


# Command
#cat ../Kruns.txt | xargs -I {} sh -c 'Rscript ~/software/lshtm_scripts/tess3r/run_tess3r_general.R -b Pf_1675_matbinT.mis.rds -c Pf_1675_coordsM.rds --prefix Pf_1675 --K {} --rep 50 --threads 10 && ~/cmdping' 
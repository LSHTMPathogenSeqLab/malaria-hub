library(SeqArray)
library(moimix)
library(optparse)

# Arguments ####
option_list = list(
  make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Specify main directory", metavar = "character"),
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "VCF file", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Tranform VCF to GDS file ####
if (file.exists(file.path(opt$workdir, opt$file))) {
  if (!file.exists(file.path(opt$workdir, gsub("vcf.gz", "gds", opt$file)))) {
    tryCatch({
      seqVCF2GDS(file.path(opt$workdir, opt$file),
                 file.path(opt$workdir, gsub("vcf.gz", "gds", opt$file)))
      print("Converted to GDS file")
    }, error = function(e) {
      print(e)
    })
  } else {
    print("GDS file already exists")
  }
} else {
  message("File does not exists")
}

# Calculate Fws ####
if (file.exists(file.path(opt$workdir, gsub("vcf.gz", "gds", opt$file)))) {
  gds_file <- file.path(opt$workdir, gsub("vcf.gz", "gds", opt$file))
  isolates <- seqOpen(gds_file)
  # Calculate Fws ####
  tryCatch({
    fws_all <- getFws(isolates)
    fws_df <- data.frame(sample = names(fws_all), fws = fws_all)
    name <- file.path(opt$workdir, stringr::str_split(gsub("vcf.gz", "gds", opt$file), "\\.")[[1]][1])
    write.table(fws_df, file.path(opt$workdir, sprintf("%s_moi_fws.tsv", name)),
                quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }, error = function(e) {
    print(e)
  })
} else {
  message("GDS file does not exists")
}

library(tess3r)
library(dplyr)
library(optparse)

option_list = list(
    make_option(c("-b", "--matbin"), type = "character", default = NULL,
              help = "Binary matrix RDS file",
              metavar = "character"),
    make_option(c("-c", "--coord"), type = "character", default = NULL,
                help = "Coordinates",
                metavar = "character"),
    make_option("--prefix", type = "character",
              default = format(Sys.Date(), format="%Y_%m_%d"),
              help = "File prefix",
              metavar = "character"),
    make_option("--K", type = "numeric", default = NULL,
                help = "K",
                metavar = "numeric"),
    make_option("--method", type = "character", default = "projected.ls",
              help = "Method projected.ls or ",
              metavar = "character"),
    make_option("--rep", type = "numeric", default = NULL,
              help = "Repetition",
              metavar = "numeric"),
    make_option("--threads", type = "numeric", default = NULL,
              help = "Threads",
              metavar = "numeric"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

matbinT <- readRDS(opt$matbin)
coordsM <- readRDS(opt$coord)

message("Loaded...")
tess3.obj <- tess3(X = matbinT,
                coord = coordsM,
                K = opt$K,
                method =  opt$method,
                ploidy = 1,
                rep = opt$rep,
                max.iteration = 100,
                openMP.core.num = opt$threads,
                verbose = TRUE)
message("Finished...")
saveRDS(tess3.obj, sprintf("%s_rep%d_K%d.rds", opt$prefix, opt$rep, opt$K))
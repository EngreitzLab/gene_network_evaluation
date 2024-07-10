# This script takes in a Seurat object that has already been run through
# testSoftPowers.R and uses hdWGCNA to run a WGCNA

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# path_rds: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.
# power: integer power to use for network construction

# Optional arguments
# name: Name to use for the WGCNA object. Default is the base filename of the path_rds.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object that has been run through testSoftPowers.R
# 2. Construct a co-expression network using hdWGCNA
# 3. Save the WGCNA object as an RDS file
# 
# Usage:
# Rscript runWGCNA.R -r <path_rds> -o <output_dir> -p <power> [-a <name>] [-s <seed>]

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--path_rds"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path for output RDS object"),
    make_option(c("-d", "--path_dendrogram"), type="character", default=NULL, help="Path to save dendrogram plot."),
    make_option(c("-l", "--path_module_sizes"), type="character", default=NULL, help="Path to save module sizes table."),
    make_option(c("-p", "--power"), type="integer", default=NULL, help="Soft power to use for network construction."),
    make_option(c("-a", "--name"), type="character", default=NULL, help="Name to use for the WGCNA object. Default is the base filename of the path_rds."),
    make_option(c("-s", "--seed"), type="integer", default=12345, help="Random seed to use for reproducibility.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
print(opt)

# Grab the arguments
cat("Parsing arguments")
path_rds <- opt$path_rds
path_output <- opt$path_output
path_dendrogram <- opt$path_dendrogram
path_module_sizes <- opt$path_module_sizes
power <- opt$power
name <- opt$name
seed <- opt$seed
cat("\n")

# 
cat("Using the following arguments\n")
cat(sprintf("Using RDS file %s\n", path_rds))
cat(sprintf("Saving object to %s\n", path_output))
cat(sprintf("Using soft power %d\n", power))
cat(sprintf("Using name %s\n", name))
cat(sprintf("Using seed %d\n", seed))
cat("\n")

# Imports
cat("Loading libraries\n")
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
theme_set(theme_cowplot())
options(repr.plot.width=6, repr.plot.height=6)
set.seed(opt$seed)
cat("\n")

# Actually read it
cat(sprintf("Reading RDS file %s\n", path_rds))
adata <- readRDS(path_rds)
cat("\n")

# Set-up the object
DefaultAssay(adata) <- "RNA"
adata <- SetActiveWGCNA(adata, wgcna_name=name)

# construct co-expression network:
cat("Constructing co-expression network\n")
setwd(dirname(path_dendrogram))
print(paste0("Working directory to throw the TOM object: ", getwd()))
adata <- ConstructNetwork(
  adata, 
  soft_power=power,
  use_metacells=TRUE,
  setDatExpr=FALSE,
  overwrite_tom=TRUE,
  tom_name = name # name of the topoligical overlap matrix written to disk
)
cat("\n")

# Overwrite the WGCNA object after adding the network
cat(sprintf("Saving object to %s\n", path_output))
saveRDS(adata, file=path_output)
cat("\n")

# plot the dendrogram of the genes in modules
options(repr.plot.width=12, repr.plot.height=12)
png(path_dendrogram, widt=600, height=600)
PlotDendrogram(adata, main=sprintf('%s Dendrogram', name))
dev.off()

# Save the module sizes to a dataframe
write.table(t(data.frame(table(get(name, adata@misc)$wgcna_net$colors))), path_module_sizes, sep="\t")
cat("\n")

# Done
cat("Done\n")

# This script takes in a Seurat object that has already been run through
# prepMetacells.R and uses hdWGCNA to find a soft power threshold to use

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# path_rds: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.
# group_by: Column of the metadata with groups in it. Do not include white space between the commas.
#           Column must have been used in making metacells
# groups: Comma-separated list of groups conatined within the group_by parameter to find soft-thresholds for.

# Optional arguments
# name: Name to use for the WGCNA object. Default is the base filename of the path_rds.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object that must have been processed in prepMetacells.R
# 2. Set the expression matrix based on the group_by and groups parameters
# 3. Run testSoftPowers to find a soft power threshold
# 4. Plot the results
# 5. Save the soft power threshold results
# 6. Overwrite the WGCNA object with the soft power threshold included

# Usage:
# Rscript 2_testSoftPowers.R -r <path_rds> -o <output_dir> -b <group_by> -g <groups> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--path_rds"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path for output RDS object"),
    make_option(c("-p", "--path_png"), type="character", default=NULL, help="Path to save soft power threshold plot."),
    make_option(c("-t", "--path_table"), type="character", default=NULL, help="Path to save soft power threshold table."),
    make_option(c("-b", "--group_by"), type="character", default=NULL, help="A variable in the metadata that contains groups. By default, all cells will be used"),
    make_option(c("-g", "--groups"), type="character", default=NULL, help="Comma-separated list of groups conatined within the group_by parameter to find soft-thresholds for. Cannot contain whitespace between commas!"),
    make_option(c("-a", "--name"), type="character", default=NULL, help="Name to use for the WGCNA object. Default is the base filename of the path_rds."),
    make_option(c("-s", "--seed"), type="integer", default=12345, help="Random seed to use for reproducibility.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
print(opt)

# Grab the arguments, convert any strings "None" to NULL
cat("Parsing arguments")
path_rds <- opt$path_rds
path_output <- opt$path_output
path_png <- opt$path_png
path_table <- opt$path_table
group_by <- opt$group_by
groups <- opt$groups
name <- opt$name
seed <- opt$seed
cat("\n")


# 
cat("Using the following arguments\n")
if (path_rds == "None") {
    stop("path_rds must be set.")
}
cat(sprintf("path_rds: %s\n", path_rds))
if (path_output == "None") {
    stop("path_output must be set.")
}
cat(sprintf("path_output: %s\n", path_output))
if (path_png == "None") {
    path_png <- NULL
}
cat(sprintf("path_png: %s\n", path_png))
if (path_table == "None") {
    path_table <- NULL
}
cat(sprintf("path_table: %s\n", path_table))
if (group_by == "None") {
    group_by <- NULL
}
cat(sprintf("group_by: %s\n", group_by))
if (groups == "None") {
    groups <- NULL
}
cat(sprintf("Using groups %s\n", groups))
if (name == "None") {
    name <- NULL
}
cat(sprintf("name: %s\n", name))
if (seed == "None") {
    seed <- NULL
}
cat(sprintf("seed: %s\n", seed))
cat("\n")

# Ensure that groups is set if groupby is set
if(!is.null(group_by) && is.null(groups)) {
    stop("If group_by is set, groups must be set.")
}

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

# transpose the matrix, taking care of the
cat("Setting the expression matrix for running testSoftPowers\n")
adata <- SetDatExpr(
    adata,
    group.by=group_by,
    group_name=groups,
    assay="RNA",
    use_metacells=TRUE,
    wgcna_name=name,
    slot="data"
)
dims <- dim(get(name, adata@misc)$datExpr)

cat(sprintf("Dimensions of the matrix are %s\n", paste0(dims, collapse="x")))
cat("\n")

# Test different soft powers:
cat("Running testSoftPowers\n")
adata <- TestSoftPowers(
  adata,
  use_metacells = TRUE,  # this is the default, I'm just being explicit
  setDatExpr = FALSE  # set this to FALSE since we did this above
)
cat("\n")

# plot the results:
cat("Plotting the soft power threshold results\n")
options(repr.plot.width=12, repr.plot.height=12)
png(path_png, widt=600, height=600)
plot_list <- PlotSoftPowers(adata)
wrap_plots(plot_list, ncol=2)
dev.off()

# Save the soft power threshold results
power_table <- GetPowerTable(adata)
power <- power_table$Power[which(power_table$SFT.R.sq > 0.80)[1]]
cat(sprintf("Automatic soft power %s\n", power))
write.table(power_table, path_table, sep="\t")
cat("\n")

# Write the WGCNA object with the soft power threshold included
cat(sprintf("Saving object to %s\n", path_output))
saveRDS(adata, file=path_output)
cat("\n")

# Done
cat("Done\n")

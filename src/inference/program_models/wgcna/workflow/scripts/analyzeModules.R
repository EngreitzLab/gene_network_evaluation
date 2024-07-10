# This script takes in a Seurat object that has already been run through
# constructNetwork.R and uses hdWGCNA analyze the co-expression modules in the network.

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# path_rds: Path to RDS file containing processed Seurat object.
# output_dir: Path to directory to output final object and intermediate files.

# Optional arguments
# harmonize_by: Comma-separated list of metadata fields to perform harmonization of module eigengenes with, no white space between commas!
# group_by: column in seurat_obj@meta.data containing grouping info, ie clusters or celltypes
# groups: Comma-separated list name of the group(s) in harmonize_by to use for kME calculation
# name: Name to use for the WGCNA object. Default is the base filename of the path_rds.
# seed: Random seed to use for reproducibility.

# The script will:
# 1. Read in the Seurat object
# 2. Compute module eigengenes
# 3. Compute kMEs
# 4. Write out the module eigengenes and kMEs to a file
# 5. Write out the Seurat object that is ready for downstream analysis

# Usage:
# Rscript 4_analyzeModules.R -r <path_rds> -o <output_dir> -b <harmonize_by> -g <harmonize_by> -a <name> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-r", "--path_rds"), type="character", default=NULL, help="Path to RDS file containing processed Seurat object."),
    make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path for output RDS object"),
    make_option(c("-p", "--path_MEs"), type="character", default=NULL, help="Path to save module eigengenes table."),
    make_option(c("-l", "--path_modules"), type="character", default=NULL, help="Path to save module membership table."),
    make_option(c("-q", "--harmonize_by"), type="character", default=NULL, help="Comma-separated list of variables to perform harmonization of module eigengenes with, no white space between commas!"),
    make_option(c("-b", "--group_by"), type="character", default=NULL, help="Column in meta.data containing grouping info, i.e. clusters or celltypes"),
    make_option(c("-g", "--groups"), type="character", default=NULL, help="Comma-separated list of groups conatined within the harmonize_by parameter to find soft-thresholds for."),
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
path_MEs <- opt$path_MEs
path_modules <- opt$path_modules
harmonize_by <- opt$harmonize_by
group_by <- opt$group_by
groups <- opt$groups
name <- opt$name
seed <- opt$seed
cat("\n")

# Parse the arguments
cat(sprintf("Using RDS file %s\n", path_rds))
cat(sprintf("Saving object to %s\n", path_output))
cat(sprintf("Saving module eigengenes to %s\n", path_MEs))
if (harmonize_by == "None") {
    harmonize_by <- NULL
}
if (!is.null(harmonize_by)) {
    harmonize_by <- strsplit(harmonize_by, ",")[[1]]
    cat(sprintf("Harmonizing by %s\n", paste0(harmonize_by, collapse=", ")))
} else {
    cat("Not harmonizing\n")
}
if (group_by == "None") {
    group_by <- NULL
}
if (groups == "None") {
    groups <- NULL
}
if(is.null(group_by)) {
    cat("Using all cells during kME calculation\n")
    if(!is.null(groups)) {
        stop("Cannot include groups argument if not using group_by for kME calculation")
    }
} else {
    cat(sprintf("Using only selected cells for kME calculation using %s column\n", group_by))
    if (is.null(groups)) {
        stop("Must provide groups if using group_by for kME calculation")
    }
    groups <- strsplit(groups, ",")[[1]]
    cat(sprintf("Using groups %s\n", paste0(groups, collapse=", ")))
}
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

# Scaling data
cat("Scaling data\n")
if (is.null(VariableFeatures(adata))) {
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
}
adata <- ScaleData(adata, features=VariableFeatures(adata))
cat("\n")

# Analyze these modules
cat("Computing module eigengenes\n")
adata <- ModuleEigengenes(
    adata, 
    group.by.vars=harmonize_by,
    verbose=TRUE,
    assay="RNA",
    wgcna_name=name
)
if (!is.null(harmonize_by)) {
    hMEs <- GetMEs(adata, harmonized=TRUE)
    path_hMEs <- gsub(".tsv", "_harmonized.tsv", path_MEs)
    write.table(hMEs, path_hMEs, sep="\t")
    cat("Computing harmonized kMEs\n")
    adata <- ModuleConnectivity(
        adata,
        group.by = group_by,
        group_name = groups,
        assay="RNA",
        slot="data",
        harmonized=TRUE,
        wgcna_name=name
    )
} else {
    cat("Computing kMEs")
    adata <- ModuleConnectivity(
        adata,
        group.by = group_by,
        group_name = groups,
        assay="RNA",
        slot="data",
        harmonized=FALSE,
        wgcna_name=name
    )
}
MEs <- GetMEs(adata, harmonized=FALSE)
write.table(MEs, path_MEs, sep="\t")
cat("\n")

# Save modules
cat("Saving module membership to file\n")
write.table(GetModules(adata), path_modules, sep="\t")

# Overwrite the WGCNA object with the module fun included
cat(sprintf("Saving object to %s\n", path_output))
saveRDS(adata, file=path_output)
cat("\n")

# Done
cat("Done\n")

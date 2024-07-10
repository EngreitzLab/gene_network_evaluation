# This script takes in a processed Seurat object and uses hdWGCNA to 
# to prepare metacells for downstream analysis.

# This script uses optparse to parse arguments for specifying fragment files, output directory, 
# and optional parameters for data processing.

# Required arguments:
# path_data: Path to h5mu file
# path_output: Path for output RDS object
# path_metacell_rds: Path to RDS file containing metacell object.
# path_metacell_stats: Path to CSV file containing metacell stats.
# group_by: Comma-separated list of variables to use for grouping cells. Do not include white space between the commas.

# Optional arguments
# nearest_neighbors: number of nearest neighbors to use for creating metacells, by default will use 25
# genes: either a float between 0 and 1 or the number of variable genes to use, by default will use 5% of genes
# reduction: dimensionality reduction to use for KNN, by default will look for "harmony" in the reductions
# name: name to use for the WGCNA object. Default is the base filename of the path_data.
# seed: random seed to use for reproducibility, by default will use 12345

# The script will:
# 1. Load in the rds file with the processed Seurat object
# 2. Plot the umap with the first grouping variable
# 3. Set-up the object for WGCNA using either variable genes or genes expressed in at least x% of cells
# 4. Construct metacells with x nearest neighbors and grouped by x. Note that the ident group is the first grouping variable.
# 5. Normalize the metacell expression matrix
# 6. Save the metacell stats
# 7. Save the metacell object
# 8. Save the object for the next step

# Usage:
# Rscript --vanilla prepMetacells.R -r <path_data> -o <path_output> -m <path_metacell_rds> -t <path_metacell_stats> -b <group_by> -n <nearest_neighbors> -g <genes> -u <reduction> -a <name> -s <seed>

# Argument parsing
suppressMessages(library(optparse))
option_list <- list(
    make_option(c("-p", "--path_data"), type="character", default=NULL, help="Path to h5mu."),
    make_option(c("-o", "--path_output"), type="character", default=NULL, help="Path for output RDS object"),
    make_option(c("-m", "--path_metacell_rds"), type="character", default=NULL, help="Path to RDS file containing metacell object."),
    make_option(c("-t", "--path_metacell_stats"), type="character", default=NULL, help="Path to CSV file containing metacell stats."),
    make_option(c("-d", "--path_umap"), type="character", default=NULL, help="Path to save UMAP plot."),
    make_option(c("-b", "--group_by"), type="character", default=NULL, help="Comma-separated list of variables to use for grouping cells. Do not include white space between the commas."),
    make_option(c("-n", "--nearest_neighbors"), type="integer", default=25, help="Number of nearest neighbors to use for creating metacells."),
    make_option(c("-g", "--genes"), type="numeric", default=0.05, help="Either a float between 0 and 1 or the number of variable genes to use."),
    make_option(c("-c", "--min_cells"), type="integer", default=100, help="Minimum number of cells in a metacell."),
    make_option(c("-q", "--target_metacells"), type="integer", default=1000, help="Target number of metacells."),
    make_option(c("-x", "--max_shared"), type="integer", default=10, help="Maximum number of shared cells between two metacells."),
    make_option(c("-u", "--reduction"), type="character", default="harmony", help="Dimensionality reduction to use for KNN."),
    make_option(c("-a", "--name"), type="character", default=NULL, help="Name to use for the WGCNA object. Default is the base filename of the path_data."),
    make_option(c("-s", "--seed"), type="integer", default=12345, help="Random seed to use for reproducibility.")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
print(opt)

# Grab the arguments
cat("Parsing arguments")
path_data <- opt$path_data
path_output <- opt$path_output
path_metacell_rds <- opt$path_metacell_rds
path_metacell_stats <- opt$path_metacell_stats
path_umap <- opt$path_umap
group_by <- strsplit(opt$group_by, ",")[[1]]
nearest_neighbors <- as.numeric(opt$nearest_neighbors)
target_metacells <- as.numeric(opt$target_metacells)
max_shared <- as.numeric(opt$max_shared)
genes <- as.numeric(opt$genes)
min_cells <- as.numeric(opt$min_cells)
reduction <- opt$reduction
name <- opt$name
seed <- opt$seed
cat("\n")

# Grab a name for this experiment
if (is.null(name)) {
    name <- sprintf("%s-genes_%s-neighbors", genes, nearest_neighbors)
}
cat(sprintf("Reading data from %s\n", path_data))
cat(sprintf("Saving output to %s\n", path_output))
cat(sprintf("Using name %s\n", name))
cat(sprintf("Grouping metacells by %s\n", paste(group_by, collapse=", ")))
cat(sprintf("Using %s nearest neighbors\n", nearest_neighbors))
cat(sprintf("Using %s genes\n", genes))
cat(sprintf("Using %s min cells\n", min_cells))
cat(sprintf("Using %s reduction\n", reduction))
cat(sprintf("Using %s target metacells\n", target_metacells))
cat(sprintf("Using %s max shared\n", max_shared))
cat(sprintf("Using %s seed\n", seed))
cat("\n")

# Imports
cat("Loading libraries\n")
suppressMessages(library(Seurat))
suppressMessages(library(MuData))
suppressMessages(library(MultiAssayExperiment))
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
cat(sprintf("Reading h5mu file %s\n", path_data))
mdata <- readH5MU(path_data)
counts <- assays(experiments(mdata)$rna)$counts
metadata <- as.data.frame(colData(mdata))
adata <- Seurat::CreateSeuratObject(
    counts = counts, 
    meta.data = metadata
)
cat("\n")

# Seurat guided clustering
cat("Running Seurat guided clustering\n")
adata <- Seurat::NormalizeData(adata)
adata <- Seurat::FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
adata <- Seurat::ScaleData(adata)
adata <- Seurat::RunPCA(adata, features = Seurat::VariableFeatures(adata))
adata <- Seurat::RunUMAP(adata, reduction = "pca", dims = 1:30)
adata <- Seurat::FindNeighbors(adata, dims = 1:30)
adata <- Seurat::FindClusters(adata, resolution = 0.5)
cat("\n")

# Plot and save umap with grouping 1
cat("Plotting umap with grouping 1\n")
png(path_umap, width=800, height=800)
DimPlot(adata, group.by=group_by[1], reduction="umap", label=FALSE) + umap_theme() + ggtitle(sprintf("%s", group_by[1]))
dev.off()
cat("\n")

# Set-up the object
DefaultAssay(adata) <- "RNA"

# Set-up a Seurat object for WGCNA
cat("Setting up Seurat object for WGCNA\n")
print(class(genes))
if (is.double(genes)) {
    print(sprintf("Using genes expressed in at least %s%% of cells", genes * 100))
    adata <- SetupForWGCNA(
        adata,
        gene_select = "fraction", # the gene selection approach
        fraction = genes, # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = name # the name of the hdWGCNA experiment
    )
} else if (genes >= 1) {
    print("Using variable genes")
    adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
    adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = genes)
    adata <- SetupForWGCNA(
        adata,
        gene_select = "variable", # the gene selection approach
        wgcna_name = name # the name of the hdWGCNA experiment
    )
}
cat(sprintf("Using %s genes\n", length(get(name, adata@misc)$wgcna_genes)))
cat("\n")

# construct metacells n each group
cat(sprintf("Constructing metacells with %s nearest neighbors and grouped by %s\n", nearest_neighbors, paste(group_by, collapse=", ")))
adata <- MetacellsByGroups(
  seurat_obj = adata,
  group.by = group_by, # specify the columns in adata@meta.data to group by
  reduction = reduction, # select the dimensionality reduction to perform KNN on
  k = nearest_neighbors, # nearest-neighbors parameter
  min_cells = min_cells, # minimum number of cells in a metacell
  max_shared = max_shared, # maximum number of shared cells between two metacells
  target_metacells = target_metacells, # target number of metacells
  ident.group = group_by[1], # set the Idents of the metacell seurat object
  assay="RNA",  # default
  slot="counts",  # default
)
cat("\n")

# normalize metacell expression matrix
adata <- NormalizeMetacells(adata)
cat("\n")

# Save the object for the next step!
cat(sprintf("Saving object to %s\n", path_output))
saveRDS(adata, file=path_output)
cat("\n")

# Save the metacell stats
cat(sprintf("Saving metacell stats to %s\n", path_metacell_stats))
write.csv(
    x=get(name, adata@misc)$wgcna_params$metacell_stats, 
    file=path_metacell_stats,
    quote=FALSE, 
    row.names=FALSE
)
cat("\n")

# Save the metacell object
cat(sprintf("Saving metacell object to %s\n", path_metacell_rds))
saveRDS(get(name, adata@misc)$wgcna_metacell_obj, file=path_metacell_rds)
cat("\n")

# Done
cat("Done!\n")

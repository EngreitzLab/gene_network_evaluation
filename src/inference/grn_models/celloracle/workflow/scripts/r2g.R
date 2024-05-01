library(cicero)
library(monocle3)
library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
binarize <- as.logical(args[8])
dim_reduction_key <- args[9]
k <- as.numeric(args[10])
window <- as.numeric(args[11])
path_all_peaks <- args[12]
path_connections <- args[13]
path_cicero_out <- args[14]
seed <- as.numeric(args[15])

# Read genome
if (organism == 'human'){
    genome <- read.table('resources/genome_sizes/human.txt')
} else {
    genome <- read.table('resources/genome_sizes/mouse.txt')
}

# Process mudata
indata <- H5Fopen(path_data)
indices <- indata$mod$atac$layers$counts$indices
indptr <- indata$mod$atac$layers$counts$indptr
data <- as.numeric(indata$mod$atac$layers$counts$data)
barcodes <- indata$mod$atac$obs$`_index`
peaks <- indata$mod$atac$var$`_index`
h5closeAll()

# Build sparse matrix and binarize if specified
indata <- Matrix::sparseMatrix(i=indices, p=indptr, x=data, index1 = FALSE)
if (binarize){
    indata@x[indata@x > 0] <- 1
}

# Format cell info
cellinfo <- data.frame(row.names=barcodes, cells=barcodes)

# Format peak info
peaks <- gsub(":", "-", peaks)
peakinfo <- data.frame(row.names=peaks, site_name=peaks)
peakinfo <- tidyr::separate(data = peakinfo, col = 'site_name', into = c("chr", "bp1", "bp2"), sep = "-", remove=FALSE)

# Add names
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(
    new_cell_data_set(
        indata,
        cell_metadata = cellinfo,
        gene_metadata = peakinfo
    )
)

# Data preprocessing
set.seed(seed)
if (!is.null(dim_reduction_key)){
    print("Aggregating into metacells")
    if (dim_reduction_key == "umap" || dim_reduction_key == "lsi"){
        print(sprintf("Running %s dimensionality reduction with Cicero", dim_reduction_key))
        input_cds <- monocle3::detect_genes(input_cds)
        input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
        input_cds <- estimate_size_factors(input_cds)
        input_cds <- preprocess_cds(input_cds, method = "LSI")
        if (dim_reduction_key == "lsi"){
            coords <- reducedDims(input_cds)$LSI
        } else if (dim_reduction_key == "umap"){
            input_cds <- reduce_dimension(
                input_cds, 
                reduction_method = 'UMAP', 
                preprocess_method = "LSI"
            )
            coords <- reducedDims(input_cds)$UMAP
        }
    } else {
        print(sprintf("Getting precomputed %s coordinates from h5", dim_reduction_key))
        indata <- H5Fopen(path_data)
        coords <- t(indata$mod$atac$obsm[[dim_reduction_key]])
        rownames(coords) <- barcodes
        h5closeAll()
    }
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = coords, k=k)
} else {
    print("No dimensionality reduction specified, not using metacells")
    cicero_cds <- input_cds
}

# Calculate distance parameter
print("Calculating distance_parameter value")
distance_parameters <- estimate_distance_parameter(
    cicero_cds,
    window=window,
    maxit=100,
    sample_num=100,
    distance_constraint = round(window/2),
    distance_parameter_convergence=1e-22,
    genomic_coords=genome
)
mean_distance_parameter <- mean(unlist(distance_parameters))

# Build models
print("Building models")
cicero_out <- generate_cicero_models(
    cicero_cds,
    distance_parameter=mean_distance_parameter,
    window=window,
    genomic_coords=genome
)

# Assemble connections
print("Assembling connections")
conns <- assemble_connections(cicero_out, silent=FALSE)

# Save
print("Saving")
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(path_all_peaks))
write.csv(x = conns, file = file.path(path_connections))
saveRDS(object = cicero_out, file = file.path(path_cicero_out))

# Done
print("Successfully finished running peak correlation with Cicero.")

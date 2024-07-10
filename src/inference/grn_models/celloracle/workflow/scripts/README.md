# Scripts

## download_genomes.py: Download specific genome versions for use with `genomepy`.

**Usage**:
```bash
python download_genomes.py -d [GENOME_DIR] -g [GENOME]

-h, --help: Displays help information and exits.
-d GENOME_DIR, --genome_dir GENOME_DIR: Specify the directory to save the genome.
-g GENOME, --genome GENOME: Specify the genome version to download (e.g., hg38, mm10).
```

**Examples**:<br>
Download the human genome version hg38: 
```bash 
python download_genomes.py --genome_dir /path/to/genome_directory --genome hg38
```
Download the mouse genome version mm10: 
```bash 
python download_genomes.py --genome_dir /path/to/genome_directory --genome mm10
```

## r2g.R: Identify coaccessible peaks with Cicero for building region to gene links

**Usage**:
```bash
Rscript cicero_analysis.R [PATH_DATA] [CHROMSIZES] [BINARIZE] [DIM_REDUCTION_KEY] [K] [WINDOW] [PATH_ALL_PEAKS] [PATH_CONNECTIONS] [PATH_CICERO_OUT] [SEED]

- PATH_DATA: Path to input MuData object, this is not modified in this step. Note, the object must contain an 'atac' modality with a 'counts' layer that is of type scipy.sparse.csr.csr_matrix. A matrix of type numpy.ndarray in that slot will not work.
- CHROMSIZES: Path to chromosome sizes file for the genome. Usually can find in same directory as genome fasta file after running `download_genomes.py`.
- BINARIZE: Logical value (TRUE or FALSE) to specify if the data should be binarized.
- DIM_REDUCTION_KEY: Dimensionality reduction key in mdata.mod['atac'].obsm to use for Cicero kNN aggregation. You can specify any key in the mdata.mod['atac'].obsm slot, e.g., 'X_pca', 'X_umap', 'X_tsne', etc. You can also specify 'umap' or 'lsi' if you want to use Cicero's internal dimensionality reduction implementation. If null is provided, Cicero will use single cells for inference of coaccessibility.
- K: Number of nearest neighbors for Cicero kNN aggregation. Ignored if DIM_REDUCTION_KEY is null.
- WINDOW: Only allow peaks within this distance to be considered for oaccessibility analysis. 
- PATH_ALL_PEAKS: Path to all peaks file output. This is a simple txt file containing all peaks used in the analysis.
- PATH_CONNECTIONS: Path to Cicero connections file output.
- PATH_CICERO_OUT: Output path for Cicero cell dataset in rds format.
- SEED: Seed for random number generation to ensure reproducibility.
```

**Example**:<br>
```bash
Rscript r2g.R /path/to/data /path/to/chromsizes.tsv FALSE X_pca 30 1000 /path/to/all_peaks.csv /path/to/cicero_connections.csv /path/to/cicero_cds.rds 1234
```

> 🚨 **Important Note:**
>
> To run this step in the same manner as the CellOracle tutorial, use the following parameters:
>
```yaml
binarize: true  # whether to binarize the peak matrix
dim_reduction_key: umap  # uses Cicero's internal dimensionality reduction implementation
k: 50  # number of neighbors for kNN
window: 500000  # window size for Cicero
```

## r2g.py: Build region to gene links using Cicero coaccessibility scores.

**Usage**:
```bash
python r2g.py -d [PATH_DATA] -a [ALL_PEAKS] -c [CONNECTIONS] -g [GENOME] -t [COACCESS_THR] -o [PATH_OUT]

-h, --help: Displays help information and exits.
-d PATH_DATA, --path_data PATH_DATA: Path to input MuData object, this is not modified in this step.
-a ALL_PEAKS, --all_peaks ALL_PEAKS: Path to scATAC-seq peak list output from r2g.R.
-c CONNECTIONS, --connections CONNECTIONS: Path to Cicero coaccessibility scores output from r2g.R.
-g GENOME, --genome GENOME: Version of the genome to use, e.g., hg38, mm10.
-t COACCESS_THR, --coaccess_thr COACCESS_THR: Filter out Cicero connections below this threshold.
-o PATH_OUT, --path_out PATH_OUT: Path to output r2g.csv containing region to gene links.
```

**Example**:<br>
```bash
python r2g.py -d /path/to/data -a /path/to/all_peaks.csv -c /path/to/cicero_connections.csv -g hg38 -t 0.5 -o /path/to/output/r2g.csv
```

## tf2r.py: Build region to TF (transcription factor) links using GimmmeMotifs motif scanning.

**Usage**:
```bash
python tf2r.py -d [PATH_DATA] -r [PATH_R2G] -gd [GENOME_DIR] -g [GENOME] -f [FPR] -b [BLEN] -t [TFB_THR] -o [PATH_OUT]

-h, --help: Displays help information and exits.
-d PATH_DATA, --path_data PATH_DATA: Path to input MuData object. This is not modified in this step.
-r PATH_R2G, --path_r2g PATH_R2G: Path to region to gene links output from r2g.py.
-gd GENOME_DIR, --genome_dir GENOME_DIR: Path to the genomepy reference genome to use, must be downloaded first using download_genomes.py.
-g GENOME, --genome GENOME: Version of the genome to use, e.g., hg38, mm10.
-f FPR, --fpr FPR: False positive rate for GimmmeMotifs motif scanning.
-b BLEN, --blen BLEN: Background sequence length for GimmmeMotifs motif scanning.
-t TFB_THR, --tfb_thr TFB_THR: Threshold for filtering TF binding predictions.
-p THREADS, --threads THREADS: Number of threads to use for motif scanning.
-ti PATH_TFINFO, --path_tfinfo PATH_TFINFO: Path to output h5 file with raw motif scanning results.
-o PATH_OUT, --path_out PATH_OUT: Path to output tf2r.csv containing TF to region links.
```

**Example**:<br>
```bash
python tf2r.py -d /path/to/data -r /path/to/r2g.csv -gd /path/to/genome_directory -g hg38 -f 0.02 -b 200 -t 5 -p 4 -o /path/to/output/tf2r.csv
```

## grn.py: Build TF to gene links for each cluster using CellOracle's bagging ridge regression.

**Usage**:
```bash
python grn.py -d [PATH_DATA] -r [PATH_R2G] -t [PATH_TF2R] -c [CLUSTER_KEY] [-l [LAYER]] -a [ALPHA] -b [BAGGING_NUMBER] -o [PATH_OUT]

-h, --help: Displays help information and exits.
-d PATH_DATA, --path_data PATH_DATA: Path to input MuData object, this is not modified in this step.
-r PATH_R2G, --path_r2g PATH_R2G: Path to region to gene links output from r2g.py.
-t PATH_TF2R, --path_tf2r PATH_TF2R: Path to TF to region links output from tf2r.py.
-b PATH_BASE_GRN, --path_base_grn PATH_BASE_GRN: Path to base GRN containing initial TF to gene links. Will ignore r2g and tf2r if provided. This is useful if you want to bring a GRN from a external pipeline. Expected to be in parquet format as specified in the CellOracle tutorial.
-c CLUSTER_KEY, --cluster_key CLUSTER_KEY: Key in mdata.obs that contains cluster information. Will be used to build GRN for each cluster.
-l LAYER, --layer LAYER (Optional): Layer in mdata.mod['rna'].layers to use for regression. If not provided, will log normalized X.
-a ALPHA, --alpha ALPHA: Alpha value for LASSO regression.
-b BAGGING_NUMBER, --bagging_number BAGGING_NUMBER: Number of bagging iterations for LASSO regression.
-s SCALE, --scale SCALE (Optional): Whether to scale data before regression. Default is True.
-o PATH_OUT, --path_out PATH_OUT: Path to output grn.csv containing TF to gene links.
```

**Example**:<br>
```bash
python grn.py -d /path/to/data -r /path/to/r2g.csv -t /path/to/tf2r.csv -c cluster_key -l counts -a 10 -b 20 -o /path/to/output/grn.csv
```

> 🚨 **Important Note:**
>
> How you normlize the matrix prior to regression can have a big impact on the results. We have found that the following parameters best replicate the CellOracle tutorial:
>
```yaml
layer: imputed_counts  # use imputed counts for regression as in the CellOracle tutorial
alpha: 10  # alpha for ridge regression regularization strength
bagging_number: 20  # number of samples for bagging
scale: true  # this will scale each cluster's expression matrix before regression
```

## post.py: Post process outputs of pipeline for evaluation.

**Usage**:
```bash
python post.py -d [PATH_DATA] -r [PATH_R2G] -t [PATH_TF2R] -g [PATH_GRN] -o [PATH_OUT]

-h, --help: Displays help information and exits.
-d PATH_DATA, --path_data PATH_DATA: Path to input MuData object. This is not modified in this step.
-r PATH_R2G, --path_r2g PATH_R2G: Path to region to gene links output from r2g.py.
-t PATH_TF2R, --path_tf2r PATH_TF2R: Path to TF to region links output from tf2r.py.
-g PATH_GRN, --path_grn PATH_GRN: Path to the GRN output from grn.py.
-o PATH_OUT, --path_out PATH_OUT: Path to output the post-processed GRN data.
```

**Example**:<br>
```bash
python post.py -d /path/to/data -r /path/to/r2g.csv -t /path/to/tf2r.csv -g /path/to/grn.csv -o /path/to/output/celloracle.h5mu
```

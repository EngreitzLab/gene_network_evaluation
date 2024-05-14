author: Adam Klie <br>
email: aklie@ucsd.edu <br>
date: 2024-02-21

# TODO


# SCENIC+ GRN inference
This directory contains a snakemake workflow for running the SCENIC+ pipeline for inference of gene regulatory networks (GRNs) from single-cell multiome data (paired scRNA-seq and scATAC-seq).

# Quick start
1. Modify the config in `config/config.yaml` or create a new one with the same structure
- [ ] Point to the correct input file (`input_loc`) (see [Expected input](#expected-input))
- [ ] Change the `organism` to the correct species (currently only human and mouse are supported)
- [ ] Change output directory (`outdir`) to where you want the output to be saved. This includes all intermediate files and the final MuData object (see [Output](#output))
- [ ] Modify the path to the SCENIC+ singularity container (`singularity_image`) (see [Environment](#environment) for more details)
- [ ] Modify the scratch directory (`scratch`) to where you want the temporary files to be saved
- [ ] Choose the number of threads to use (`threads`) based on your system
2. Run the pipeline
```bash
snakemake --cores 1 --configfile /path/to/config.yaml tri.csv  # need to figure out how to not have to pass in tri.csv
snakemake --use-singularity --cores 1 tri.csv  # Use singularity container
```

# Expected input
* `mdata.h5mu` — MuData object in h5mu format containing single-cell multiome data. (see [MuData documentation](https://mudata.readthedocs.io/en/latest/))
    * The MuData object MUST contain the following:
        * `atac` in `mod` -- h5ad for scATAC-seq data
            * `layers["counts"]` — a sparse matrix of raw fragment counts
            * `var_names` — a list of region names in 'chr-start-end' format
        * `rna` in `mod` -- h5ad for scRNA-seq data
            * `layers["counts"]` — a sparse matrix of raw UMI counts
        * `obs` — a dataframe of cell metadata
            * `"cell_identity"` — a column of cell type annotations (specify `cluster_key` to match in the config file)
        * `var` — a dataframe of gene metadata

# Workflow
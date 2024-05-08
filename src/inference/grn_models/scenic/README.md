author: Adam Klie <br>
email: aklie@ucsd.edu <br>
date: 2024-05-02

# TODO
- [ ] Test on a diverse set of datasets
- [ ] Add to IGVF docs

# SCENIC GRN inference

# Quick start
1. Modify the config in `config/config.yaml` or create a new one with the same structure
- [ ] Point to the correct input file (`input_loc`) (see [Expected input](#expected-input))
- [ ] Change the `genome` to the correct species and version (e.g. `hg38`)
- [ ] Change output directory (`outdir`) to where you want the output to be saved. This will include all intermediate files and the final MuData object (see [Output](#output))
- [ ] Only if you are using the `--use-singularity` flag. Modify the path to the scenic singularity container (`singularity_image`) (see [Environment](#environment) for more details)
- [ ] Modify the scratch directory (`scratchdir`) to where you want the temporary genome files o be saved
- [ ] Choose the number of threads to use (`threads`) based on your system
- [ ] Modify other parameters as needed (see [Parameters](#parameters) for more details

2. Run the pipeline
```bash
snakemake --cores 1 outdir/scenic.h5mu --configfile /path/to/config.yaml
snakemake --use-singularity --cores 1 outdir/scenic.h5mu --configfile /path/to/config.yaml # Use singularity container
```

# Expected input
* `mdata.h5mu` — MuData object in h5mu format containing single-cell multiome data. (see [MuData documentation](https://mudata.readthedocs.io/en/latest/))
    * The MuData object MUST contain the following:
        * `rna` in `mod` -- h5ad for scRNA-seq data
            * `layers["counts"]` — a sparse matrix of raw UMI counts
        * `obs` — a dataframe of cell metadata
            * `"cell_identity"` — a column of cell type annotations (specify `cluster_key` to match in the config file)
        * `var` — a dataframe of gene metadata

# Workflow
The SCENIC method consists of the following steps:

# Pipeline outputs
* `tf2r.csv` — TF to region links
```bash
cre,tf,score,pval
chr1-107964928-107965436,SP4,306.97722557821425,
chr1-107964928-107965436,E2F2,39.214405253215055,
chr1-107964928-107965436,KLF6,33.173333506316666,
chr1-107964928-107965436,TCF7L2,21.085326175707472,
chr1-107964928-107965436,SP1,20.0,
```
> **Note**
> The `score` column is the motif score for the TF binding to the region. The `pval` column is NaN for scenic as GimmmeMotifs does not provide p-values. Rigth now, I believe we are allowing for both direct and indirect TF binding to the region. We may want to keep track of this in the future.

* `grn.csv` — TF to gene links
```bash
tf,gene,score,pval,cluster
ARID5B,FAM13A,0.3710782,1.2700890408360003e-05,HSC
ARID5B,SAMD9,0.13423152,0.024985336730666132,Proerythroblast
ARID5B,FAM13A,0.11488605,0.04676196383997795,Proerythroblast
ARID5B,DST,0.006444159,0.860277418467874,MK/E prog
ARID5B,SAMD9,-0.02657091,0.6167909876627851,HSC
```
> **Note**
> The `score` column is the coefficient of the TF in the Bayesian Ridge regression. The `pval` column is the p-value of the coefficient across the bagging. The `cluster` column is the cell type cluster that the GRN was built in.

# More details

## Parameter recommendations
Coming soon

## Environment
You can find a `.def` file for building the singularity container in the `envs/` directory:
```bash
singularity build --remote envs/scenic.sif workflow/envs/scenic.def 
```

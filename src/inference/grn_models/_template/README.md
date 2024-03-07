author: Adam Klie <br>
email: aklie@ucsd.edu <br>
date: 2024-02-28

# Adding a new GRN model

# TODO
- [ ] How are we handling fragment files?
- [ ] What exactly should the format of the edges be?
- [ ] How do we adapt this to better facilitate a jamboree?

## Implementing your model
If your model is simple to run (i.e. can be run with a few functions), we recommend formatting your implementation as described in `src/inference/program_models/_template`. 

However, GRN methods are typically more complex and require multiple steps and scripts (sometimes across programming languages) to run. In this case, we recommend building your model as a [snakemake pipeline](##Building-a-snakemake-pipeline-for-your-model)

Eiter way, your model should be implemented in a subdirectory of `src/inference/grn_models` and should adhere to the following:

1. GRN models are expected to take in a [MuData object in `h5mu` format](https://mudata.readthedocs.io/en/latest/). You can expect this object to contain one or more single-cell data modalities with the following properties:
    * The MuData object *WILL ALWAYS* contain the following:
        * `obs` — a dataframe of cell metadata
            * `"celltype"` — a column of cell type annotations
        * `var` — a dataframe of feature metadata for the modalties contained in the MuData object

    * The MuData object WILL contain some combination of the following depending on how the data was generated:
        * `atac` in `mod` -- h5ad for scATAC-seq data
            * `layers["counts"]` — a sparse matrix of raw fragment counts
            * `var_names` — a list of region names in 'chr-start-end' format
        * `rna` in `mod` -- h5ad for scRNA-seq data
            * `layers["counts"]` — a sparse matrix of raw UMI counts
    * TODO: fragment files

Your method should be able to load data from this format.

2. Upon successful completion, your model MUST output a NEW `h5mu` that contains some combination of the following
    * uns["TF2r] -- a dataframe of inferred TF-to-region edges with optional weights (THE FORMAT OF THIS IS STILL A WORK IN PROGRESS)
    * uns["r2g"] -- a dataframe of inferred region-to-gene edges with optional weights and p-values
    ```
   source target weight pval (other)
      ch1-1-10     G1    2.5  0.02       X
      ch3-3-27     G2   -1.7  0.05       X
      ch9-8-42     G1    1.2  0.01       X
      ...    ...    ...   ...     ...
    ```
    * uns["grn"] -- a dataframe of inferred TF-to-gene edges with optional weights and p-values
    ```
   source target weight pval (other)
      TF1     G1    2.5  0.02       X
      TF1     G2   -1.7  0.05       X
      TF2     G1    1.2  0.01       X
      ...    ...    ...   ...     ...
    ```
    * uns["tri"] -- a dataframe of inferred TF-gene edges that act through a region with optional weights and p-values. These triplets are often generated using some combination of the above edges. We leave it up to the modeler to decide how to generate these edges.
    ```
   source target   region weight pval (other)
      TF1     G1 ch1-1-10    2.5  0.02       X
      TF1     G2 ch3-3-27   -1.7  0.05       X
      TF2     G1 ch9-8-42    1.2  0.01       X
      ...    ...      ...    ...   ...     ...
    ```
    * uns["tri_filtered"] -- the same as uns["tri"] but with edges that are likely false removed. This is optional, but some methods include a pruning step and it is useful to keep the full set of edges for benchmarking as well as the pruned set.
    ```

    * Your pipeline can also output any other files that are useful for downstream analysis.

Note, the exact combination and whether or not you include both weights and p-values is up to the modeler. Either weights or p-values are required for any outputs however, in order to be able to rank and filter the edges. Please use the same column names as above for weights and p-values to ensure compatibility with downstream benchmarking and analysis.

## Building a snakemake pipeline for your model
If you are building a snakemake pipeline for your model, we recommend following these steps:
1. Copy the template directory (`/src/inference/grn_models/_template`) to a new subdirectory substituting the name of your model (e.g. `src/inference/grn_models/my_grn_model/`).
2. Add any code (e.g. modules, scripts, etc.) you need to run your model to `src/inference/grn_models/my_grn_model/workflow/scripts`
    - Make sure that
3. Build a config file for your model in `src/inference/grn_models/my_grn_model/workflow/config.yaml`.
4. Build a snakemake rule to run your these scripts in `src/inference/grn_models/my_grn_model/workflow/rules`. Uses config 
5. Optional add a
5. Test your pipeline by running it on a small dataset.
        
## Considerations for jamboree
1. Environment is going to be tricky since each method likely needs it's own
2. 









## Writing rules to run your model
1. Write a set of snakemake rule to run your model in `src/inference/grn_models/my_grn_model/workflow/rules`.
2. There is no strict format for these rules, but we recommend following the format of the other GRN models in this directory (e.g. download.smk downloads the appropriate reference files, my_grn_model.smk runs the model).
3. Create a general rule to run your model that includes the above rules (e.g. `src/inference/grn_models/my_grn_model/workflow/general.smk`) (TODO: is this necessary?)
4. Make a Snakefile in `src/inference/grn_models/my_grn_model/workflow` that includes the general rule and any other rules you need to run your model (e.g. `src/inference/grn_models/my_grn_model/workflow/Snakefile`) (see above, is this necessary?)

## Add a environment (TODO)

## Add a config.yaml file (TODO)

## Add a README file

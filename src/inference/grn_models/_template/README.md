author: Adam Klie <br>
email: aklie@ucsd.edu <br>
date: 2024-02-28

# Adding a new GRN model

# Expected input

# Formatting outputs

# 
If your model is simple to run (i.e. can be run with a few functions), we recommend formatting it as described in `src/inference/program_models/_template`. 

However, GRN methods are typically more complex and require multiple steps and scripts (sometimes across programming languages) to run. In this case, we recommend building your model as a snakemake pipeline.

## Setting up a new model directory
1. Copy the template directory (`/src/inference/grn_models/_template`) to a new subdirectory substituting the name of your model (e.g. `src/inference/grn_models/my_grn_model/`).

## Adding code needed to run your model
1. Modify/write your code to take in the expected input [Expected input](#expected-input) and format your outputs as described in [Formatting outputs](#formatting-outputs). See scenic+ and celloracle for examples on how to do this.
2. Add any code (e.g. modules, scripts, etc.) you need to run your model to `src/inference/grn_models/my_grn_model/workflow/scripts` 

## Writing rules to run your model
1. Write a set of snakemake rule to run your model in `src/inference/grn_models/my_grn_model/workflow/rules`.
2. There is no strict format for these rules, but we recommend following the format of the other GRN models in this directory (e.g. download.smk downloads the appropriate reference files, my_grn_model.smk runs the model).
3. Create a general rule to run your model that includes the above rules (e.g. `src/inference/grn_models/my_grn_model/workflow/general.smk`) (TODO: is this necessary?)
4. Make a Snakefile in `src/inference/grn_models/my_grn_model/workflow` that includes the general rule and any other rules you need to run your model (e.g. `src/inference/grn_models/my_grn_model/workflow/Snakefile`) (see above, is this necessary?)

## Add a environment (TODO)

## Add a config.yaml file (TODO)

## Add a README file

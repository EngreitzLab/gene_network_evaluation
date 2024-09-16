# Topyfic Snakemake workflow

This directory contains a Snakemake pipeline for running the Topyfic automatically.

The snakemake will run training (Train) and building model (topModel, Analysis). 

**Note**: Please make sure to install necessary packages and set up your Snakemake appropriately.

**Note**: pipeline is tested for Snakemake >= 8.X ([more info](https://snakemake.readthedocs.io/en/stable/index.html))

## Getting started

### 1. setting up environment

Build your environment and install necessary packages
- [Suggested environment](workflow/envs/Topyfic_env.yml)

### 2. Setting up config file

Modify the [config file](config/config.yaml) or create a new one with the same structure.

1. **names**
   - Contains name of the input dataset(s). 
   - Name will be used as a name of train and topModel models
   - If there is multiple names, Topyfic will normalize the models across names using [harmony](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884693/).
   - list of name: `[parse, 10x]`

2. **count_data**
   - Contains path of each input data
   - Name of each path should match name in `names`
   - Recommended to use full path rather than relative path

3. **n_topics**
   - Contains list of number of initial topics you wish to train model base on them
   - list of int: `[5, 10, 15, 20, 25, 30, 35, 40, 45, 50]`

4. **organism**
   - Indicate spices which will be used for downstream analysis
   - Example: human or mouse

5. **workdir**
   - Directory to put the outputs
   - Make sure to have write access.
   - It will create one folder per dataset.

6. **train**
   - most of the item is an input of `train_model()`
   - n_runs: number of run to define rLDA model (default: 100)
   - random_states: list of random state, we used to run LDA models (default: range(n_runs))

7. **top_model**
   - n_top_genes (int): Number of highly-variable genes to keep (default: 50)
   - resolution (int): A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. (default: 1)
   - max_iter_harmony (int): Number of iteration for running harmony (default: 10)
   - min_cell_participation (float): Minimum cell participation across for each topic to keep them, when is `None`, it will keep topics with cell participation more than 1% of #cells (#cells / 100)

8. **merge**
   - Indicate if you want to also get a model for all data together.


### 3. Run snakemake

First run it with `-n` to make sure the steps that it plans to run are reasonable. 
After it finishes, run the same command without the `-n` option.

`snakemake -n`

For SLURM:

```
snakemake \
-j 1000 \
--latency-wait 300 \
--use-conda \
--rerun-triggers mtime \
--executor cluster-generic \
--cluster-generic-submit-cmd \
"sbatch -A model-ad_lab \
  --partition=highmem \
  --cpus-per-task 16 \
  --mail-user=nargesr@uci.edu \
  --mail-type=START,END,FAIL \
  --time=72:00:00" \
-n \
-p \
--verbose
```
highmem
standard

Development hints: If you ran to any error `-p --verbose` would give you more detail about each run and will help you to debug your code.

Once you get all the three main objects (Train, TopModel, Analysis), we could move to the evaluation pipeline.




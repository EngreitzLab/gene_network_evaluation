# `evaluation`

This directory contains the code and outputs for running gene program evaluation on the objects contained in the [`inference` directory](../inference/). Subdirectories are first broken up by dataset, and then by `prog_key`. 

Each of these program_key subdirectories contains the outputs of the full evaluation pipeline. See the [main README.md](../../README.md) for more information on the input and output formats for evaluation.

## Configuration files
Running the evaluation pipeline can be done in several ways make use of a configuration file. `evaluation_pipeline.yaml` contains an example evaluation pipeline configuration file. This file can be copied and modified to specify the evaluation pipeline for a given dataset and set of gene programs. It has 7 sections

1. `io`: Contains the paths to the input and output files for the evaluation pipeline.
2. `categorical_association`: Contains the parameters for the categorical association analysis.
3. `perturbation_association`: Contains the parameters for the perturbation association analysis.
4. `gene_set_enrichment`: Contains the parameters for the gene set enrichment analysis.
5. `trait_enrichment`: Contains the parameters for the trait enrichment analysis.
6. `motif_enrichment`: Contains the parameters for the motif enrichment analysis.
7. `explained_variance`: Contains the parameters for the explained variance analysis.

For an example of how to generate multiple configuration files for different program runs on the same dataset, see the [`examples/evaluation/iPSC_EC/make_configs.ipynb`](/examples/evaluation/iPSC_EC/make_configs.ipynb) notebook.

This configuration file can be used as input to the `evaluation_pipeline.py` script as such:

```bash
python evaluation_pipeline.py --config evaluation_pipeline.yaml
```

Or can be loaded into the `evaluation_pipeline.ipynb` notebook for interactive use.

We have included an example SLURM script for running the evaluation pipeline in parallel on a cluster: `evaluation_pipelines.sh`

Example configuration files can be found at:
-  [`examples/evaluation/Endothelial/cNMF/evaluation_pipeline.yml`](/examples/evaluation/Endothelial/cNMF/evaluation_pipeline.yml)
-  [`examples/evaluation/iPSC_EC/cNMF/cNMF_30/evaluation_pipeline.yml`](/examples/evaluation/iPSC_EC/cNMF/cNMF_30/evaluation_pipeline.yml)
## Selection of k
Most program inference methods require the prespecification of the number of gene programs prior to inference. This is a key hyperparameter that can have a large impact on the quality of the inferred programs. We currently do not have an automated method for selecting k, but provide several example notebooks for each method in the [`examples/evaluation/k_selection`](/examples/evaluation/k_selection) directory. These notebooks can be used to explore the effect of different k values on the quality of the inferred programs.

# Note:
This pipeline will eventually be integrated into a snakemake workflow for easier reproducibility. This is a known issue and will be addressed in future releases. This represents a stopgap solution for the purposes of the jamboree.

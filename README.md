# gene_program_evaluation
* Evaluation framework for gene program inferred from single-cell omics data. 
* Cell x program scores are expected in the anndata format with metadata stored according to specification. 
* single-cell omics data and program scores are supplied jointly in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).

![image](https://github.com/EngreitzLab/gene_program_evaluation/assets/25486108/946dd101-97de-485c-9eab-baffb0f7fd8a)

## Roadmap
#### 1. Reimplement cNMF pipeline evaluations
* ~~Variance explained~~
* ~~Batch association~~
* ~~GO Term/gene-set~~
* ~~motif enrichment~~
* Perturbation significance (V2G2P paper)
#### 2. Implement snakemake pipeline
* Rules for evaluations
* Plotting functions
* Report generation
* Containers for rule dependencies
#### 3. Implement evaluation [ideas from sub-group](https://docs.google.com/spreadsheets/d/15a9xLCvqBuh5mUtXj8hq6JD55qPCIf5b6cgCYZKZDUI/edit#gid=1041024840)
* Cell-type specificity
* Cross-modality prediction
* Indirect perturbation effect sensitivity
#### 4. Maintenance & robustness routines
* Finish functional & unit tests
* Continuous integration
* Code linting
      

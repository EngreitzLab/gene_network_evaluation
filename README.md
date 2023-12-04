# gene_program_benchmarking
* Evaluation framework for gene program inferred from single-cell omics data. 
* Cell x program scores are expected in the anndata format with metadata stored according to specification. 
* single-cell omics data and program scores are supplied jointly in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).

![image](https://github.com/EngreitzLab/gene_program_benchmarking/assets/25486108/6d194ca5-ceba-40e2-89db-f67115120c01)
  
## Roadmap
1. Reimplement cNMF pipeline benchmarks
    * ~~Variance explained~~
    * ~~Batch association~~
    * GO Term/gene-set/motif enrichment (GSEA, ORA...)
    * Perturbation significance (V2G2P paper)
2. Implement benchmark [ideas from sub-group](https://docs.google.com/spreadsheets/d/15a9xLCvqBuh5mUtXj8hq6JD55qPCIf5b6cgCYZKZDUI/edit#gid=1041024840)
    * Cell-type specifcity
    * Cross-modality prediction
    * Indirect perturbation effect sensitivity
      

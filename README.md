# gene_program_benchmarking
* Evaluation framework for gene program inferred from single-cell omics data. 
* Cell x program scores are expected in the anndata format with metadata stored according to specification. 
* single-cell omics data and program scores are supplied jointly in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).

![Screenshot 2023-11-29 153708](https://github.com/EngreitzLab/gene_program_benchmarking/assets/25486108/60384810-3ea8-4aba-b73d-0fde171bd595)
  
## Roadmap
1. Reimplement cNMF pipeline benchmarks
    * Variance explained
    * ~~Batch association~~
    * GO Term/gene-set/motif enrichment (GSEA, ORA...)
    * Perturbation significance (V2G2P paper)
2. Implement benchmark ideas from sub-group
    * Cell-type specifcity
    * Cross-modality prediction
    * Indirect perturbation effect sensitivity
      

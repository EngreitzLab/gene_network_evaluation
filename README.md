# Modelling gene programs as latent variables
Gene program inference methods are "factor analysis" style latent variable models that decompose single-cell data into a $cell \times program$ and a $program \times feature$ matrices. The models can generally be described with the following expression. For single-cell data $\textbf{X} \in \mathbb{R}^{c \times f}$ with $c$ cells and $f$ features,

$$\Large \textbf{X} = P \cdot W + \epsilon$$

where $\textbf{P} \in \mathbb{R}^{c \times k}$ and $\textbf{W} \in \mathbb{R}^{k \times f}$ are the cell-wise program scores and feature weights per program respectively for $k$ programs. The programs represent the coordinated activity of genes towards a biological function. These models rely on the correlation structure of the data to estimate these matrices. If high spurious correlation is present in the data, the models can report latent components that do not reflect biological activity. Furthermore, while the models are purely correlative (not causal), training on interventional data implicitly benefits the model due to the alterations in the correlation structure of the data.

A few examples of such models are [cNMF](https://github.com/dylkot/cNMF), [LDVAE](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html), [f-scLVM](https://github.com/scfurl/f-scLVM). More advanced models consider additional inputs or work on multi-omic data such as [Spectra](https://github.com/dpeerlab/spectra/), [muVI](https://github.com/MLO-lab/MuVI). 

### <ins>Wrap methods of choice with inference/methods/_template.py and create a pull request</ins> 

# Evaluating the biological basis of learnt programs
The goal is to develop evaluation criteria covering various technical and biological facets. This process includes both coming up with a concept and then an implementation. Following is a list of implemented or planned criteria. 

### <ins>Implement evaluations using evaluations/methods/_template.py and create a pull request.</ins>

| Criterion    | Implementation | External resource | Interpretation | Caveats |
| -------- | ------- | -------- | ------- | ------- |
| Goodness of fit  | [Explained variance](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html) per program    | None | A program explaining more variance in the data might represent dominant biological variation. | Technical variation might be the highest source of variance (e.g. batch effects). |
| Variation across batch levels | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are high in a batch then the component likely is modelling technical noise. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Variation across cell-types | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are high in a cell-type then the component likely is modelling a process specific to the cell-type. | Rare cell-types may not contribute enough co-variation of genes to be learned effectively. |
| Variation across conditions | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are high in a condition then the component likely is modelling a process specific to the condition. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Gene-set enrichment | [GSEA](https://gseapy.readthedocs.io/en/latest/introduction.html) using program x feature scores | MsigDB, Enrichr | If a program is significantly associated with a gene-set then it could explain the biological process the program represents | |
| Motif enrichment | Pearson correlation of motif counts per gene and program x gene scores | HOCOMOCO v12 | If genes with high contributions to a program are also enriched with same enhancer/promoter motifs they could be co-regulated | A biological pathway could involve genes with different regulation but still contribute to a common function | 
| Perturbation sensitivity |  | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 
| Cross-modality prediction |  | Multi-omic data | If program x feature scores learnt from one modality can lead to a good fit for another modality/dataset indicate a robust biological connection b/w genes in the program | Technical variation b/w datasets and mapping features b/w modalities are practical challenges | 

# Structure of the evaluation pipeline
* Cell x program scores are expected in the anndata format with metadata stored according to specification. 
* single-cell omics data and program scores are supplied jointly in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).
* 
## Roadmap
#### 1. Reimplement cNMF pipeline evaluations
* ~~Variance explained~~
* ~~Batch association~~
* ~~GO Term/gene-set~~
* ~~motif enrichment~~
* Perturbation significance (V2G2P paper)
* Add doc-strings
#### 2. Implement snakemake pipeline
* Rules for evaluations
* Plotting functions
* Report generation
* Containers for rule dependencies
#### 4. Maintenance & robustness routines
* Finish functional & unit tests
* Continuous integration
* Code linting
      

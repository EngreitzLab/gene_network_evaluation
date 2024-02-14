## Assessing the biological relevance of data-driven gene networks

Gene networks inferred from single-cell genomic data (scRNASeq., scATACseq., multi-omics and Perturb-seq.) are useful in discovering contextual biological mechanisms. These networks can be viewed as data-driven hypotheses of gene interactions. We aim to implement a flexible framework to evaluate the plausibility of networks inferred by computational methods. The assessment is broken down into themes such as goodness if fit (ability to explain the data), co-regulation, mechanistic interactions etc. Under each theme, multiple evaluation tasks are conceptualised and implemented using appropriate statistical tests.

Gene network inference methods are further classified as gene program inference, gene regulatory network (GRN) inference and enchancer-gene (E2G) linking methods. An overview of each method class is given below.

### Modelling gene programs as latent variables
Gene program inference methods are "factor analysis" style latent variable models that decompose single-cell data into a $cell \times program$ and a $program \times feature$ matrices. The models can generally be described with the following expression. For single-cell data $\textbf{X} \in \mathbb{R}^{c \times f}$ with $c$ cells and $f$ features,

$$\Large \textbf{X} = P \cdot W + \epsilon$$

where $\textbf{P} \in \mathbb{R}^{c \times k}$ and $\textbf{W} \in \mathbb{R}^{k \times f}$ are the cell-wise program scores and feature weights per program respectively for $k$ programs. The programs represent the coordinated activity of genes towards a biological function. These models rely on the correlation structure of the data to estimate these matrices. If high spurious correlation is present in the data, the models can report latent components that do not reflect biological activity. Furthermore, while the models are purely correlative (not causal), training on interventional data implicitly benefits the model due to the alterations in the correlation structure of the data.

A few examples of such models are [cNMF](https://github.com/dylkot/cNMF), [LDVAE](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html), [f-scLVM](https://github.com/scfurl/f-scLVM). More advanced models consider additional inputs or work on multi-omic data such as [Spectra](https://github.com/dpeerlab/spectra/), [muVI](https://github.com/MLO-lab/MuVI). 

### Gene regulatory network inference

Gene regulatory networks are tri-partite graphs connecting transcription factors (TFs) - candidate regulatory elements (CREs) - genes. Such methods aim to model the regulatory process governing gene expression and are typically trained using expression and genomic data together. Examples of such methods are [SCENIC+](https://doi.org/10.1038/s41592-023-01938-4).

### Enhancer-gene linking

Enhancer-gene linking methods aim to identify and link enhancer sequences to genes. The methods require genomic data such as ATACseq. but may integrate expression data as well. The methods aim to quantify the regulatory impact of enhancer sequences on downstream genes (and gene expression). Examples of such methods are [scE2G](https://github.com/EngreitzLab/sc-E2G), [ABC](https://www.nature.com/articles/s41588-019-0538-0) and [Enhlink](https://www.biorxiv.org/content/10.1101/2023.05.11.540453v1).

## Structure of the evaluation pipeline
* Single-cell omics data and outputs from computational inference are stored in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).
* Templates for implementing evaluations or method-wrappers can be found in ```src/evaluation/``` and ```src/inference/``` respectively.
* Evaluations and methods implemented under ```src/``` are stitched together in an evaluation pipeline ```smk/```
  
![image](https://github.com/EngreitzLab/gene_program_evaluation/assets/25486108/eb14e159-de07-47bd-87d4-12d5b94102fa)

### Evaluating the biological basis of inferred gene programs

| Criterion    | Implementation | External resource | Interpretation | Caveats |
| -------- | ------- | -------- | ------- | ------- |
| Goodness of fit  | [Explained variance](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html) per program    | None | A program explaining more variance in the data might represent dominant biological variation. | Technical variation might be the highest source of variance (e.g. batch effects). |
| Variation across category levels | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are variable between batch levels then the component likely is modelling technical noise. Alternatively, if program scores are variable between a biological category like cell-type or condition then the program is likely modelling a biological process specific to the category. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Gene-set enrichment | [GSEA](https://gseapy.readthedocs.io/en/latest/introduction.html) using program x feature scores | MsigDB, Enrichr | If a program is significantly associated with a gene-set then it could explain the biological process the program represents | |
| Motif enrichment | Pearson correlation of motif counts per gene (promoter or enchancer) and program x gene scores | HOCOMOCO v12 | If genes with high contributions to a program are also enriched with same enhancer/promoter motifs they could be co-regulated | A biological pathway could involve genes with different regulation but still contribute to a common function | 
| Co-regulation |  | TF-gene links | If a program contains genes that are regulated by TFs that form a regulatory module then it indicates mechansitic commonality  |  | 
| Perturbation sensitivity |  | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 
| Cross-modality prediction |  | Multi-omic data | If program x feature scores learnt from one modality can lead to a good fit for another modality/dataset indicate a robust biological connection b/w genes in the program | Technical variation b/w datasets and mapping features b/w modalities are practical challenges | 

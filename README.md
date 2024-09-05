# Assessing the biological relevance of data-driven gene programs

Gene programs inferred from single-cell genomic data (scRNASeq., scATACseq., multi-omics and Perturb-seq.) are useful in discovering contextual biological mechanisms. These programs can be viewed as data-driven hypotheses of gene interactions. We aim to implement a flexible framework to evaluate the plausibility of programs inferred by computational methods. The assessment is broken down into themes such as goodness if fit (ability to explain the data), co-regulation, mechanistic interactions etc. Under each theme, multiple evaluation tasks are conceptualised and implemented using appropriate statistical tests.

## Structure of the pipeline
![image](https://github.com/EngreitzLab/gene_program_evaluation/assets/25486108/eb14e159-de07-47bd-87d4-12d5b94102fa)

Gene program inference and analysis can be broken down into 3 major steps:
1. [**Modelling gene programs as latent variables**](#modelling-gene-programs-as-latent-variables)
2. [**Evaluating the biological basis of inferred gene programs**](#evaluating-the-biological-basis-of-inferred-gene-programs)
3. [**Analysis and annotation of gene programs**](#analysis-and-annotation-of-gene-programs)

The following sections describe the key components of each step.

### Modelling gene programs as latent variables

#### Quick start:
* Several gene program inference methods are implemented in [`src/inference/`](src/inference/).
* Examples showing how to run inference using these methods can be found in [`examples/inference`](examples/inference).
* Templates for implementing inference method-wrappers can be found in [`src/inference/_template`](src/inference/_template).

#### Expected inputs and outputs of inference methods:
Running the inference methods in this repository requires that the original single-cell data is stored in the MuData format. The format is described in detail in the [MuData documentation](https://mudata.readthedocs.io/en/latest/). The requirments for this MuData object are minimal:

1. `mdata.mod[data_key]` - The original single-cell data should be stored in an AnnData object under a user specified key (e.g. `rna`).
    - The AnnData should contain a n_obs x n_vars `.X` matrix of single-cell data where n_obs is the number of cells and n_vars is the number of genes/peaks/features.
    - The AnnData can also contain other layers with transformed data (e.g. normalized counts, log-transformed counts, etc.) in `.layers`. This is recommended for most inference methods.
    - Some additional inputs are required for running specific evaluations (e.g. guide assignments for Perturb-seq evaluations, enhancer-to-gene links for motif enrichment analysis). These are covered in the next section.

Running the evaluations on gene programs in this repository requires that the inferred gene programs be added to *this same* MuData object. We require the following to be present:

1. `mdata.mod[data_key]` - The original single-cell data should be stored in an AnnData object under a user specified key (e.g. 'rna').
    - This AnnData should be the unmodified object from above that was used to run the inference methods.
2. `mdata.mod[program_key]` - The inferred gene programs should be stored in an AnnData object under a user specified key (e.g. 'cNMF').
    - The AnnData should contain a n_obs x n_vars `.X` matrix where n_obs is the number of cells and n_vars is the number of gene programs inferred. These values are the program scores/loadings for each cell and is sometimes called a cell membership matrix.
    - The `.varm` of this AnnData should contain a key called `loadings` which is NumPy array of shape n_vars x n_features where n_features is the number of genes/peaks/features used for inference. This matrix contains the weights of each gene/peak/feature in each program and is sometimes called a gene loadings matrix.
    The `.uns` of this AnnData should contain a key called `var_names` which is a 1D NumPy array of containing the names of the genes/peaks/features used for inference.
    - This object *does not* need to contain any metadata about the cells in `.obs` (though it can).
    
For examples inference methods inputs, see [`src/examples/datasets`](src/examples).
For examples on how to run inference methods and method outputs, see [`src/inference/inference`](src/examples/inference).

#### Some background on gene program inference methods:
Gene program inference methods are "factor analysis" style latent variable models that decompose single-cell data into a $cell \times program$ and a $program \times feature$ matrices. The models can generally be described with the following expression. For single-cell data $\textbf{X} \in \mathbb{R}^{c \times f}$ with $c$ cells and $f$ features,

$$\Large \textbf{X} = P \cdot W + \epsilon$$

where $\textbf{P} \in \mathbb{R}^{c \times k}$ and $\textbf{W} \in \mathbb{R}^{k \times f}$ are the cell-wise program scores and feature weights per program respectively for $k$ programs. The programs represent the coordinated activity of genes towards a biological function. These models rely on the correlation structure of the data to estimate these matrices. If high spurious correlation is present in the data, the models can report latent components that do not reflect biological activity. Furthermore, while the models are purely correlative (not causal), training on interventional data implicitly benefits the model due to the alterations in the correlation structure of the data.

A few examples of such models are [cNMF](https://github.com/dylkot/cNMF), [Topyfic](https://github.com/mortazavilab/Topyfic), [LDVAE](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html), [f-scLVM](https://github.com/scfurl/f-scLVM). More advanced models consider additional inputs or work on multi-omic data such as [Spectra](https://github.com/dpeerlab/spectra/), [muVI](https://github.com/MLO-lab/MuVI). 

### Evaluating the biological basis of inferred gene programs

#### Quick start:
* Several gene program evaluation methods are implemented in [`src/evaluation/`](src/evaluation/).
* Templates for implementing evaluations can be found in [`src/evaluation/_template`](src/evaluation/_template).
* Examples showing how to run evaluations using these methods can be found in [`examples/evaluation`](examples/evaluation).

#### Evaluation criteria:
| Criterion    | Implementation | External resource | Interpretation | Caveats |
| -------- | ------- | -------- | ------- | ------- |
| Goodness of fit  | [Explained variance](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html) per program | None | A program explaining more variance in the data might represent dominant biological variation. | Technical variation might be the highest source of variance (e.g. batch effects). |
| Variation across category levels | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are variable between batch levels then the component likely is modelling technical noise. Alternatively, if program scores are variable between a biological category like cell-type or condition then the program is likely modelling a biological process specific to the category. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Gene-set enrichment | [GSEA](https://gseapy.readthedocs.io/en/latest/introduction.html) using program x feature scores | MsigDB, Enrichr | If a program is significantly associated with a gene-set then it could explain the biological process the program represents | |
| Motif enrichment | Pearson correlation of motif counts per gene (promoter or enchancer) and program x gene scores | HOCOMOCO v12 | If genes with high contributions to a program are also enriched with same enhancer/promoter motifs they could be co-regulated | A biological pathway could involve genes with different regulation but still contribute to a common function | 
| Perturbation sensitivity |  | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 

#### Expected inputs and outputs of evaluation methods:
The expected inputs to the evaluation pipeline follow the same format as the outputs to the [inference methods](#expected-inputs-and-outputs-of-inference-methods). For certain evaluation methods, however, additional inputs are required:

1. Variation across category levels - This evaluation requires that the cell metadata be stored in the MuData object under a user specified column in `mdata.mod[data_key].obs`. This evaluation can also be run at the pseudobulk level by providing a separate column name that specifies the pseudobulk grouping.
2. Motif enrichment - This evaluation requires that separate enhancer-to-gene links be provided. Enhancer-gene linking methods aim to quantify the regulatory impact of enhancer sequences on downstream genes (and gene expression). Examples of such methods are [scE2G](https://github.com/EngreitzLab/sc-E2G), [ABC](https://www.nature.com/articles/s41588-019-0538-0) and [Enhlink](https://www.biorxiv.org/content/10.1101/2023.05.11.540453v1). Enhancer-to-gene links should be stored in the following format:

```plaintext
chromosome      start   end     seq_name        seq_class       seq_score       gene_name
chr12   52299565        52301301        promoter|chr12:52299415-52301451        promoter        56.971568       ACVRL1
chr12   52299565        52301301        promoter|chr12:52299415-52301451        promoter        56.971568       ANKRD33
chr12   57028808        57030874        promoter|chr12:57028658-57031024        promoter        41.6775 BAZ2A
chr12   57028808        57030874        promoter|chr12:57028658-57031024        promoter        41.6775 ATP5B
```
3. Perturbation sensitivity - This evaluation requires that the guide assignments for each cell be stored in the MuData object as follows:
    - `.uns["guide_names"]` - A 1D NumPy array containing the names of the guides used in the perturbation experiment.
    - `.uns["guide_targets"]` - A 1D NumPy array containing the gene targets of the guides used in the perturbation experiment
    - `.obsm["guide_assignment"]` - A 2D NumPy array of shape n_obs x n_guides containing the guide assignments for each cell. The values should be 0 or 1 indicating whether the cell was assigned to the guide or not.

For examples of evaluation method inputs see [`src/examples/evaluation`](src/examples/inference).
For examples of evaluation method outputs and how to run evaluation methods see [`src/evaluation/evaluation`](src/examples/evaluation).


#### Selection of the best k:
TODO

#### Details on evaluation methods:

##### Explained variance
TODO

##### Variation across category levels


##### Gene set enrichment
TODO

##### Motif enrichment
TODO

##### Perturbation sensitivity
TODO

### Analysis and annotation of gene programs

#### Quick start:
Once you have run the evaluation pipeline on your set of gene programs, you can visualize and annotate your programs with our dashboard:
```bash
python app/app.py --config /path/to/config.yaml
```


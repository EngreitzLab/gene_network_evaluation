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
    - The AnnData can also contain other layers with transformed data (e.g. normalized counts, log-transformed counts, etc.) in `.layers`. This is recommended for use with most inference methods.
    - Some additional inputs are required for running specific evaluations (e.g. guide assignments for Perturb-seq evaluations, enhancer-to-gene links for motif enrichment analysis). These are covered in the [evaluation section](#expected-inputs-and-outputs-of-evaluation-methods).

Running the evaluations on gene programs in this repository requires that the inferred gene programs be added to *this same* MuData object. We require the following to be present:

1. `mdata.mod[data_key]` - The original single-cell data should be stored in an AnnData object under a user specified key (e.g. 'rna').
    - This AnnData should be the unmodified object from above that was used to run the inference method.
2. `mdata.mod[program_key]` - The inferred gene programs should be stored in an AnnData object under a user specified key (e.g. 'cNMF').
    - The AnnData should contain a n_obs x n_vars `.X` matrix where n_obs is the number of cells and n_vars is the number of programs inferred. We call these values the program scores for each cell (sometimes called a cell membership matrix).
    - The `.varm` of this AnnData should contain a key called `loadings` which is NumPy array of shape n_vars x n_features where n_features is the number of genes/peaks/features used for inference. This matrix contains the weights of each gene/peak/feature in each program.
    - The `.uns` of this AnnData should contain a key called `var_names` which is a 1D NumPy array containing the names of the genes/peaks/features used for inference.
    - This object *does not* need to contain any metadata about the cells in `.obs` (though it can).
    
For examples of inference methods inputs, see [`examples/datasets`](examples).
For examples on how to run inference methods and method outputs, see [`examples/inference`](examples/inference).

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
| Trait enrichment | [Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) | OpenTargets database | If a program is significantly associated with a trait then it could explain the biological process the program represents | |
| Perturbation sensitivity | Mann-Whitney U test of program scores between perturbed cells and non-targeted/reference cells | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 

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
    - `mdata[data_key].uns["guide_names"]` - A 1D NumPy array containing the names of the guides used in the perturbation experiment.
    - `mdata[data_key].uns["guide_targets"]` - A 1D NumPy array containing the gene targets of the guides used in the perturbation experiment
    - `mdata[data_key].obsm["guide_assignment"]` - A 2D NumPy array of shape n_obs x n_guides containing the guide assignments for each cell. The values should be 0 or 1 indicating whether the cell was assigned to the guide or not.

For examples of evaluation method inputs see [`examples/evaluation`](examples/inference).
For examples of evaluation method outputs and how to run evaluation methods see [`evaluation/evaluation`](src/examples/evaluation).

#### Selection of the best k:
A crucial step in program analysis is the selection of the number of programs (k) to infer. This is a hyperparameter for almost all methods, and is highly dependent on dataset (size, complexity, etc.). Many methods have their own criteria for selecting k, and **we are in the process of determining the best practices for selecting k across methods**. Worked examples (in notebooks) of how we select k for `cNMF` and `Topyfic` can be found in the `examples/evaluation/k_selection` directory.

#### Details on evaluation methods (including expected output file formats):
Each evaluation method is implemented as a separate Python file in the `src/evaluation/` directory. The ideas behind the implementation of each evaluationm as well as expected output file formats, are described below:

##### Explained variance
The explained variance per program is the proportion of the variance in the original data that is modeled by that program. We calculate this by taking the outer product between each program's score vector and its feature loadings vector. This gives us a cell x feature matrix that can be thought of as a reconstruction of the original data using only that program. Calculating a ratio of the variance in this reconstruction to the variance in the original data gives us the explained variance. This can give a sense of the program's 'importance' in the overall model, however, technical variation might be the highest source of this variance (e.g., batch effects)

We use the `sklearn.metrics.explained_variance_score` function to calculate the explained variance ratio.

Output file format that is compatible with the dashboard:
```plaintext
program_name    explained_variance_ratio
program_1       0.5
program_2       0.3
program_3       0.2
```
This should saved in a file called: `{prog_key}_variance_explained_ratio.txt` where `prog_key` is the key used to store the programs in the MuData object.

##### Variation across category levels
We have implemented several tests to associate programs with a categorical covariate. For example, if we have a program that represents a biological process that is specific to a cell type, we would expect the program scores to be significantly different between cell types. Our goals here are to 1) determine the overall association of each program with the categorical coveriate and 2) identify which levels of the categorical covariate are driving the association.

By default, we binarize the categorical covariate for each level (e.g. 1 if categorical_key == level and 0 otherwise) and then perform a Pearson correlation between the program scores and the binarized covariate. We then take the maximum correlation value as the test statistic for the program's overall association with the covariate.

We have also implemented a Kruskal-Wallis test paired with a Dunn's posthoc test. The Kruskal-Wallis test is a non-parametric ANOVA that tests the null hypothesis that the population median of all groups are equal. The Kruskal-Wallis test only tells us if there is a significant difference between groups, but does not tell us which groups are different from each other. By defualt, we use a Dunn's posthoc test to identify which groups are significantly different from each other, although other posthoc methods are available.

This evaluation outputs two files:

1. `{prog_key}_{categorical_key}_association_results.txt` - This file contains each program's overall association with the categorical covariate.
```plaintext
program_name    {categroical_key}_{test}_stat	{categroical_key}_{test}_pval
program_1       10.0    0.001
program_2       5.0     0.01
program_3       2.0     0.1
```

`e.g. cNMF_batch_association_results.txt`

2. `{prog_key}_{categorical_key}_association_posthoc.txt` - This file contains the individual level's association with the program.
```plaintext
program_name    {categroical_key}_{level_1}_stat    {categroical_key}_{level_1}_pval    {categroical_key}_{level_1}_adj_pval    {categroical_key}_{level_2}_log2FC    {categroical_key}_{level_2}_pval    {categroical_key}_{level_2}_adj_pval ...
program_1       10.0    0.001   0.001   2.0     0.01    0.01
program_2       5.0     0.01    0.01    1.0     0.1     0.1
program_3       2.0     0.1     0.1     0.5     0.5     0.5
```

e.g. `cNMF_batch_association_posthoc.txt`

##### Gene set enrichment
To tease apart the biological processes that the program represents, we look for enrichment of each program with published gene sets. Here we use two primary methods for gene set enrichment analysis: GSEA and a Fisher's exact test.

In GSEA, we rank the genes in the program by their program loadings and then test for enrichment of gene sets in the ranked list. We use the `gseapy` package to run GSEA.

For the Fisher's exact test, we use a contingency table to test for enrichment of gene sets in the program. We use `gseapy` to run the Fisher's exact test.

Output file format that is compatible with the dashboard:
```plaintext
program_name	term	pval	adj_pval	enrichment
program_1	term_1	0.001	0.001	1.5
program_1	term_2	0.01	0.01	1.2
program_1	term_3	0.1	0.1	0.8
```
This should saved in a file called: `{prog_key}_{library}_{test}_enrichment.txt` where `prog_key` is the key used to store the programs in the MuData object, `library` is the gene set library used for the enrichment analysis (e.g. "Reactome_2022"), and `test` is the method used for the enrichment analysis (e.g. "gsea" or "fisher").

e.g. `cNMF_Reactome_2022_fisher_geneset_enrichment.txt`

##### Motif enrichment
To identify potential regulators of each program, we look for enrichment of motifs in the enhancers or promoters of genes that are highly weighted in a program. Specifically, we scan the promoters and enhancers of genes of all genes in a program for the presence of motifs from a user provided database using the FIMO algorithm. We then aggregate the number of significant motif hits for each gene and correlate this counts vector with the program scores (by default with Pearson correlation), computing a p-value for the correlation. At both stages (motif scanning and correlation calcualtion) p-values are corrected for multiple testing using the Benjamini-Hochberg procedure.

We use the `tangermeme` implementation of FIMO, which is very fast when running on a GPU, hence we recommend running this evaluation on a GPU when possible.

Output file format that is compatible with the dashboard:
```plaintext
program_name	motif	stat	pval	adj_pval
program_1	motif_1	0.5	0.001	0.001
program_1	motif_2	0.3	0.01	0.01
program_1	motif_3	0.2	0.1	0.1
```
This should saved in a file called: `{prog_key}_{E/P_type}_{motif_database}_{test}_motif_enrichment.txt` where `prog_key` is the key used to store the programs in the MuData object, `E/P_type` is the type of region type scanned for motifs (e.g. "enhancer" or "promoter"), `motif_database` is the motif database used for the enrichment analysis (e.g. "hocomoco_v12"), and `test` is the method used for the enrichment analysis (e.g. "ttest").

e.g. `cNMF_promoter_hocomoco_ttest_motif_enrichment.txt`

It is often more XX to run this evaluation in a stratified manner across a categroical covariate of interest (e.g. cell type, batch, day of differentation). In this case, we expect a file for each level of the categorical covariate. The file should be named: `{prog_key}_{E/P_type}_{motif_database}_{test}_{stratification_key}_{level_key}_motif_enrichment.txt`.

e.g. `cNMF_30_enhancer_test_pearsonr_sample_D0_motif_enrichment.txt`

##### Trait enrichment
To identify potential associations of each program with traits, we look for enrichment of traits in the genes that are highly weighted in a program. Specifically, we query the OpenTargets database for genes associated with traits. By default, we use a locus2gene (L2G) threshold of 0.2 to determine the set of genes associated with each trait. We then perform a Fisher's exact test to test for enrichment of traits in the program. We use the Benjamini-Hochberg procedure to correct for multiple testing.

Output file format that is compatible with the dashboard:
```plaintext
term	adj_pval	trait_efos	trait_category	program_name	enrichment	trait_reported	genes	study_id	pmid	-log10(adj_pval)
term_1	0.001	term_1	term_category_1	program_1	1.5	term_reported_1	genes_1	study_id_1	pmid_1	3
term_2	0.01	term_2	term_category_2	program_1	1.2	term_reported_2	genes_2	study_id_2	pmid_2	2
term_3	0.1	term_3	term_category_3	program_1	0.8	term_reported_3	genes_3	study_id_3	pmid_3	1
```

This should saved in the similar format to geneset enrichment: `{prog_key}_{trait_database}_{test}_trait_enrichment.txt` where `prog_key` is the key used to store the programs in the MuData object, `trait_database` is the trait database used for the enrichment analysis (e.g. "OpenTargets"), and `test` is the method used for the enrichment analysis (e.g. "fisher").

e.g. `cNMF_OT_GWAS_fisher_trait_enrichment.txt`

##### Perturbation sensitivity
To test the sensitivity of a program to perturbations, we compare the program scores of cells that were perturbed with cells that were not perturbed. We use a Mann-Whitney U test to test for differences in program scores between the two groups. We use the Benjamini-Hochberg procedure to correct for multiple testing. If a program is sensitive to perturbations, we would expect the program scores to be significantly different between the two groups.

Output file format that is compatible with the dashboard:
```plaintext
target_name	program_name	stat	log2FC	pval	adj_pval
target_1	program_1	10.0	1.5	0.001	0.001
target_2	program_1	5.0	1.2	0.01	0.01
target_3	program_1	2.0	0.8	0.1	0.1
```

This should saved in a file called: `{prog_key}_{gene_guide}_perturbation_association.txt` where `prog_key` is the key used to store the programs in the MuData object and`gene_guide` indicates the resolution at which the perturbations where tested (e.g. "gene" or "guide").

e.g. `cNMF_gene_perturbation_association.txt`

Like motif enrichment, this is another analysis that we often want to run in a stratified manner across a categroical covariate of interest (e.g. cell type, batch, day of differentation). In this case, we expect a file for each level of the categorical covariate. The file should be named: `{prog_key}_{gene_guide}_{stratification_key}_{level_key}_perturbation_association.txt`.

e.g. `cNMF_30_gene_sample_D0_perturbation_association.txt`

### Analysis and annotation of gene programs

#### Quick start:
Once you have run the evaluation pipeline on your set of gene programs, you can visualize and annotate your programs with our dashboard:
```bash
python app/app.py --config /path/to/config.yaml
```

Minimally, the configuration file should contain the following information:
```yaml
path_evaluation_outs:  # List of paths to the evaluation output directories
  - "/path/to/evaluation/output/directory_1"
  - "/path/to/evaluation/output/directory_2"
path_evaluation_config: "/path/to/evaluation/config/evaluation.yaml"  # Path to the evaluation configuration file
path_mdata: "/path/to/eval.h5mu"  # Path to the MuData object containing the data_key and program_key
path_report_out: "/path/to/report/output/directory"  # Path to the directory where the report will be saved
data_key: "rna"  # Key used to store the original single-cell data in the MuData object
categorical_keys:   # List of categorical keys to use for dashboard
  - "sample"
continuous_keys:  # List of continuous keys to use for dashboard
  - "n_counts"
annotations_loc: "annotations.csv"  # Name of file in the report directory to dump annotations
```

See the [app README](app/README.md) for more information on how to generate the dashboard.


# Modelling gene programs as latent variables
Gene program inference methods are "factor analysis" style latent variable models that decompose single-cell data into a $cell \times program$ and a $program \times feature$ matrices. The models can generally be described with the following expression. For single-cell data $\textbf{X} \in \mathbb{R}^{c \times f}$ with $c$ cells and $f$ features,

$$\Large \textbf{X} = P \cdot W + \epsilon$$

where $\textbf{P} \in \mathbb{R}^{c \times k}$ and $\textbf{W} \in \mathbb{R}^{k \times f}$ are the cell-wise program scores and feature weights per program respectively for $k$ programs. The programs represent the coordinated activity of genes towards a biological function. These models rely on the correlation structure of the data to estimate these matrices. If high spurious correlation is present in the data, the models can report latent components that do not reflect biological activity. Furthermore, while the models are purely correlative (not causal), training on interventional data implicitly benefits the model due to the alterations in the correlation structure of the data.

A few examples of such models are [cNMF](https://github.com/dylkot/cNMF), [LDVAE](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html), [f-scLVM](https://github.com/scfurl/f-scLVM). More advanced models consider additional inputs or work on multi-omic data such as [Spectra](https://github.com/dpeerlab/spectra/), [muVI](https://github.com/MLO-lab/MuVI). 

# Evaluating the biological basis of learnt programs
The goal is to develop evaluation criteria covering various technical and biological facets for programs inferred from these models. This process includes both coming up with a concept and then an implementation. A list of implemented or planned criteria can be found below.

# Structure of the evaluation pipeline
* Single-cell omics data and program scores are stored in the mudata format (see [mudata documentation](https://mudata.readthedocs.io/en/latest/)).
* <ins>Templates for implementing methods or evaluations can be found in **src/evaluation/** or **src/inference/program_methods/**
  
![image](https://github.com/EngreitzLab/gene_program_evaluation/assets/25486108/eb14e159-de07-47bd-87d4-12d5b94102fa)

| Criterion    | Implementation | External resource | Interpretation | Caveats |
| -------- | ------- | -------- | ------- | ------- |
| Goodness of fit  | [Explained variance](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html) per program    | None | A program explaining more variance in the data might represent dominant biological variation. | Technical variation might be the highest source of variance (e.g. batch effects). |
| Variation across category levels | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are variable between batch levels then the component likely is modelling technical noise. Alternatively, if program scores are variable between a biological category like cell-type or condition then the program is likely modelling a biological process specific to the category. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Gene-set enrichment | [GSEA](https://gseapy.readthedocs.io/en/latest/introduction.html) using program x feature scores | MsigDB, Enrichr | If a program is significantly associated with a gene-set then it could explain the biological process the program represents | |
| Motif enrichment | Pearson correlation of motif counts per gene (promoter or enchancer) and program x gene scores | HOCOMOCO v12 | If genes with high contributions to a program are also enriched with same enhancer/promoter motifs they could be co-regulated | A biological pathway could involve genes with different regulation but still contribute to a common function | 
| Co-regulation |  | TF-gene links | If a program contains genes that are regulated by TFs that form a regulatory module then it indicates mechansitic commonality  |  | 
| Perturbation sensitivity |  | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 
| Cross-modality prediction |  | Multi-omic data | If program x feature scores learnt from one modality can lead to a good fit for another modality/dataset indicate a robust biological connection b/w genes in the program | Technical variation b/w datasets and mapping features b/w modalities are practical challenges | 

import os
import pandas as pd


def parse_methods(mdata, data_key="rna"):
    methods = {}
    n_components = {}
    for key in mdata.mod.keys():
        method_split = key.split("_")
        if len(method_split) > 1:
            method = "_".join(method_split[:-1])
        else:
            method = method_split[0]
        if method != data_key:
            methods[key] = method
            n_components[method] = mdata.mod[key].X.shape[1]
    return methods, n_components


def parse_explained_variance(subdirs):
    explained_variance_ratios = {}
    cumulative_explained_variance = {}
    for subdir in subdirs:
        try:
            run_name = os.path.basename(subdir)
            df = pd.read_csv(os.path.join(subdir, "explained_variance_ratio.txt"), sep="\t")
            df.columns = ["program_name", "explained_variance_ratio"]
            explained_variance_ratios[run_name] = df
            cumulative_explained_variance[run_name] = df["explained_variance_ratio"].sum()
        except FileNotFoundError:
            print(f"File not found: {subdir}")
    return explained_variance_ratios, cumulative_explained_variance


def parse_loadings(mdata, data_key="rna"):
    loadings = {}
    for key in mdata.mod.keys():
        if key != data_key:
            loadings[key] = pd.DataFrame(
                data=mdata.mod[key].varm["loadings"],
                index=mdata.mod[key].var_names,
                columns=mdata.mod[key].uns["var_names"]
            )
            loadings[key].index.name = "program_name"
            loadings[key].columns.name = "gene_name"
    return loadings


def parse_geneset_enrichments(subdirs):
    geneset_enrichments = {}
    for subdir in subdirs:
        try:
            run_name = os.path.basename(subdir)
            gene_set_enrichment_file = os.path.join(subdir, "geneset_enrichment.txt")
            gene_set_enrichment_df = pd.read_csv(gene_set_enrichment_file, sep="\t")
            geneset_enrichments[run_name] = gene_set_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {gene_set_enrichment_file}")
            continue
    return geneset_enrichments


def parse_motif_enrichments(subdirs):
    motif_enrichments = {}
    for subdir in subdirs:
        try:
            run_name = os.path.basename(subdir)
            motif_enrichment_file = os.path.join(subdir, "motif_enrichment.txt")
            motif_enrichment_df = pd.read_csv(motif_enrichment_file, sep="\t")
            motif_enrichments[run_name] = motif_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {motif_enrichment_file}")
            continue
    return motif_enrichments


def parse_trait_enrichments(subdirs):
    trait_enrichments = {}
    for subdir in subdirs:
        try:
            run_name = os.path.basename(subdir)
            trait_enrichment_file = os.path.join(subdir, "trait_enrichment.txt")
            trait_enrichment_df = pd.read_csv(trait_enrichment_file, sep="\t")
            trait_enrichments[run_name] = trait_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {trait_enrichment_file}")
            continue
    return trait_enrichments


def parse_perturbation_associations(subdirs):
    perturbation_associations = {}
    for subdir in subdirs:
        try:
            run_name = os.path.basename(subdir)
            perturbation_association_file = os.path.join(subdir, "perturbation_association_results.txt")
            perturbation_association_df = pd.read_csv(perturbation_association_file, sep="\t")
            perturbation_associations[run_name] = perturbation_association_df
        except FileNotFoundError:
            print(f"File not found: {perturbation_association_file}")
            continue
    return perturbation_associations


def parse(
    mdata,
    subdirs,
    data_key="rna",
):
    methods, n_components = parse_methods(mdata, data_key)
    explained_variance_ratios, cumulative_explained_variance = parse_explained_variance(subdirs)
    loadings = parse_loadings(mdata, data_key)
    geneset_enrichments = parse_geneset_enrichments(subdirs)
    motif_enrichments = parse_motif_enrichments(subdirs)
    trait_enrichments = parse_trait_enrichments(subdirs)
    perturbation_associations = parse_perturbation_associations(subdirs)
    return {
        "methods": methods,
        "n_components": n_components,
        "explained_variance_ratios": explained_variance_ratios,
        "cumulative_explained_variance": cumulative_explained_variance,
        "loadings": loadings,
        "geneset_enrichments": geneset_enrichments,
        "motif_enrichments": motif_enrichments,
        "trait_enrichments": trait_enrichments,
        "perturbation_associations": perturbation_associations,
    }
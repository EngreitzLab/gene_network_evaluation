import os
import pandas as pd
from utils import load_config


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


def parse_obs(mdata, data_key="rna"):
    obs = {}
    for key in mdata.mod.keys():
        if key != data_key:
            obs[key] = mdata.mod[key].obs
    return obs


def parse_obs_memberships(mdata, data_key="rna"):
    obs_memberships = {}
    for key in mdata.mod.keys():
        if key != data_key:
            obs_memberships[key] = pd.DataFrame(
                data=mdata.mod[key].X,
                index=mdata.mod[key].obs_names,
                columns=mdata.mod[key].var_names
            )
            obs_memberships[key].index.name = "obs_name"
            obs_memberships[key].columns.name = "program_name"
    return obs_memberships


def parse_explained_variance(dirs):
    explained_variance_ratios = {}
    cumulative_explained_variance = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            explained_variance_ratio_file = os.path.join(dir, "explained_variance_ratio.txt")
            df = pd.read_csv(explained_variance_ratio_file, sep="\t")
            df.columns = ["program_name", "explained_variance_ratio"]
            explained_variance_ratios[run_name] = df
            cumulative_explained_variance[run_name] = df["explained_variance_ratio"].sum()
        except FileNotFoundError:
            print(f"File not found: {explained_variance_ratio_file}")
    return explained_variance_ratios, cumulative_explained_variance


def parse_categprocal_associations(dirs):
    categorical_associations = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            categorical_association_file = os.path.join(dir, "categorical_association_results.txt")
            categorical_association_df = pd.read_csv(categorical_association_file, sep="\t")
            categorical_associations[run_name] = categorical_association_df
        except FileNotFoundError:
            print(f"File not found: {categorical_association_file}")
            continue
    return categorical_associations


def parse_geneset_enrichments(dirs):
    geneset_enrichments = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            gene_set_enrichment_file = os.path.join(dir, "geneset_enrichment.txt")
            gene_set_enrichment_df = pd.read_csv(gene_set_enrichment_file, sep="\t")
            geneset_enrichments[run_name] = gene_set_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {gene_set_enrichment_file}")
            continue
    return geneset_enrichments


def parse_motif_enrichments(dirs):
    motif_enrichments = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            motif_enrichment_file = os.path.join(dir, "motif_enrichment.txt")
            motif_enrichment_df = pd.read_csv(motif_enrichment_file, sep="\t")
            motif_enrichments[run_name] = motif_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {motif_enrichment_file}")
            continue
    return motif_enrichments


def parse_trait_enrichments(dirs):
    trait_enrichments = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            trait_enrichment_file = os.path.join(dir, "trait_enrichment_processed.txt")
            trait_enrichment_df = pd.read_csv(trait_enrichment_file, sep="\t")
            trait_enrichments[run_name] = trait_enrichment_df
        except FileNotFoundError:
            print(f"File not found: {trait_enrichment_file}")
            continue
    return trait_enrichments


def parse_perturbation_associations(dirs):
    perturbation_associations = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            perturbation_association_file = os.path.join(dir, "perturbation_association_results.txt")
            perturbation_association_df = pd.read_csv(perturbation_association_file, sep="\t")
            perturbation_associations[run_name] = perturbation_association_df
        except FileNotFoundError:
            print(f"File not found: {perturbation_association_file}")
            continue
    return perturbation_associations


def parse_software_versions(dirs):
    software_versions = {}
    for dir in dirs:
        try:
            run_name = os.path.basename(dir)
            software_versions_file = os.path.join(dir, "software_versions.yml")
            software_versions_df = load_config(software_versions_file)
            software_versions[run_name] = software_versions_df
        except FileNotFoundError:
            print(f"File not found: {software_versions_file}")
            continue
    return software_versions


def parse(
    mdata,
    dirs,
    data_key="rna",
):
    methods, n_components = parse_methods(mdata, data_key)
    loadings = parse_loadings(mdata, data_key)
    obs = parse_obs(mdata, data_key)
    obs_memberships = parse_obs_memberships(mdata, data_key)
    explained_variance_ratios, cumulative_explained_variance = parse_explained_variance(dirs)
    categorical_associations = parse_categprocal_associations(dirs)
    geneset_enrichments = parse_geneset_enrichments(dirs)
    motif_enrichments = parse_motif_enrichments(dirs)
    trait_enrichments = parse_trait_enrichments(dirs)
    perturbation_associations = parse_perturbation_associations(dirs)
    software_versions = parse_software_versions(dirs)
    return {
        "methods": methods,
        "n_components": n_components,
        "cumulative_explained_variance": cumulative_explained_variance,
        "loadings": loadings,
        "obs": obs,
        "obs_memberships": obs_memberships,
        "explained_variance_ratios": explained_variance_ratios,
        "categorical_associations": categorical_associations,
        "geneset_enrichments": geneset_enrichments,
        "motif_enrichments": motif_enrichments,
        "trait_enrichments": trait_enrichments,
        "perturbation_associations": perturbation_associations,
        "software_versions": software_versions,
    }

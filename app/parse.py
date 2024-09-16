import os
import glob
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
            n_components[key] = mdata.mod[key].X.shape[1]
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


def parse_categorical_associations(dirs, prog_keys):
    categorical_associations_results = {}
    categorical_associations_posthoc = {}
    for dir, prog_key in zip(dirs, prog_keys):
        
        # Initialize dictionaries 
        categorical_associations_results[prog_key] = {}
        categorical_associations_posthoc[prog_key] = {}

        # Load association results
        categorical_association_files = glob.glob(os.path.join(dir, f"{prog_key}_*_association_results.txt"))
        for categorical_association_file in categorical_association_files:
            categorical_association_df = pd.read_csv(categorical_association_file, sep="\t")
            categorical_key = categorical_association_file.split(f"{prog_key}_")[1].split("_association_results.txt")[0]
            categorical_associations_results[prog_key][categorical_key] = categorical_association_df
        
        # Load association posthoc results
        categorical_association_posthoc_files = glob.glob(os.path.join(dir, f"{prog_key}_*_association_posthoc.txt"))
        for categorical_association_posthoc_file in categorical_association_posthoc_files:
            categorical_association_posthoc_df = pd.read_csv(categorical_association_posthoc_file, sep="\t")
            categorical_key = categorical_association_posthoc_file.split(f"{prog_key}_")[1].split("_association_posthoc.txt")[0]
            categorical_associations_posthoc[prog_key][categorical_key] = categorical_association_posthoc_df
            
    return categorical_associations_results, categorical_associations_posthoc


def parse_perturbation_associations(
    dirs, 
    prog_keys,
    stratification_key=None
):
    perturbation_associations = {}
    for dir, prog_key in zip(dirs, prog_keys):
        if prog_key not in perturbation_associations:
            perturbation_associations[prog_key] = {}
        perturbation_associations[prog_key]["results"] = {}
        perturbation_associations[prog_key]["gene_guides"] = []
        perturbation_associations[prog_key]["stratification_keys"] = []
        perturbation_associations[prog_key]["level_keys"] = []
        perturbation_association_files = glob.glob(os.path.join(dir, f"{prog_key}_*_perturbation_association.txt"))
        for perturbation_association_file in perturbation_association_files:
            gene_guide = perturbation_association_file.split(f"{prog_key}_")[1].split("_")[0]
            if not stratification_key:
                curr_stratification_key = "global"
                curr_level_key = "global"
            else:
                curr_stratification_key = stratification_key
                curr_level_key = perturbation_association_file.split(f"{prog_key}_{gene_guide}_{curr_stratification_key}_")[1].split("_perturbation_association.txt")[0]
            print(f"Gene/guide: {gene_guide}, Stratification key: {curr_stratification_key}, Level key: {curr_level_key}")
            df = pd.read_csv(perturbation_association_file, sep="\t")
            perturbation_associations[prog_key]["results"][f"{gene_guide}_{curr_stratification_key}_{curr_level_key}"] = df
            perturbation_associations[prog_key]["gene_guides"].append(gene_guide)
            perturbation_associations[prog_key]["stratification_keys"].append(curr_stratification_key)
            perturbation_associations[prog_key]["level_keys"].append(curr_level_key)
    return perturbation_associations


def parse_geneset_enrichments(dirs, prog_keys):
    geneset_enrichments = {}
    for dir, prog_key in zip(dirs, prog_keys):
        if prog_key not in geneset_enrichments:
            geneset_enrichments[prog_key] = {}
        geneset_enrichments[prog_key]["results"] = {}
        geneset_enrichments[prog_key]["libraries"] = []
        geneset_enrichments[prog_key]["methods"] = []
        geneset_enrichment_files = glob.glob(os.path.join(dir, f"{prog_key}_*_geneset_enrichment.txt"))
        for gene_set_enrichment_file in geneset_enrichment_files:
            method = gene_set_enrichment_file.split("_geneset_enrichment.txt")[0].split("_")[-1]
            library = gene_set_enrichment_file.split(f"{prog_key}_")[1].split(f"_{method}_geneset_enrichment.txt")[0]
            print(f"Library: {library}, Method: {method}")
            df = pd.read_csv(gene_set_enrichment_file, sep="\t")
            geneset_enrichments[prog_key]["results"][f"{library}_{method}"] = df
            geneset_enrichments[prog_key]["libraries"].append(library)
            geneset_enrichments[prog_key]["methods"].append(method)
    return geneset_enrichments


def parse_trait_enrichments(dirs, prog_keys):
    trait_enrichments = {}
    for dir, prog_key in zip(dirs, prog_keys):
        if prog_key not in trait_enrichments:
            trait_enrichments[prog_key] = {}
        trait_enrichments[prog_key]["results"] = {}
        trait_enrichments[prog_key]["databases"] = []
        trait_enrichments[prog_key]["methods"] = []
        trait_enrichment_files = glob.glob(os.path.join(dir, f"{prog_key}_*_trait_enrichment.txt"))
        for trait_enrichment_file in trait_enrichment_files:
            method = trait_enrichment_file.split("_trait_enrichment.txt")[0].split("_")[-1]
            database = trait_enrichment_file.split(f"{prog_key}_")[1].split(f"_{method}_trait_enrichment.txt")[0]
            print(f"Database: {database}, Method: {method}")
            df = pd.read_csv(trait_enrichment_file, sep="\t")
            trait_enrichments[prog_key]["results"][f"{database}_{method}"] = df
            trait_enrichments[prog_key]["databases"].append(database)
            trait_enrichments[prog_key]["methods"].append(method)
    return trait_enrichments
    
    
def parse_motif_enrichments(
    dirs, 
    prog_keys,
    stratification_key=None
):
    motif_enrichments = {}
    for dir, prog_key in zip(dirs, prog_keys):
        if prog_key not in motif_enrichments:
            motif_enrichments[prog_key] = {}
        motif_enrichments[prog_key]["results"] = {}
        motif_enrichments[prog_key]["E_P_types"] = []
        motif_enrichments[prog_key]["databases"] = []
        motif_enrichments[prog_key]["test_types"] = []
        motif_enrichments[prog_key]["stratification_keys"] = []
        motif_enrichments[prog_key]["level_keys"] = []
        motif_enrichment_files = glob.glob(os.path.join(dir, f"{prog_key}_*_motif_enrichment.txt"))
        for motif_enrichment_file in motif_enrichment_files:
            E_P_type = motif_enrichment_file.split(f"{prog_key}_")[1].split("_")[0]
            database = motif_enrichment_file.split(f"{prog_key}_{E_P_type}_")[1].split("_")[0]
            test_type = motif_enrichment_file.split(f"{prog_key}_{E_P_type}_{database}_")[1].split("_")[0]
            if not stratification_key:
                curr_stratification_key = "global"
                curr_level_key = "global"
            else:
                curr_level_key = motif_enrichment_file.split(f"{prog_key}_{E_P_type}_{database}_{test_type}_{curr_stratification_key}_")[1].split("_motif_enrichment.txt")[0]
                curr_stratification_key = stratification_key
            print(f"E_P_type: {E_P_type}, Database: {database}, Test type: {test_type}, Stratification key: {curr_stratification_key}, Level key: {curr_level_key}")
            df = pd.read_csv(motif_enrichment_file, sep="\t")
            motif_enrichments[prog_key]["results"][f"{E_P_type}_{database}_{test_type}_{curr_stratification_key}_{curr_level_key}"] = df
            motif_enrichments[prog_key]["E_P_types"].append(E_P_type)
            motif_enrichments[prog_key]["databases"].append(database)
            motif_enrichments[prog_key]["test_types"].append(test_type)
            motif_enrichments[prog_key]["stratification_keys"].append(curr_stratification_key)
            motif_enrichments[prog_key]["level_keys"].append(curr_level_key)
    
    return motif_enrichments


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


def parse_explained_variance(dirs, prog_keys):
    explained_variance_ratios = {}
    for dir, prog_key in zip(dirs, prog_keys):
        try:
            explained_variance_ratio_file = os.path.join(dir, f"{prog_key}_variance_explained_ratio.txt")
            df = pd.read_csv(explained_variance_ratio_file, sep="\t")
            explained_variance_ratios[prog_key] = df
        except FileNotFoundError:
            print(f"File not found: {explained_variance_ratio_file}")
    return explained_variance_ratios


def parse(
    mdata,
    dirs,
    data_key="rna",
    perturbation_association_stratification_key=None,
    motif_enrichment_stratification_key=None
):
    methods, n_components = parse_methods(mdata, data_key)
    loadings = parse_loadings(mdata, data_key)
    obs_memberships = parse_obs_memberships(mdata, data_key)
    categorical_associations_results, categorical_associations_posthoc = parse_categorical_associations(dirs, methods.keys())
    perturbation_associations = parse_perturbation_associations(dirs, methods.keys(), perturbation_association_stratification_key)
    geneset_enrichments = parse_geneset_enrichments(dirs, methods.keys())
    trait_enrichments = parse_trait_enrichments(dirs, methods.keys())
    motif_enrichments = parse_motif_enrichments(dirs, methods.keys(), motif_enrichment_stratification_key)
    explained_variance_ratios = parse_explained_variance(dirs, methods.keys())
    software_versions = parse_software_versions(dirs)
    return {
        "methods": methods,
        "n_components": n_components,
        "loadings": loadings,
        "obs_memberships": obs_memberships,
        "categorical_associations_results": categorical_associations_results,
        "categorical_associations_posthoc": categorical_associations_posthoc,
        "perturbation_associations": perturbation_associations,
        "geneset_enrichments": geneset_enrichments,
        "trait_enrichments": trait_enrichments,
        "motif_enrichments": motif_enrichments,
        "explained_variance_ratios": explained_variance_ratios,
        "software_versions": software_versions,
    }

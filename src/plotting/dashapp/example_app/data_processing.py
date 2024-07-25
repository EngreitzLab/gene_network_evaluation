import numpy as np
import pandas as pd


def count(categorical_var, count_var, dataframe):
    counts_df = dataframe.value_counts([categorical_var, count_var])
    counts_df = counts_df.groupby(categorical_var).sum()
    counts_df = counts_df.sort_values(ascending=False)
    counts_df = pd.DataFrame(counts_df.reset_index().values,
                             columns=[categorical_var, count_var])
    return counts_df


def count_unique(categorical_var, count_var, dataframe, cummul=False, unique=False):
    counts_df = count(categorical_var, count_var, dataframe)
    new_df = []
    terms = []
    for prog in counts_df[categorical_var].unique():
        terms_ = dataframe.loc[dataframe[categorical_var] == prog, count_var].unique()
        unique_terms = [term for term in terms_ if term not in terms]
        terms.extend(unique_terms)
        new_df.append([prog, len(unique_terms)])
    new_df = pd.DataFrame(new_df, columns=[categorical_var, count_var])
    if cummul:
        new_df[count_var] = new_df[count_var].cumsum()
    if unique:
        return counts_df
    else:
        return new_df


def filter_data(data: pd.DataFrame, filters: dict):
    """Filter data based on a dictionary of filter criteria.
    
    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing the data to be filtered.
    filters : dict
        Dictionary of filter criteria where keys are column names and values are the filter values.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    filtered_data = data.copy()
    for column, value in filters.items():
        if isinstance(value, tuple) and len(value) == 2:
            filtered_data = filtered_data[(filtered_data[column] >= value[0]) & (filtered_data[column] <= value[1])]
        else:
            filtered_data = filtered_data[filtered_data[column] <= value]
    return filtered_data


def extract_total_unique_counts(
    cross_run,
    summary,
    keys=["enrichment_gsea", "enrichment_motif", "enrichment_gwas"],
    id_keys=["ID", "EPType-TFMotif", "Term"],
    filter_cols=["qvalue", "FDR", "Adjusted_P-value"],
    filter_thresholds=[0.05, 0.05, 0.05],
    verbose=False
):
    """Extract the total number of unique IDs for each modality and enrichment type.
    
    This function will count the number of unique IDs for each modality and enrichment type.
    It will also filter the data based on the filter_cols and filter_thresholds.
    """

    # Create a dictionary to store the total unique counts
    total_unique_ids = {}
    
    # For each enrichment type, count the number of unique IDs
    for i in range(len(keys)):
        key = keys[i]
        id_key = id_keys[i]
        unique_ids = [cross_run[mod][key][id_key].nunique() for mod in cross_run.keys()]
        mod_names = list(cross_run.keys())
        order = np.argsort(unique_ids)
        mod_names = [mod_names[i] for i in order]
        unique_ids = [unique_ids[i] for i in order]
        total_unique_ids[key] = unique_ids
    n_programs = [summary[mod]["n_programs"] for mod in mod_names]
    total_unique_ids["prog_names"] = mod_names
    total_unique_ids["n_programs"] = n_programs
    
    df = pd.DataFrame(total_unique_ids, index=mod_names)
    return df


def extract_enrichment(
    adata, 
    library="GSEA", 
    geneset_index="Term", 
    program_index="program_name",
    varmap_name_prefix="gsea_varmap",
    verbose=False
):
    # If verbose
    if verbose:
        print(f"Extracting enrichment data for {library} library "
                f"with geneset index {geneset_index} and program index {program_index} "
                f"using varmap prefix {varmap_name_prefix}")
    
    # Create a mudata key to column name mapping dictionary
    mudata_keys_dict = {}
    for key in adata.varm.keys():
        if library in key:
            colname = key.replace(f"_{library}", "")
            mudata_keys_dict[colname] = key
    
    # Unpivot the programs and genesets
    programs = adata.var.index.tolist()
    genesets = adata.uns[f"{varmap_name_prefix}_{library}"].tolist()
    
    # Create an empty dataframe
    df = pd.DataFrame(index=genesets, columns=programs)

    # Melting the dataframe
    df = df.melt(value_vars=programs, var_name=program_index, value_name="value", ignore_index=False).reset_index()
    df.columns = [geneset_index, program_index, "value"]
    df = df.drop(columns="value")
    
    # For each key in the dictionary, extract the values from the mudata object and put them into the dataframe
    for colname, key in mudata_keys_dict.items():
        
        # Extract the values from the mudata object
        all_progs_array = adata.varm[key].flatten()
        
        # Add the values to the dataframe
        df[colname] = all_progs_array

    # Drop any rows with NaN values
    df = df.dropna()
    
    return df


def parse_mdata_summary(
    mdata,
    mod_keys=None,
    verbose=False
):
    """Extract the summary data from the mdata object.
    
    Will parse all the modalities and return a nested dictionary:
    - First level keys are the modalities
    - Second level keys are the summary stats for the data:
        - n_genes
        - n_cells
        - n_programs
        - n_counts for each .obs column
    """
    
    # If no keys are provided, extract all keys
    if mod_keys is None:
        mod_keys = list(mdata.mod.keys())

    # Create a dictionary to store the dataframes
    data_dict = {}

    # Grab number of cells
    n_cells = mdata.n_obs
    data_dict["n_cells"] = n_cells
    
    # Grab obs
    obs = mdata.obs
    data_dict["obs"] = obs

    # For each modality, extract the summary data
    for mod in mod_keys:
        if verbose:
            print(f"Extracting summary data for modality {mod}")
        
        # Extract the adata object
        adata = mdata[mod]
        
        # Extract the summary data
        n_programs = adata.n_vars
        
        # Store the dataframes in a dictionary
        data_dict[mod] = {
            "n_programs": n_programs,
        }
    
    return data_dict


def parse_mdata_cross_run(
    mdata,
    blacklist=["rna", "atac"],
    geneset_library="GSEA",
    motif_library="Motif",
    trait_library="GWAS",
    verbose=False
):
    """Parse the mdata object for cross run evaluation.
    
    This function will extract all the modalities from the mdata object not part of the blacklist
    and return a nested dictionary:
    - First level keys are the modalities
    - Second level keys are the evaluation types (explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait)
    - Values are the corresponding dataframes

    It will attempt to extract variance explained from .var
    It will attempt to extract enrichment data from .varm:
    1. geneset enrichment
    2. motif enrichment 
    3. trait enrichment
    """

    # Extract the modalities from the mdata object
    modalities = [mod for mod in mdata.mod.keys() if not any([x in mod for x in blacklist])]
    if verbose:
        print(f"Extracting data for the following modalities: {modalities}")
        
    # Create a dictionary to store the dataframes
    data_dict = {}

    for mod in modalities:
        if verbose:
            print(f"Extracting data for modality {mod}")

        # Extract the adata object
        adata = mdata[mod]

        # Extract the explained variance
        if "VarianceExplained" in adata.var.keys():
            explained_variance = pd.DataFrame({
                "ProgramID": adata.var.index,
                "VarianceExplained": adata.var["VarianceExplained"]
            })
        else:
            explained_variance = None

        # Extract the enrichment data
        enrichment_gsea = extract_enrichment(
            adata, 
            library=geneset_library,
            geneset_index="ID", 
            program_index="ProgramID",
            varmap_name_prefix="gsea_varmap",
            verbose=verbose
        )
        enrichment_motif = extract_enrichment(
            adata, 
            library=motif_library,
            geneset_index="EPType-TFMotif",
            program_index="ProgramID",
            varmap_name_prefix="motif_varmap",
            verbose=verbose
        )
        enrichment_trait = extract_enrichment(
            adata,
            library=trait_library,
            geneset_index="Term",
            program_index="ProgramID",
            varmap_name_prefix="gwas_varmap",
            verbose=verbose
        )

        # Store the dataframes in a dictionary
        data_dict[mod] = {
            "explained_variance": explained_variance,
            "enrichment_gsea": enrichment_gsea,
            "enrichment_motif": enrichment_motif,
            "enrichment_gwas": enrichment_trait
        }

    return data_dict


# Deprecated
def load_example_data(file_path):
    evaluation_output = pd.ExcelFile(file_path)
    explained_variance = pd.read_excel(evaluation_output, sheet_name='explained_variance')
    enrichment_gsea = pd.read_excel(evaluation_output, sheet_name='enrichment_gsea')
    enrichment_motif = pd.read_excel(evaluation_output, sheet_name='enrichment_motif')
    enrichment_trait = pd.read_excel(evaluation_output, sheet_name='enrichment_trait')
    return explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait

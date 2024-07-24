import pandas as pd


def load_data(file_path):
    evaluation_output = pd.ExcelFile(file_path)
    explained_variance = pd.read_excel(evaluation_output, sheet_name='explained_variance')
    enrichment_gsea = pd.read_excel(evaluation_output, sheet_name='enrichment_gsea')
    enrichment_motif = pd.read_excel(evaluation_output, sheet_name='enrichment_motif')
    enrichment_trait = pd.read_excel(evaluation_output, sheet_name='enrichment_trait')
    return explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait


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

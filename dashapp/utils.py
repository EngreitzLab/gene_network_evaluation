from typing import List, Dict
import re
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


def infer_dashboard_type(keys: List[str]) -> str:
    if len(keys) == 1:
        return "single_run"
    
    base_names = {}
    for key in keys:
        match = re.match(r'([a-zA-Z]+)_?(\d+)?', key)
        if match:
            base_name, num = match.groups()
            if base_name not in base_names:
                base_names[base_name] = []
            if num:
                base_names[base_name].append(int(num))
    
    if len(base_names) == 1:
        return "cross_k"
    elif all(len(nums) == 0 for nums in base_names.values()):
        return "cross_method"
    else:
        return "mixed"
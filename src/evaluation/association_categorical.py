import os
import argparse

import mudata
import numpy as np
import pandas as pd

from scipy import stats, sparse
from scikit_posthocs import posthoc_dscf, posthoc_conover, posthoc_dunn
from statsmodels.stats.multitest import fdrcorrection

from typing import Literal
from joblib import Parallel, delayed
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)

# Perform non-parametric test (~one way ANOVA) to compare prog scores across categorical levels
def perform_kruskall_wallis(mdata, prog_key='prog', prog_name=None,  
                            categorical_key='batch', pseudobulk_key='sample'):
    """
    Performs non-parametric one-way ANOVA using specified categorical levels.
    If pseudobulk key is provided then data is averaged over those levels.
    Mudata object is updated in place with KW test results.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        prog_name: str
            index of the feature (mdata[prog_key][:, prog_name]) which is the response variable.
        categorical_key: str
            index of the categorical levels (mdata[prog_key].obs[categorical_key]) being tested.
        pseudobulk_key: str (optional)
            index of the feature summarisation (mean) levels (mdata[prog_key].obs[pseudobulk_key]).
    
    UPDATES
        if pseudobulk_key is None:
            store_key = categorical_key
        else:
            store_key = '_'.join(categorical_key, pseudobulk_key)
        mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_stat'.format(store_key)]  
        mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_pval'.format(store_key)] 
        mdata[prog_key].uns['{}_association_categories'.format(categorical_key)]

    """
    samples = []
    if pseudobulk_key is None:
        # Run with single-cells are individual data points
        for category in mdata[prog_key].obs[categorical_key].astype(str).unique():
            sample_ = mdata[prog_key][mdata[prog_key].obs[categorical_key].astype(str)==category, 
                                        prog_name].X[:,0]
            if sparse.issparse(sample_):
                sample_ = sample_.toarray().flatten()
            samples.append(sample_)
    else:
        # Pseudobulk the data before performing KW test (ideally these are biological replicates)
        if sparse.issparse(mdata[prog_key][:, prog_name].X[:,0]):
            prog_data = mdata[prog_key][:, prog_name].X[:,0].toarray()
        else:
            prog_data = mdata[prog_key][:, prog_name].X[:,0]
        prog_data = pd.DataFrame(prog_data, columns=[prog_name], index=mdata[prog_key].obs.index)
        prog_data[categorical_key] = mdata[prog_key].obs[categorical_key]
        prog_data[pseudobulk_key] = mdata[prog_key].obs[pseudobulk_key]
        prog_data = prog_data.groupby([pseudobulk_key, categorical_key]).mean().dropna().reset_index()

        for category in mdata[prog_key].obs[categorical_key].astype(str).unique():
            if len(prog_data.loc[prog_data[categorical_key]==category, pseudobulk_key].unique()) < 3:
                raise ValueError('Cannot use less than 3 replicates per categorical level')
            sample_ = prog_data.loc[prog_data[categorical_key].astype(str)==category, prog_name].values.flatten()
            samples.append(sample_)

    stat, pval = stats.kruskal(*samples, nan_policy='propagate')

    # Store results
    if pseudobulk_key is None:
        store_key = categorical_key
    else:
        store_key = '_'.join((categorical_key, pseudobulk_key))

    mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_stat'.format(store_key)] = stat
    mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_pval'.format(store_key)] = pval

    mdata[prog_key].uns['{}_association_categories'.format(categorical_key)] = \
    mdata[prog_key].obs[categorical_key].astype(str).unique()


def perform_correlation(
    prog_df: pd.DataFrame,
    group_col: str,
    val_col: str,
    low_threshold: float=0.0,
    correlation: Literal['pearsonr', 'spearmanr', 'kendalltau']='pearsonr',
    mode: Literal['one_vs_all', 'one_vs_one']='one_vs_all',
    df = [],
):
    """Perform post-hoc analysis with correlation tests

    Parameters
    ----------
    prog_df : pd.DataFrame
        DataFrame containing program scores and categorical information.
    group_col : str
        Column name of the categorical information.
    val_col : str
        Column name of the program scores.
    low_threshold : float
        Threshold to remove cells with low program scores.
    correlation : str
        Type of correlation test to perform.
    mode : str
        Type of correlation test to perform.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of the correlation tests:
        - program_name
        - {group_col}_{group_col_category_1}_{correlation}_stat
        - {group_col}_{group_col_category_1}_{correlation}_pval
        - {group_col}_{group_col_category_1}_{correlation}_adj_pval
        - {group_col}_{group_col_category_1}_log2FC
    """

    # Get unique categories
    categories = prog_df[group_col].unique().tolist()

    # Remove cell memberships that are below the threshold
    prog_df = prog_df.loc[prog_df[val_col] > low_threshold]
    
    # Perform correlation tests
    if mode=='one_vs_all':
        stats_df = pd.DataFrame(index=categories, columns=['stat', 'pval', 'adj_pval', "log2FC"], dtype=float)
        for idx in stats_df.index.values:
            if correlation=='pearsonr':
                prog_df['binarized'] = prog_df[group_col].apply(lambda x: 1 if x==idx else 0)
                stat, pval = stats.pearsonr(prog_df[val_col], prog_df['binarized'])
                _, adj_pval = fdrcorrection([pval])
                log2FC = np.log2(prog_df.loc[prog_df[group_col]==idx, val_col].mean() / prog_df.loc[prog_df[group_col]!=idx, val_col].mean())
                stats_df.loc[idx, 'stat'] = stat
                stats_df.loc[idx, 'pval'] = pval
                stats_df.loc[idx, 'adj_pval'] = adj_pval
                stats_df.loc[idx, 'log2FC'] = log2FC
            elif correlation=='spearmanr':
                stat, pval = stats.spearmanr(prog_df[val_col], prog_df[group_col].astype('category').cat.codes)
                _, adj_pval = fdrcorrection([pval])
                log2FC = np.log2(prog_df.loc[prog_df[group_col]==idx, val_col].mean() / prog_df.loc[prog_df[group_col]!=idx, val_col].mean())
                stats_df.loc[idx, 'stat'] = stat
                stats_df.loc[idx, 'pval'] = pval
                stats_df.loc[idx, 'adj_pval'] = adj_pval
                stats_df.loc[idx, 'log2FC'] = log2FC
            elif correlation=='kendalltau':
                stat, pval = stats.kendalltau(prog_df[val_col], prog_df[group_col].astype('category').cat.codes)
                _, adj_pval = fdrcorrection([pval])
                log2FC = np.log2(prog_df.loc[prog_df[group_col]==idx, val_col].mean() / prog_df.loc[prog_df[group_col]!=idx, val_col].mean())
                stats_df.loc[idx, 'stat'] = stat
                stats_df.loc[idx, 'pval'] = pval
                stats_df.loc[idx, 'adj_pval'] = adj_pval
                stats_df.loc[idx, 'log2FC'] = log2FC
            
        # Format
        wide_df = pd.DataFrame()
        columns = ['stat', 'pval', 'adj_pval', 'log2FC']
        for idx, row in stats_df.iterrows():
            for col in columns:
                wide_df[f'{group_col}_{idx}_{correlation}_{col}'] = [row[col]]
        wide_df.index = [val_col]
        wide_df.index.name = 'program_name'
        df.append(wide_df)
    
    elif mode=='one_vs_one':
        pvals = pd.DataFrame(index=categories, columns=categories)

        for idx in pvals.index.values:
            for col in pvals.columns.values:
                if idx==col:
                    continue
                else:
                    test_df = prog_df.loc[prog_df[group_col].isin([idx, col])]
                    if correlation=='pearsonr':
                        stat, pval = stats.pearsonr(test_df[val_col], 
                                                    test_df[group_col].astype('category').cat.codes)
                        pvals.loc[idx, col] = pval
                    elif correlation=='kendalltau':
                        stat, pval = stats.kendalltau(test_df[val_col], 
                                                      test_df[group_col].astype('category').cat.codes)
                        pvals.loc[idx, col] = pval

        return pvals


# Perfom posthoc test and compute categorical-program score
def perform_posthoc(mdata, prog_key='prog', prog_name=None,
                    categorical_key='batch', pseudobulk_key='sample',
                    test='dunn', mode='one_vs_one', df=[]):

    """
    Performs post-hoc tests for Kruskall-Wallis test and 
    pairwise p-vals b/w categorical levels are reduced via mean and min.
    If pseudobulk key is provided then data is averaged over those levels.
    Mudata object is updated in place with KW test results. 

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        prog_name: str
            index of the feature (mdata[prog_key][:, prog_name]) which is the response variable.
        categorical_key: str
            index of the categorical levels (mdata[prog_key].obs[categorical_key]) being tested.
        pseudobulk_key: str (optional)
            index of the feature summarisation (mean) levels (mdata[prog_key].obs[pseudobulk_key]).
        test: {'conover','dunn', 'dscf', 'pearsonr', 'kendalltau'}
            posthoc test to use to test categorical levels.
    
    UPDATES
        if pseudobulk_key is None:
            store_key = categorical_key
        else:
            store_key = '_'.join(categorical_key, pseudobulk_key)
        mdata[prog_key].varm['{}_association_{}_min_pval'.format(store_key, test)]
        mdata[prog_key].varm['{}_association_{}_mean_pval'.format(store_key, test)] 
        mdata[prog_key].uns['{}_association_{}_pvals'.format(store_key, test)]

    """
    
    store_key = categorical_key

    prog_data_ = mdata[prog_key][:, prog_name].X
    if sparse.issparse(prog_data_):
        prog_data_ = prog_data_.toarray()

    prog_df = pd.DataFrame(prog_data_, 
                           columns=[prog_name], 
                           index=mdata[prog_key].obs.index)
    prog_df[categorical_key] = mdata[prog_key].obs[categorical_key].astype(str).values
    if pseudobulk_key is not None:
        store_key = '_'.join((categorical_key, pseudobulk_key))
        prog_df[pseudobulk_key] = pseudobulk_key
        prog_df = prog_df.groupby([pseudobulk_key, categorical_key]).mean().reset_index()

    # Peform post hoc tests
    if test=='dunn':
        p_vals = posthoc_dunn(prog_df, val_col=prog_name, group_col=categorical_key)
    elif test=='conover':
        p_vals = posthoc_conover(prog_df, val_col=prog_name, group_col=categorical_key)
    elif test=='dscf':
        p_vals = posthoc_dscf(prog_df, val_col=prog_name, group_col=categorical_key)
    elif test=='pearsonr':
        p_vals = perform_correlation(prog_df, val_col=prog_name, group_col=categorical_key,
                                     correlation='pearsonr', mode=mode, df=df)
        if mode == 'one_vs_all':
            return
    elif test=='kendalltau':
        p_vals = perform_correlation(prog_df, val_col=prog_name, group_col=categorical_key,
                                     correlation='kendalltau')

    # Create combined statistic with p-vals for each categorical level
    for i, category in enumerate(mdata[prog_key].obs[categorical_key].astype(str).unique()):
        min_pval = np.min(p_vals.loc[p_vals.index!=category, category])
        mean_pval = np.mean(p_vals.loc[p_vals.index!=category, category])

        prog_idx = mdata[prog_key].var.index.get_loc(prog_name)

        # Store results
        mdata[prog_key].varm['{}_association_{}_min_pval'.format(store_key, test)][prog_idx,i] = min_pval
        mdata[prog_key].varm['{}_association_{}_mean_pval'.format(store_key, test)][prog_idx,i] = mean_pval
         
    mdata[prog_key].uns['{}_association_{}_pvals'.format(store_key, test)][prog_name] = p_vals


# TODO: Add one way ANOVA, MANOVA, multi-variate KW, logistic-regression
def compute_categorical_association(
    mdata, 
    prog_key='prog', 
    categorical_key='batch', 
    pseudobulk_key=None, 
    test='dunn', 
    mode='one_vs_one',
    n_jobs=1, 
    inplace=True, 
    **kwargs
):

    """
    Compute association of continous features with categorical levels.
    Currently only the Kruskall-Wallis test is implemented to determine if
    any categorical level has a significant association and posthoc tests
    determine which categorical levels have significant associations.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        categorical_key: str
            index of the categorical levels (mdata[prog_key].obs[categorical_key]) being tested.
        pseudobulk_key: str (optional)
            index of the feature summarisation (mean) levels (mdata[prog_key].obs[pseudobulk_key]).
        test: {'conover','dunn', 'dscf', 'pearsonr', 'kendalltau'}
            posthoc test to use to test categorical levels.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
       
       ---
        if inplace:
            UPDATES
            mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_stat'.format(store_key)]  
            mdata[prog_key].var.loc[prog_name, '{}_kruskall_wallis_pval'.format(store_key)] 
            mdata[prog_key].varm['{}_association_{}_min_pval'.format(store_key, test)]
            mdata[prog_key].varm['{}_association_{}_mean_pval'.format(store_key, test)] 
            mdata[prog_key].uns['{}_association_{}_pvals'.format(store_key, test)]
            mdata[prog_key].uns['{}_association_categories'.format(categorical_key)]    
        else:
            RETURNS
            results_df, posthoc_df       

    """
    # Read in mudata if it is provided as a path
    frompath=False
    if isinstance(mdata, str):
        if os.path.exists(mdata):
            mdata = mudata.read(mdata)
            if inplace:
                logging.warning('Changed to inplace=False since path was provided')
                inplace=False
            frompath=True
        else: raise ValueError('Incorrect mudata specification.')

    if not inplace and not frompath:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy()})
    
    if pseudobulk_key is None:
        logging.info('Performing tests at single-cell level. Significance will likely be inflated')
        store_key = categorical_key
    else:
        logging.info('Perform testing by averaging over {}'.format(pseudobulk_key))
        store_key = '_'.join((categorical_key, pseudobulk_key))

    mdata[prog_key].var['{}_kruskall_wallis_stat'.format(store_key)] = None
    mdata[prog_key].var['{}_kruskall_wallis_pval'.format(store_key)] = None
    
    # Run KW test
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(perform_kruskall_wallis)(mdata, 
                                                                   prog_key=prog_key,
                                                                   prog_name=prog_name, 
                                                                   categorical_key=categorical_key,
                                                                   pseudobulk_key=pseudobulk_key) \
                                                            for prog_name in tqdm(mdata[prog_key].var_names,
                                                            desc='Testing {} association'.format(categorical_key), 
                                                            unit='programs'))

    # Convert to float to prevent error when saving mudata
    mdata[prog_key].var['{}_kruskall_wallis_stat'.format(store_key)] = \
    mdata[prog_key].var['{}_kruskall_wallis_stat'.format(store_key)].astype(float)
    mdata[prog_key].var['{}_kruskall_wallis_pval'.format(store_key)] = \
    mdata[prog_key].var['{}_kruskall_wallis_pval'.format(store_key)].astype(float)


    # If test = 'pearsonr' with mode = 'one_vs_all' we are going to do something special until we can more properly integrate this
    if test=='pearsonr' and mode=='one_vs_all':
        logging.info('Running jamboree specific version of posthoc with pearsonr, this is not yet integrated into the main pipeline')
        res = []        
        # Run posthoc-tests and append results
        Parallel(n_jobs=n_jobs,
                    backend='threading')(delayed(perform_posthoc)(mdata,
                                                                prog_key=prog_key,
                                                                prog_name=prog_name, 
                                                                categorical_key=categorical_key,
                                                                pseudobulk_key=pseudobulk_key,
                                                                test=test,
                                                                mode=mode,
                                                                df=res) \
                                                                for prog_name in tqdm(mdata[prog_key].var_names,
                                                                desc='Identifying differential {}'.format(categorical_key), 
                                                                unit='programs'))
        
        posthoc_df = pd.concat(res, axis=0)
        results_df = mdata[prog_key].var.loc[:, ['{}_kruskall_wallis_stat'.format(store_key), 
                                                 '{}_kruskall_wallis_pval'.format(store_key)]]
        posthoc_df["program_name"] = mdata[prog_key].var_names
        results_df["program_name"] = mdata[prog_key].var_names
        return (results_df, posthoc_df)

    else:

        mdata[prog_key].varm['{}_association_{}_min_pval'.format(store_key, test)] = \
        np.zeros((mdata[prog_key].shape[1], mdata[prog_key].obs[categorical_key].unique().shape[0]))
        mdata[prog_key].varm['{}_association_{}_mean_pval'.format(store_key, test)] = \
        np.ones((mdata[prog_key].shape[1], mdata[prog_key].obs[categorical_key].unique().shape[0]))
        mdata[prog_key].uns['{}_association_{}_pvals'.format(store_key, test)] = {}
        
        # Run posthoc-tests
        Parallel(n_jobs=n_jobs, 
                backend='threading')(delayed(perform_posthoc)(mdata, 
                                                            prog_key=prog_key,
                                                            prog_name=prog_name, 
                                                            categorical_key=categorical_key,
                                                            pseudobulk_key=pseudobulk_key,
                                                            test=test) \
                                                            for prog_name in tqdm(mdata[prog_key].var_names,
                                                            desc='Identifying differential {}'.format(categorical_key), 
                                                            unit='programs'))
        
        # Returning test results only
        if not inplace: 

            results_df = mdata[prog_key].var.loc[:, ['{}_kruskall_wallis_stat'.format(store_key), 
                                                    '{}_kruskall_wallis_pval'.format(store_key)]]
            min_pval_df = pd.DataFrame(mdata[prog_key].varm['{}_association_{}_min_pval'.format(store_key, test)],
                                    index=mdata[prog_key].var.index,
                                    columns=['{}_{}_association_{}_min_pval'.format(col, store_key, test) \
                                                for col in mdata[prog_key].uns['{}_association_categories'.format(categorical_key)]]
                                    )
            results_df = results_df.merge(min_pval_df, left_index=True, right_index=True)
            
            mean_pval_df = pd.DataFrame(mdata[prog_key].varm['{}_association_{}_mean_pval'.format(store_key, test)],
                                    index=mdata[prog_key].var.index,
                                    columns=['{}_{}_association_{}_mean_pval'.format(col, store_key, test) \
                                                for col in mdata[prog_key].uns['{}_association_categories'.format(categorical_key)]]
                                    )
            results_df = results_df.merge(mean_pval_df, left_index=True, right_index=True)

            posthoc_df = []
            dict_ = mdata[prog_key].uns['{}_association_{}_pvals'.format(store_key, test)]
            for key, values in dict_.items():

                tmp_df = values.where(np.triu(np.ones(values.shape), k=1).astype(bool))
                tmp_df = tmp_df.stack().reset_index()
                tmp_df.columns = ['{}_a'.format(store_key),
                                '{}_b'.format(store_key),
                                'p_value']
                tmp_df['program_name'] = key
                posthoc_df.append(tmp_df)

            posthoc_df = pd.concat(posthoc_df, ignore_index=True)

            return (results_df, posthoc_df)
        
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('categorical_key', type=str)
    parser.add_argument('--pseudobulk_key', default=None, type=str) 
    parser.add_argument('--test', default='dunn', type=str) 
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_categorical_association(mdata, prog_key=args.prog_key, categorical_key=args.categorical_key,
                                    pseudobulk_key=args.pseudobulk_key, test=args.test, n_jobs=args.n_jobs, 
                                    inplace=args.output)




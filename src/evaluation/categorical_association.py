import os
import argparse

import numpy as np
import pandas as pd

from scipy import stats, sparse
from scikit_posthocs import posthoc_dscf, posthoc_conover, posthoc_dunn

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# Perform non-parametric test (~one way ANOVA) to compare prog scores across categorical levels
def perform_kruskall_wallis(mdata, prog_nam=None, 
                            prog_key='prog', 
                            categorical_key='batch'):
    
    samples = []
    for category in mdata[prog_key].obs[categorical_key].astype(str).unique():
        sample_ = mdata[prog_key][mdata[prog_key].obs[categorical_key].astype(str)==category, 
                                  prog_nam].X[:,0]
        if sparse.issparse(sample_):
            sample_ = sample_.toarray().flatten()
        samples.append(sample_)

    stat, pval = stats.kruskal(*samples, nan_policy='propagate')

    mdata[prog_key].var.loc[prog_nam, '{}_kruskall_wallis_stat'.format(categorical_key)] = stat
    mdata[prog_key].var.loc[prog_nam, '{}_kruskall_wallis_pval'.format(categorical_key)] = pval

    mdata[prog_key].uns['{}_association_categories'.format(categorical_key)] = \
    mdata[prog_key].obs[categorical_key].astype(str).unique()

# Perfom posthoc test and compute categorical-program score
def perform_posthoc(mdata, prog_nam=None,
                    prog_key='prog', 
                    categorical_key='batch',
                    test='dunn'):
    
    prog_df = pd.DataFrame(mdata[prog_key][:, prog_nam].X, columns=[prog_nam])
    prog_df[categorical_key] = mdata[prog_key].obs[categorical_key].astype(str).values

    # Create combined statistic with p-vals for each categorical level
    if test=='dunn':
        p_vals = posthoc_dunn(prog_df, val_col=prog_nam, group_col=categorical_key)
    elif test=='conover':
        p_vals = posthoc_conover(prog_df, val_col=prog_nam, group_col=categorical_key)
    elif test=='dscf':
        p_vals = posthoc_dscf(prog_df, val_col=prog_nam, group_col=categorical_key)

    for i, category in enumerate(mdata[prog_key].obs[categorical_key].astype(str).unique()):
        min_pval = np.min(p_vals.loc[p_vals.index!=category, category])
        mean_pval = np.mean(p_vals.loc[p_vals.index!=category, category])

        prog_idx = mdata[prog_key].var.index.get_loc(prog_nam)
        mdata[prog_key].varm['{}_association_min_pval'.format(categorical_key)][prog_idx,i] = min_pval
        mdata[prog_key].varm['{}_association_mean_pval'.format(categorical_key)][prog_idx,i] = mean_pval

    mdata[prog_key].uns['{}_association_pvals'.format(categorical_key)][prog_nam] = p_vals

# TODO: Add one way ANOVA, MANOVA, multi-variate KW
def compute_categorical_association(mdata, categorical_key='batch', n_jobs=1,
	                                prog_key='prog', inplace=True, **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata
    
    mdata[prog_key].var['{}_kruskall_wallis_stat'.format(categorical_key)] = None
    mdata[prog_key].var['{}_kruskall_wallis_pval'.format(categorical_key)] = None
    
    # Run in parallel (max n_jobs=num_progs)
    Parallel(n_jobs=n_jobs, backend='threading')(delayed(perform_kruskall_wallis)(mdata, 
                                                             prog_key=prog_key,
                                                             prog_nam=prog_nam, 
                                                             categorical_key=categorical_key) \
                                                             for prog_nam in tqdm(mdata[prog_key].var_names,
                                                             desc='Testing {} association'.format(categorical_key), 
                                                             unit='programs'))
    
    mdata[prog_key].varm['{}_association_min_pval'.format(categorical_key)] = np.zeros((mdata[prog_key].shape[1],
                                                               mdata[prog_key].obs[categorical_key].unique().shape[0]))
    mdata[prog_key].varm['{}_association_mean_pval'.format(categorical_key)] = np.ones((mdata[prog_key].shape[1],
                                                              mdata[prog_key].obs[categorical_key].unique().shape[0]))
    mdata[prog_key].uns['{}_association_pvals'.format(categorical_key)] = {}
    
    Parallel(n_jobs=n_jobs, backend='threading')(delayed(perform_posthoc)(mdata, 
                                                                          prog_key=prog_key,
                                                                          prog_nam=prog_nam, 
                                                                          categorical_key=categorical_key,
                                                                          **kwargs) \
                                                             for prog_nam in tqdm(mdata[prog_key].var_names,
                                                             desc='Identifying differential {}'.format(categorical_key), 
                                                             unit='programs'))
    
    if not inplace: return (mdata[prog_key].var.loc[:, ['{}_kruskall_wallis_stat'.format(categorical_key), 
                                                        '{}_kruskall_wallis_pval'.format(categorical_key)]],
                           mdata[prog_key].uns['{}_association_categories'.format(categorical_key)],
                           mdata[prog_key].varm['{}_association_stat'.format(categorical_key)], 
                           mdata[prog_key].varm['{}_association_pval'.format(categorical_key)],
                           mdata[prog_key].uns['{}_association_pvals'.format(categorical_key)]
                           )
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('categorical_key') 
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    compute_categorical_association(mdata, categorical_key=args.categorical_key, n_jobs=args.n_jobs,
                                    prog_key=args.prog_key, inplace=args.output)




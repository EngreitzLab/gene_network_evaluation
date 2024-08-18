import os
import argparse

import numpy as np
import pandas as pd

from scipy import stats

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import mudata

import logging
logging.basicConfig(level = logging.INFO)

# Convenience function for parallelization
# TODO: Make this usable as a standalone function with mdata as input
# TODO: Update mdata keys and use those to generate not inplace output
def compute_perturbation_association_(test_data, reference_data, program, level_name, test_stats_df):

    test_data_ = test_data[:, program].X.toarray()
    # TODO: Resample reference data to match number of obs
    reference_data_ = reference_data[:, program].X.toarray()

    results = stats.mannwhitneyu(test_data_, reference_data_)

    test_stats_df.append([level_name, program, results[0][0], results[1][0]])

# TODO: Add support for stratification by categorical levels
# Rename function to compute_{eval_measure}
def compute_perturbation_association(mdata, prog_key='prog',
                                     collapse_targets=True, 
                                     pseudobulk=False,
                                     reference_targets=('non-targeting'),
                                     n_jobs=1, inplace=True):
    """
    Computes eval_measure.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the gene program anndata object (mdata[prog_key]) in the mudata object.
        collapse_targets: bool (default:True)
            If target gene per guide is provided, perform tests on target levels. 
            Exclsuive with pseudobulk.
        pseudobulk: bool (default: False)
            If multiple non-targeting guides are available - optionally test at pseudobulk level.
            Exclusive with collapse_targets.
        reference_targets: tuple
            List of target values to use as reference distribution.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
    
        --- 
        if inplace:
            UPDATES
            mdata[prog_key].varm['perturbation_association_pval']
            mdata[prog_key].varm['perturbation_association_stat'] 
        else:
            RETURNS
            returns test_stats_df
            
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

    # Guide metadata
    guide_metadata = pd.DataFrame(index=mdata[prog_key].uns['guide_names'], columns=['Target'])
    guide_metadata['Target'] = mdata[prog_key].uns['guide_targets']

    # Run tests on pseudobulk if multiple non-targeting guides otherwise at single-cell level
    if collapse_targets:
        if pseudobulk:
            pseudobulk=False
            logging.info('Setting pseudobulk to False since collapse_targets is True')

    if pseudobulk:
        pseudobulked_scores = pd.DataFrame(index=mdata[prog_key].uns['guide_names'], 
                                           columns=mdata[prog_key].obs_names)

        #TODO: Parallelize
        for i, guide in enumerate(mdata[prog_key].uns['guide_names']):
            for program in mdata[prog_key].obs_names:
                pseudobulked_scores.loc[guide, program] = \
                mdata[prog_key][mdata[prog_key].obsm['guide_assignment'][:,i].astype(bool), program].X.mean()
        
        pseudobulked_scores['Target'] = mdata[prog_key].uns['guide_targets']

    # Create reference data
    if type(reference_targets) is str:
        reference_targets = [reference_targets]
    if pseudobulk:
        reference_data = pseudobulked_scores.loc[pseudobulked_scores.Target.isin(reference_targets)]
    else:
        reference_guide_idx = guide_metadata.index.get_indexer(guide_metadata.loc[guide_metadata.Target.isin(reference_targets)].index.values)
        reference_data = mdata[prog_key][np.any(mdata[prog_key].obsm['guide_assignment'][:,reference_guide_idx].astype(bool), axis=1)]

    # Run tests
    test_guides = guide_metadata.loc[~guide_metadata.Target.isin(reference_targets)].index.values

    test_stats_df = []
    if collapse_targets:
        level_key='target'

        # Run tests for each guide target
        for target in tqdm(guide_metadata['Target'].unique(), unit='targets'):
            test_guide_indices = guide_metadata.index.get_indexer(guide_metadata.index[guide_metadata['Target']==target])
            test_data = mdata[prog_key][mdata[prog_key].obsm['guide_assignment'][:,test_guide_indices].any(-1).astype(bool)]

            Parallel(n_jobs=n_jobs, backend='threading')(delayed(compute_perturbation_association_)(test_data, reference_data, program, target, test_stats_df) \
                                                        for program in mdata[prog_key].var_names)

    else:
        level_key='guide'

        # Run tests for each guide
        for guide in tqdm(test_guides, desc='Testing perturbation association', unit='guides'):
            # Run test at using pseudobulks at guide level -> requires multiple non-targeting guides as reference
            if pseudobulk:
                #TODO: Implement per guide probability under reference approach
                raise NotImplementedError()
            else:
                test_guide_idx = guide_metadata.index.get_loc(guide)
                test_data = mdata[prog_key][mdata[prog_key].obsm['guide_assignment'][:,test_guide_idx].astype(bool)]

                Parallel(n_jobs=n_jobs, backend='threading')(delayed(compute_perturbation_association_)(test_data, reference_data, program, guide, test_stats_df) \
                                                            for program in mdata[prog_key].var_names)

    test_stats_df = pd.DataFrame(test_stats_df, columns=['{}_name'.format(level_key), 'program_name', 'stat', 'pval'])

    # Return only the evaluations if not in inplace mode
    if not inplace: return test_stats_df
    else:
        init_array = np.empty(mdata[prog_key].shape[1], len(test_stats_df['{}_name'.format(level_key)].unique()))
        init_array[:] = np.nan

        stats, pvals = init_array.copy(), init_array.copy()

        for level_idx, level_name in enumerate(test_stats_df['{}_name'.format(level_key)].unique()):
            stats[:, level_idx] = test_stats_df.loc[test_stats_df['{}_name'.format(level_key)]==level_name, 'stat'][mdata[prog_key].var_names].values
            pvals[:, level_idx] = test_stats_df.loc[test_stats_df['{}_name'.format(level_key)]==level_name, 'pval'][mdata[prog_key].var_names].values
        mdata[prog_key].varm['perturbation_association_{}_stat'.format(level_key)] = stats
        mdata[prog_key].varm['perturbation_association_{}_pval'.format(level_key)] = pvals
        mdata[prog_key].uns['perturbation_association_{}_names'.format(level_key)] = test_stats_df['{}_name'.format(level_key)].unique().values
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-rt', '--reference_targets', nargs='+')
    parser.add_argument('--pseudobulk', default=False, action='store_true') 
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_perturbation_association(mdata, prog_key=args.prog_key, pseudobulk=args.pseudobulk, 
                                     reference_targets=args.reference_targets, n_jobs=args.n_jobs, 
                                     inplace=args.output)


                                     
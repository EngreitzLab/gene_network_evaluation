import os
import argparse

import numpy as np
import pandas as pd

from scipy import stats
from statsmodels.stats.multitest import multipletests

from joblib import Parallel, delayed
from typing import List, Dict, Tuple, Union, Optional, Literal
from tqdm.auto import tqdm

import mudata

import logging
logging.basicConfig(level = logging.INFO)


def get_guide_metadata(
    mdata: mudata.MuData,
    prog_key: str,
    guide_names_key: str = 'guide_names',
    guide_targets_key: str = 'guide_targets'
) -> pd.DataFrame:
    """
    Get guide metadata from the mudata object.

    Parameters
    ----------
    mdata : mudata.MuData
        The mudata object.
    prog_key : str
        The key of the program anndata in the mudata object.

    Returns
    -------
    pd.DataFrame
        The guide metadata.
    """
    guide_metadata = pd.DataFrame(index=mdata[prog_key].uns[guide_names_key], columns=['Target'])
    guide_metadata['Target'] = mdata[prog_key].uns[guide_targets_key]

    return guide_metadata


# Convenience function for parallelization
# TODO: Make this usable as a standalone function with mdata as input
# TODO: Update mdata keys and use those to generate not inplace output
def compute_perturbation_association_(
    test_data: mudata.MuData,
    reference_data: mudata.MuData,
    program: str,
    level_name: str,
    test_stats_df: List[List],
    balanced: bool = False
):
    # Get perturbation data
    test_data_ = test_data[:, program].X.toarray()

    # TODO: Resample reference data to match number of obs
    reference_data_ = reference_data[:, program].X.toarray()
    if balanced:
        reference_data_ = reference_data_[np.random.choice(reference_data_.shape[0], test_data_.shape[0], replace=False)]

    # Calculate log2FC
    ref_mean = np.mean(reference_data_)
    test_mean = np.mean(test_data_)
    log2fc = np.log2(test_mean / ref_mean)

    # Compute Mann-Whitney U test
    results = stats.mannwhitneyu(test_data_, reference_data_)

    # Append to test stats df
    test_stats_df.append([level_name, program, ref_mean, test_mean, log2fc, results[0][0], results[1][0]])


# TODO: Add support for stratification by categorical levels
def compute_perturbation_association(
    mdata: Union[os.PathLike, mudata.MuData],
    prog_key: str,
    guide_names_key: str = 'guide_names',
    guide_targets_key: str = 'guide_targets',
    guide_assignments_key: str = 'guide_assignment',
    collapse_targets: bool = True,
    pseudobulk: bool = False,
    reference_targets: Union[str, List[str]] = ['non-targeting'],
    balanced: bool = False,
    n_jobs: int = 1,
    inplace: bool = True
):
    """Compute Mann-Whitney U test for cells targeted with a perturbation against a reference set of cells.

    Parameters
    ----------
    mdata : MuData
        mudata object containing anndata of program scores and cell-level metadata.
    prog_key : str
        index for the gene program anndata object (mdata[prog_key]) in the mudata object.
    guide_names_key : str (default: 'guide_names')
        key in mdata[prog_key].uns for 1D np.array of guide names.
    guide_targets_key : str (default: 'guide_targets')
        key in mdata[prog_key].uns for 1D np.array of target gene names. Should be in the same order as guide_names.
    guide_assignments_key : str (default: 'guide_assignment')
        key in mdata[prog_key].obsm for 2D cell x guide assignment matrix. Should correspond to guide_names and guide_targets.
    collapse_targets : bool (default:True)
        If target gene per guide is provided, perform tests on target levels. 
        Mutually exclusive with pseudobulk.
    pseudobulk : bool (default: False)
        If multiple non-targeting guides are available - optionally test at pseudobulk level.
        Mutually exclusive with collapse_targets.
    reference_targets : tuple
        List of target values to use as reference distribution.
    balanced : bool (default: False)
        If True, resample reference data to match number of observations in test
    n_jobs: int (default: 1)
        number of threads to run processes on.
    inplace: Bool (default: True)
        update the mudata object inplace or return a copy

    Returns
    -------
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
        mdata = mudata.MuData({
            prog_key: mdata[prog_key].copy()
        })

    # Guide metadata
    guide_metadata = get_guide_metadata(
        mdata, 
        prog_key=prog_key,
        guide_names_key=guide_names_key, 
        guide_targets_key=guide_targets_key
    )

    # Run tests on pseudobulk if multiple non-targeting guides otherwise at single-cell level
    if collapse_targets:
        if pseudobulk:
            pseudobulk=False
            logging.info('Setting pseudobulk to False since collapse_targets is True')

    if pseudobulk:
        pseudobulked_scores = pd.DataFrame(index=mdata[prog_key].uns[guide_names_key], columns=mdata[prog_key].obs_names)

        #TODO: Parallelize
        for i, guide in enumerate(mdata[prog_key].uns['guide_names']):
            for program in mdata[prog_key].obs_names:
                pseudobulked_scores.loc[guide, program] = \
                mdata[prog_key][mdata[prog_key].obsm[guide_assignments_key][:,i].astype(bool), program].X.mean()
        
        pseudobulked_scores['Target'] = mdata[prog_key].uns['guide_targets']

    # Create reference data
    if type(reference_targets) is str:
        reference_targets = [reference_targets]

    # If pseudobulk, grab a DataFrame of all non-targeting guides
    if pseudobulk:
        reference_data = pseudobulked_scores.loc[pseudobulked_scores.Target.isin(reference_targets)]

    # Otherwise, grab a MuData for all those guides assigned to the reference targets
    else:
        reference_guide_idx = guide_metadata.index.get_indexer(guide_metadata.loc[guide_metadata.Target.isin(reference_targets)].index.values)
        reference_data = mdata[prog_key][np.any(mdata[prog_key].obsm[guide_assignments_key][:,reference_guide_idx].astype(bool), axis=1)]

    # Get the rest of the guides to test against the reference
    test_guides = guide_metadata.loc[~guide_metadata.Target.isin(reference_targets)].index.values

    # Run tests
    test_stats_df = []

    # If collapse_targets, run tests for each target
    if collapse_targets:

        # Set level key
        level_key = 'target'

        # Run tests for each target (with multiple guides)
        targets_no_ref = [targ for targ in guide_metadata['Target'].unique() if targ not in reference_targets]
        for target in tqdm(targets_no_ref, desc='Testing perturbation association', unit='targets'):

            # Grab data for the current target
            test_guide_idx = guide_metadata.index.get_indexer(guide_metadata.index[guide_metadata['Target']==target])
            test_data = mdata[prog_key][mdata[prog_key].obsm[guide_assignments_key][:,test_guide_idx].any(-1).astype(bool)]

            # Run test for every program
            Parallel(n_jobs=n_jobs, backend='threading')(delayed(compute_perturbation_association_)(
                test_data, 
                reference_data, 
                program, 
                target, 
                test_stats_df,
                balanced) for program in mdata[prog_key].var_names)

    else:

        # Set level key
        level_key='guide'

        # Run tests for each guide
        for guide in tqdm(test_guides, desc='Testing perturbation association', unit='guides'):
            
            # Run test at using pseudobulks at guide level -> requires multiple non-targeting guides as reference
            if pseudobulk:
                # TODO: Implement per guide probability under reference approach
                raise NotImplementedError()
            else:
                # Grab data for the current guide
                test_guide_idx = guide_metadata.index.get_loc(guide)
                test_data = mdata[prog_key][mdata[prog_key].obsm[guide_assignments_key][:,test_guide_idx].astype(bool)]

                # Run test for every program
                Parallel(n_jobs=n_jobs, backend='threading')(delayed(compute_perturbation_association_)(
                    test_data, 
                    reference_data, 
                    program, 
                    guide, 
                    test_stats_df,
                    balanced) for program in mdata[prog_key].var_names)
    
    # Create DataFrame
    test_stats_df = pd.DataFrame(test_stats_df, columns=['{}_name'.format(level_key), 'program_name', 'ref_mean', 'test_mean', 'log2FC', 'stat', 'pval'])

    # Correct for multiple testing
    test_stats_df['adj_pval'] = multipletests(test_stats_df['pval'], method='fdr_bh')[1]

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

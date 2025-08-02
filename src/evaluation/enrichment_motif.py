import os
import argparse

import mudata
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.multitest import fdrcorrection

from tangermeme.io import extract_loci
from memelite.io import read_meme
from memelite import fimo

from .enrichment_geneset import get_idconversion

from joblib import Parallel, delayed
from tqdm.auto import tqdm
from typing import List, Dict, Tuple, Union, Optional, Literal, Mapping


def read_loci(
    path_loci: os.PathLike
) -> pd.DataFrame:
    """Read loci to gene links file

    Read promoter/enhancer to gene links tab delimited file with the following column headers:
        - chromosome: chromosome of the locus
        - start: start position of the locus
        - end: end position of the locus
        - seq_name: name of the locus
        - seq_class {promoter, enhancer}: class of the locus
        - seq_score: score of the locus
        - gene_name: name of the gene linked to the locus

    Parameters
    ----------
    path_loci : str
        Path to the coordinates file.

    Returns
    -------
    loci : pd.DataFrame
        DataFrame containing the coordinates.
    """

    # Read formatted coords file
    loci = pd.read_csv(path_loci, sep='\t')

    # Enforce header format
    expected_headers = [
        'chromosome', 
        'start', 
        'end', 
        'seq_name', 
        'seq_class', 
        'seq_score', 
        'gene_name'
    ]

    # Check for expected headers
    for col in expected_headers:
        try: assert col in loci.columns
        except: raise ValueError('Coordinate file is not formatted correctly')

    # TODO: Support multiple classes at once
    try: assert len(loci['seq_class']==1)
    except: raise ValueError('Coordinate file contains multiple sequence classes')

    return loci

def perform_motif_match(
    loci: pd.DataFrame,
    sequences: os.PathLike,
    pwms: Mapping[str, np.ndarray],
    in_window: int=1000,
    sig: float=1e-4,
    eps: float=1e-4,
    reverse_complement: bool=True,
    output_loc: os.PathLike=None
):
        """Score motif matches

        Perform motif matching on sequences linked to genes
        via enhancer/promoter coordinates.

        Parameters
        ----------
        loci : pd.DataFrame
            DataFrame containing sequence coordinates.
        sequences : os.PathLike
            Path to FASTA formatted genomic sequence.
        pwms : Mapping[str, np.ndarray]
            Dictionary of PWMs where keys are motif names and values are PWMs.
        in_window : int
            Window size to extract sequences around center of loci.
        sig : float
            Threshold for motif matching.
        output_loc : os.PathLike
            Path to directory to store motif matches for individual motifs.
        """
        # Perform motif matching
        X = extract_loci(loci, sequences, in_window=in_window).float()
        hits = fimo(pwms, X, threshold=sig, eps=eps, reverse_complement=reverse_complement)

        # Create motif match dataframe
        motif_match_df = pd.DataFrame()
        for i, hit in enumerate(hits):
            annotated_hit = hit.merge(loci[["chromosome", "seq_name", "seq_class", "gene_name"]], left_on="sequence_name", right_index=True).drop(columns=["sequence_name"])
            annotated_hit["motif_name"] = list(pwms.keys())[i]
            annotated_hit = annotated_hit[["chromosome", "start", "end", "strand", "motif_name", "score", "p-value", "seq_name", "seq_class", "gene_name"]]
            if output_loc is not None:
                annotated_hit.to_csv(os.path.join(output_loc, f"motif_match_{list(pwms.keys())[i]}.txt"), sep='\t', index=False)
            motif_match_df = pd.concat([motif_match_df, annotated_hit])
        motif_match_df["adj_pval"] = multipletests(motif_match_df["p-value"], method="fdr_bh")[1]
        motif_match_df.reset_index(drop=True, inplace=True)
        return motif_match_df


def compute_motif_instances(
    motif_match_df: pd.DataFrame,
    motif_var: str = 'motif_name',
    gene_names: Optional[np.ndarray]=None
):
    """Count motif instances per gene (via enahncer/promoter linking)

    Parameters
    ----------
    motif_match_df : pd.DataFrame
        DataFrame containing motif matches.
    motif_var : str
        Column name for motif names. Default is 'motif_name'.
    gene_names : np.ndarray
        Array of gene names. Default is None.
    
    Returns
    -------
    motif_count_df : pd.DataFrame
        DataFrame containing
    """

    # Count up significant occurences of motif
    motif_match_df_ = motif_match_df.value_counts(subset=['gene_name', motif_var]).reset_index()
    motif_match_df_.columns = ['gene_name', motif_var, 'motif_count']
    motif_match_df_ = motif_match_df_.pivot(index='gene_name', columns=motif_var, values='motif_count')

    motif_count_df = pd.DataFrame(index=gene_names, columns=motif_match_df_.columns)
    motif_count_df.loc[motif_match_df_.index.values] = motif_match_df_ # Gene names should match as this point

    return motif_count_df


def perform_correlation(
    motif_count_df: pd.DataFrame,
    prog_genes: pd.DataFrame,
    motif_enrich_stat_df: pd.DataFrame,
    motif_enrich_pval_df: pd.DataFrame,
    motif_idx: int,
    prog_idx: int,
    correlation: str='pearsonr',
    low_cutoff: float = -np.inf,
    n_top: int = None,
):
    """Compute motif enrichment as correlation b/w gene weights/ranks and motif counts
    
    Perform pearson correlation test for motif count enrichment vs gene loadings
    If loadings are dichotomized then this is equivalent to a point biserial correlation test.

    Parameters
    ----------
    motif_count_df : pd.DataFrame
        DataFrame containing motif counts.
    prog_genes : pd.DataFrame
        DataFrame containing gene program loadings.
    motif_enrich_stat_df : pd.DataFrame
        DataFrame to store correlation statistics.
    motif_enrich_pval_df : pd.DataFrame
        DataFrame to store correlation p-values.
    motif_idx : int
        Index of motif.
    prog_idx : int
        Index of gene program.
    correlation : {'pearsonr','spearmanr','kendalltau'}
        Type of correlation to perform. Default is 'pearsonr'.
    low_cutoff : float
        Remove features with loadings at or below this value.
    n_top : int
        Take the top n features with the highest loadings.
    """

    # Make unique index for subsetting
    prog_genes_ = prog_genes.copy().T
    prog_genes_['unique_id'] = prog_genes_.index.astype(str) + '_' + \
                               prog_genes_.groupby(prog_genes_.index).cumcount().astype(str)

    motif_count_df_ = motif_count_df.copy()
    motif_count_df_.set_index(prog_genes_['unique_id'], inplace=True, drop=True)

    prog_genes_.set_index('unique_id', inplace=True, drop=True)

    loadings = prog_genes_.iloc[:, prog_idx]

    # Filter out low loadings
    loadings = loadings[(loadings > low_cutoff)]

    # Take top n features if specified
    if n_top is not None:
        loadings = loadings.sort_values(ascending=False).head(n_top)
        if len(loadings) < n_top:
            logging.warning(f"Program {i} has less than {n_top} features after filtering. Only {len(loadings)} features will be used.")

    counts = motif_count_df_.T.loc[:, loadings.index]

    loadings = loadings.values.flatten()
    counts = counts.iloc[motif_idx].fillna(0).values.flatten()
    
    if correlation=='pearsonr':
        stat, pval = pearsonr(loadings, counts)
    elif correlation=='spearmanr':
        stat, pval = spearmanr(loadings, counts)
    elif correlation=='kendalltau':
        stat, pval = kendalltau(loadings, counts)

    motif_enrich_stat_df.iloc[prog_idx, motif_idx]  = stat
    motif_enrich_pval_df.iloc[prog_idx, motif_idx]  = pval

def compute_motif_enrichment_(
    mdata: mudata.MuData,
    motif_count_df: pd.DataFrame,
    prog_key: str='prog',
    gene_names: Optional[np.ndarray]=None,
    weighted: bool=True,
    num_genes: Optional[int]=None,
    correlation: str='pearsonr',
    low_cutoff: float = -np.inf,
    n_top: int = None,
    n_jobs: int=1
):
    """Count up motif ocurrences and perform diff. test

    Perform motif enrichment using gene program loadings and
    motif counts linked to genes via motif scanning of
    linked enhancer/promoter sequences.
    
    Parameters
    ----------
    mdata : MuData
        MuData object containing anndata of program scores and cell-level metadata.
    motif_count_df : pd.DataFrame
        DataFrame containing motif counts.
    prog_key : str
        Key for the anndata object in the mudata object. Default is 'prog'.
    gene_names : np.ndarray
        Array of gene names. Default is None.
    weighted : bool
        Use weighted loadings. Default is True.
    num_genes : int
        Number of genes threshold to dichotomize loadings. Default is None.
    correlation : {'pearsonr','spearmanr','kendalltau'}
        Type of correlation to perform. Default is 'pearsonr'.
    low_cutoff : float
        Remove features with loadings at or below this value.
    n_top : int
        Take the top n features with the highest loadings.
    n_jobs : int
        Number of threads to run processes on. Default is 1.
    """

    # Both weighted and num_genes cannot be set
    if num_genes is not None and weighted:
        raise ValueError('Will not use weighted when num_genes specified.')
    
    loadings = pd.DataFrame(data=mdata[prog_key].varm['loadings'],
                            index=mdata[prog_key].var_names,
                            columns=gene_names)

    # FIXME: Causes expansion due to duplications in gene_names
    # Ensure index matches b/w loadings and counts
    # loadings = loadings.loc[:, motif_count_df.index.values]

    # Binary matrix 
    if not weighted:
        prog_genes = (loadings.rank(axis=1)<num_genes).astype(int)
    elif weighted and num_genes is None:
        prog_genes = loadings
    elif weighted and isinstance(num_genes, int):
        prog_genes = (loadings.rank(axis=1)<num_genes).astype(int)
        prog_genes *= loadings
    else:
        raise ValueError('num_genes not specified correctly.')

    # Use pearson correlation
    # If dichotomized then this is equivalent to point biserial correlation
    motif_enrich_stat_df = pd.DataFrame(index=mdata[prog_key].var_names,
                                        columns=motif_count_df.columns.values)
    motif_enrich_pval_df = pd.DataFrame(index=mdata[prog_key].var_names,
                                        columns=motif_count_df.columns.values)

    # FIXME: If n_jobs>1 then parallel processes dont seem to terminate.
    # Perform test in parallel across motifs and programs
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(perform_correlation)(motif_count_df,
                                                               prog_genes, 
                                                               motif_enrich_stat_df,
                                                               motif_enrich_pval_df,
                                                               motif_idx, 
                                                               prog_idx,
                                                               correlation=correlation) \
                                                            for motif_idx in tqdm(range(motif_count_df.columns.values.shape[-1]),
                                                                                     desc='Computing motif enrichment',
                                                                                     unit='motifs') \
                                                            for prog_idx in range(mdata[prog_key].var_names.shape[0]))

    return motif_enrich_stat_df, motif_enrich_pval_df


def compute_motif_enrichment(
    mdata: Union[mudata.MuData, os.PathLike],
    prog_key: str='prog',
    data_key: str='data',
    organism: Literal['human', 'mouse'] = 'human',
    motif_file: Optional[os.PathLike]=None,
    seq_file: Optional[os.PathLike]=None,
    loci_file: Optional[os.PathLike]=None,
    output_loc: Optional[os.PathLike]=None,
    window: int=1000,
    sig: float=1e-4,
    eps: float=1e-4,
    reverse_complement: bool=True,
    num_genes: Optional[int]=None,
    correlation: str='pearsonr',
    low_cutoff: float = -np.inf,
    n_top: int = 2000,
    n_jobs: int=1,
    inplace: bool=True,
    **kwargs
):
    
    """Compute motif enrichment in enhancers or promoters associated with a gene
    
    Perform motif enrichment using gene program loadings and
    motif counts linked to genes via motif scanning of
    linked enhancer/promoter sequences.

    Parameters
    ----------
    mdata : MuData
        mudata object containing anndata of program scores and cell-level metadata.
    prog_key: 
        index for the anndata object (mdata[prog_key]) in the mudata object.
    data_key: str
        index of the genomic data anndata object (mdata[data_key]) in the mudata object.
    motif_file: str
        path to motif file formatted in MEME format.
    seq_file: str
        path to FASTA formatted genomic sequence.
    loci_file: str
        path to enhancer/promoter gene links file with sequence coordinates.
        Tab delimited file with the following column headers:
        chr, start, end, seq_name, seq_class {promoter, enhancer}, 
        seq_score, gene_name.
    organism : {'human', 'mouse'} (default: 'human')
        species to which the sequencing data was aligned to.
    output_loc: str
        path to directory to store motif - gene counts.
    sig: (0,1] (default: 0.05)
        significance level for inferring a motif match.
    num_genes: int (default: None)
        number of genes threshold to dichtomize loadings.
    correlation: {'pearsonr','spearmanr','kendalltau'} (default: 'peasronsr')
        correlation type to use to compute motif enirchments.
        Use kendalltau when expecting enrichment/de-enrichment at both ends.
    low_cutoff : float
        Remove features with loadings at or below this value.
    n_top : int
        Take the top n features with the highest loadings.
    use_previous: bool (default: True)
        if outplot is provided try to load motif matches from previous run.
    n_jobs: int (default: 1)
        number of threads to run processes on.
    inplace: Bool (default: True)
        update the mudata object inplace or return a copy
       
    Returns
    -------
    if not inplace:
        motif_match_df,
        motif_count_df.loc[gene_names].values,
        motif_enrichment_df
    else:
        None, edits mdata in place
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
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})

    # Get gene names in MuData
    if 'var_names' in mdata[prog_key].uns.keys():
        gene_names = get_idconversion(mdata[prog_key].uns['var_names'], organism=organism)
    else:
        try: assert mdata[prog_key].varm['loadings'].shape[1] == mdata[data_key].var.shape[0]
        except: raise ValueError('Different number of genes present in data and program loadings')
        gene_names = get_idconversion(mdata[data_key].var_names, organism=organism)

    # # 
    # if ':ens' in gene_names[0].lower():
    #     gene_names = [name.split(':')[0] for name in gene_names]

    # Check if output loc exists
    if output_loc is not None:
        try: os.makedirs(output_loc, exist_ok=True)
        except: raise ValueError('Output location does not exist.')

    # If num_genes specified then cannot be weighted
    if num_genes is None:
        weighted=True
    elif isinstance(num_genes, int):
        weighted=False

    # Intake motif file path or in memory
    if os.path.exists(motif_file):
        pwms = read_meme(motif_file)
    else:
        raise ValueError('Motif file not found.')
    
    # Intake coord file path or in memory
    if os.path.exists(loci_file):
        loci = read_loci(loci_file)
    else:
        raise ValueError('Coordinate file not found.')

    # Valid genes
    matching_gene_names = np.intersect1d(loci['gene_name'].unique(), gene_names)
    print(f'Number of matching genes: {len(matching_gene_names)}')
    try: assert len(matching_gene_names) > 0
    except: raise ValueError('No matching genes b/w data and coordinate files')

    # Compute motif matching
    loci_ = loci[loci['gene_name'].isin(matching_gene_names)]
    print(f'Number of loci: {loci_.shape[0]}')

    motif_match_df = perform_motif_match(
        loci=loci_,
        sequences=seq_file,
        pwms=pwms,
        in_window=window,
        sig=sig,
        eps=eps,
        reverse_complement=reverse_complement,
        output_loc=output_loc
    )

    # Count motif enrichment
    motif_count_df = compute_motif_instances(
        motif_match_df,
        motif_var='motif_name',
        gene_names=gene_names
    )
    motif_enrich_stat_df, motif_enrich_pval_df = compute_motif_enrichment_(
        mdata,
        motif_count_df,
        prog_key=prog_key,
        gene_names=gene_names,
        weighted=weighted,
        num_genes=num_genes,
        n_jobs=n_jobs)
                         
    # Store motif counts
    if inplace:
        mdata[prog_key].uns['motif_counts'] = motif_count_df.loc[gene_names].values
        mdata[prog_key].uns['motif_names'] = motif_count_df.columns.values

        mdata[prog_key].varm['motif_enrich_{}_stat'.format(correlation)] = motif_enrich_stat_df.values
        mdata[prog_key].varm['motif_enrich_{}_pval'.format(correlation)] = motif_enrich_pval_df.values
        mdata[prog_key].uns['motif_names'] = motif_count_df.columns.values

    else:

        motif_enrich_stat_df = motif_enrich_stat_df.reset_index().melt(id_vars='index',
                                                                       var_name='motif', 
                                                                       value_name='stat')
        motif_enrich_stat_df = motif_enrich_stat_df.set_index(['index', 'motif'])

        motif_enrich_pval_df = motif_enrich_pval_df.reset_index().melt(id_vars='index',
                                                                       var_name='motif', 
                                                                       value_name='pval')
        motif_enrich_pval_df = motif_enrich_pval_df.set_index(['index', 'motif'])

        motif_enrichment_df = motif_enrich_stat_df.merge(motif_enrich_pval_df,
                                                        left_index=True, 
                                                        right_index=True)
        motif_enrichment_df = motif_enrichment_df.reset_index()
        motif_enrichment_df['program_name'] = motif_enrichment_df['index']
        motif_enrichment_df.drop('index', axis=1, inplace=True)
        motif_enrichment_df = motif_enrichment_df.sort_values(['program_name', 'pval'])
        motif_enrichment_df["adj_pval"] = multipletests(motif_enrichment_df["pval"], method="fdr_bh")[1]

        motif_count_df.columns.name=''

        return (motif_match_df,
                motif_count_df.loc[gene_names].fillna(0),
                motif_enrichment_df)


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-mf', '--motif_file', default=None, type=str, required=True) 
    parser.add_argument('-sf', '--seq_file',  default=None, type=str, required=True) 
    parser.add_argument('-lf', '--loci_file', default=None, type=str, required=True) 
    parser.add_argument('--organism', default='human', type=str, choices=['human', 'mouse']) 
    parser.add_argument('--store_files', default=None, type=str)
    parser.add_argument('--significance', default=0.05, choices=range(0,1), 
                                          metavar="(0,1)", type=float)
    parser.add_argument('--num_genes', default=None, type=int)
    parser.add_argument('--correlation', default='pearsonr', type=str)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_motif_enrichment(mdata, prog_key=args.prog_key, data_key=args.data_key, 
                             motif_file=args.motif_file, seq_file=args.seq_file, 
                             loci_file=args.loci_file, organism=args.organism,
                             output_loc=args.store_files, sig=args.significance, 
                             num_genes=args.num_genes, n_jobs=args.n_jobs, 
                             correlation=args.correlation, inplace=args.output)





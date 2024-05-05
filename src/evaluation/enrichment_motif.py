import os
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymemesuite.common import MotifFile, Sequence, Background, Alphabet
from pymemesuite.fimo import FIMO

import mudata
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr, kendalltau

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# Compute background freq based on selected sequences
def compute_background_freq(selected_sequences):

    background = Background.from_sequences(Alphabet.dna(), 
                                           *selected_sequences)
    return background

# Read motif database 
def read_motif_file(motif_file_loc):

    """
    Read MEME formatted motif file.

    """

    # Motifs in meme format
    records = []
    with MotifFile(motif_file_loc) as motif_file:
        for motif in motif_file:
            records.append(motif)

    # If not in meme format then nothing is read
    try: assert len(records)>0
    except: raise ValueError('No motifs detected. Check if meme formatted.')

    return records

# Read sequence database
def read_sequence_file(seq_file_loc, format_='fasta'):

    """
    Read genomic sequence to scan for motifs.

    """

    # parse sequence file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(seq_file_loc), format_))

    return records

# Read coordinates - tab formatted file
def read_coords_file(coords_file_loc):

    """
    Read promoter/enhancer to gene links tab 
    delimited file with the following column headers:
    chr, start, end, seq_name, seq_class {promoter, enhancer}, 
    seq_score, gene_name

    """

    # Read formatted coords file
    records = pd.read_csv(coords_file_loc, sep='\t')

    # Enforced header format
    expected_headers = ['chr', 'start', 'end', 'seq_name', 
                        'seq_class', 'seq_score', 'gene_name']
    for col in expected_headers:
        try: assert col in records.columns
        except: raise ValueError('Coordinate file is not formatted correctly')

    # TODO: Support multiple classes at once
    try: assert len(records['seq_class']==1)
    except: raise ValueError('Coordinate file contains multiple sequence classes')

    return records

# Return relevant sequences
# TODO: Parallelize
def get_sequences(sequences_record, coords):

    """
    Extract sequences from sequence file to be scanned
    for motifs loaded from the motif file.

    """

    # search for short sequences
    selected_sequences = []
    for chr_ in coords['chr'].unique():
        coords_ = coords.loc[coords['chr']==chr_]
        long_seq_record = sequences_record[chr_]
        long_seq = long_seq_record.seq

        for idx in coords_.index.values:
            short_seq = str(long_seq)[coords_.loc[idx, 'start']-1:coords_.loc[idx, 'end']]
            short_seq_record = SeqRecord(Seq(short_seq), 
                                         id=coords_.loc[idx, 'seq_name'], 
                                         description='_'.join((coords_.loc[idx, 'seq_class'],
                                                              coords_.loc[idx, 'gene_name'])))
            selected_sequences.append(short_seq_record)

    # Convert to pymemesuite Sequence class
    selected_sequences = [Sequence(str(record.seq), name=record.id.encode()) \
                          for record in selected_sequences]

    return selected_sequences

# Score motif matches
def perform_motif_match(sequences_record, motifs_record, coords,
                        motif_match_dict, gene_seq_num, output_loc, 
                        motif_idx, gene_name):
    
    """
    Use FIMO to identify motif instances in sequences.
    
    """

    # Subselect coords to gene
    coords_ = coords.loc[coords['gene_name']==gene_name]

    # Get gene associated sequences
    gene_assoc_sequences = get_sequences(sequences_record, coords_)
    gene_seq_num[gene_name] = len(gene_assoc_sequences)

    # Compute background
    background = compute_background_freq(gene_assoc_sequences)

    # Compute motif matching
    motif = motifs_record[motif_idx]

    pattern_matches = {}
    fimo = FIMO(both_strands=True)
    pattern = fimo.score_motif(motif, gene_assoc_sequences, background)
    
    # Store in dataframe
    motif_match_data_ = [[m.source.accession.decode(), m.start, m.stop, m.strand,
                          m.score, m.pvalue, m.qvalue] for m in pattern.matched_elements]

    motif_match_df_ = pd.DataFrame(data=motif_match_data_, 
                                columns=['seq_name', 'start', 'end', 
                                         'strand', 'score', 'pvalue', 'qvalue'])

    motif_match_dict[(gene_name, motif.accession.decode())] = motif_match_df_

    if output_loc is not None and np.unique(motif_match_data_).shape[0]!=0:
        motif_match_df_.to_csv(os.path.join(output_loc, '_'.join((gene_name, 
                                                                 motif.accession.decode(),
                                                                 '.txt'))), sep='\t',
                                                                 index=False)

# Count up motif occurences per gene
def compute_motif_instances(mdata, motif_match_df, sig=0.05, gene_names=None):

    """
    Count motif instances per gene (via enahncer/promoter linking)

    """

    # Count up significant occurences of motif
    motif_match_df_ = motif_match_df.loc[motif_match_df.qvalue<=sig]
    motif_match_df_ = motif_match_df.value_counts(subset=['gene_name', 'motif_name']).reset_index()
    motif_match_df_ = motif_match_df_.pivot(index='gene_name', columns='motif_name', values='count')

    motif_count_df = pd.DataFrame(index=gene_names, columns=motif_match_df_.columns)
    motif_count_df.loc[motif_match_df_.index.values] = motif_match_df_ # Gene names should match as this point

    return motif_count_df

# Perform pearson correlation test for motif count enrichment vs gene loadings
# If loadings are dichotomized then this is equivalent to a point biserial correlation test.
def perform_correlation(motif_count_df, prog_genes, 
                        motif_enrich_stat_df, 
                        motif_enrich_pval_df, 
                        motif_idx, prog_idx,
                        correlation='pearsonr'):
    """
    Compute motif enrichment as correlation b/w gene weights/ranks and motif counts
    
    """

    loadings = prog_genes.iloc[prog_idx].values.flatten()
    counts = motif_count_df.T.iloc[motif_idx].fillna(0).values.flatten()
    
    if correlation=='pearsonr':
        stat, pval = pearsonr(loadings, counts)
    elif correlation=='spearmanr':
        stat, pval = spearmanr(loadings, counts)
    elif correlation=='kendalltau':
        stat, pval = kendalltau(loadings, counts)

    motif_enrich_stat_df.iloc[prog_idx, motif_idx]  = stat
    motif_enrich_pval_df.iloc[prog_idx, motif_idx]  = pval

# Count up motif ocurrences and perform diff. test
def compute_motif_enrichment_(mdata, motif_count_df, prog_key='prog',  
                              gene_names=None, weighted=True, num_genes=None, 
                              correlation='pearsonr', n_jobs=1):
    """
    Perform motif enrichment using gene program loadings and
    motif counts linked to genes via motif scanning of
    linked enhancer/promoter sequences.
    
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

# Compute motif enrichment in enhancers or promoters associated with a gene
def compute_motif_enrichment(mdata, prog_key='prog', data_key='rna', motif_file=None, 
                             seq_file=None, coords_file=None, output_loc=None, sig=0.05, 
                             num_genes=None, correlation='pearsonr', n_jobs=1, inplace=True, 
                             **kwargs):
    
    """
    Perform motif enrichment using gene program loadings and
    motif counts linked to genes via motif scanning of
    linked enhancer/promoter sequences.

    ARGS
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
        coords_file: str
            path to enhancer/promoter gene links file with sequence coordinates.
            Tab delimited file with the following column headers:
            chr, start, end, seq_name, seq_class {promoter, enhancer}, 
            seq_score, gene_name.
        output_loc: str
            path to directory to store motif - gene counts.
        sig: (0,1] (default: 0.05)
            significance level for inferring a motif match.
        num_genes: int (default: None)
            number of genes threshold to dichtomize loadings.
        correlation: {'pearsonr','spearmanr','kendalltau'} (default: 'peasronsr')
            correlation type to use to compute motif enirchments.
            Use kendalltau when expecting enrichment/de-enrichment at both ends.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
       
    RETURNS 
        if not inplace:
           mdata[prog_key].uns['motif_matching'],
           mdata[prog_key].uns['motif_counts'],
           mdata[prog_key].uns['motif_names'],
           mdata[prog_key].varm['motif_enrich_{}_stat'.format(correlation)],
           mdata[prog_key].varm['motif_enrich_{}_pval'.format(correlation)],
           mdata[prog_key].uns['motif_names']     

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

    if 'var_names' in mdata[prog_key].uns.keys():
        gene_names = mdata[prog_key].uns['var_names']
    else:
        try: assert mdata[prog_key].varm['loadings'].shape[1]==mdata[data_key].var.shape[0]
        except: raise ValueError('Different number of genes present in data and program loadings')
        gene_names = mdata[data_key].var_names
    
    # 
    if ':ens' in gene_names[0].lower():
        gene_names = [name.split(':')[0] for name in gene_names]

    # Check if output loc exists
    if output_loc is not None:
        try: os.makedirs(output_loc, exist_ok=True)
        except: raise ValueError('Output location does not exist.')

    # If num_genes specified then cannot be weighted
    if num_genes is None:
        weighted=True
    elif isinstance(num_genes, int):
        weighted=False

    # Set kwargs
    if 'sequence_kwargs' not in kwargs.keys():
        kwargs['sequence_kwargs'] = {}

    # Intake motif file path or in memory
    if isinstance(motif_file, str) and os.path.exists(motif_file):
        motifs_record = read_motif_file(motif_file)
    else:
        motifs_record = motif_file

    # Intake sequence file path or in memory
    if isinstance(seq_file, str) and os.path.exists(seq_file):
        sequences_record = read_sequence_file(seq_file, 
                                              **kwargs['sequence_kwargs'])
    else:
        sequences_record = seq_file
    
    # Intake coord file path or in memory
    if isinstance(coords_file, str) and os.path.exists(coords_file):
        coords = read_coords_file(coords_file)
    else:
        coords = coords_file

    # Valid genes
    matching_gene_names = np.intersect1d(coords['gene_name'].unique(), 
                                         gene_names)

    try: assert len(matching_gene_names) > 0
    except: raise ValueError('No matching genes b/w data and coordinate files')

    # Compute motif matching
    motif_match_dict = {}
    gene_seq_num = {}
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(perform_motif_match)(sequences_record, 
                                                               motifs_record, 
                                                               coords,
                                                               motif_match_dict,
                                                               gene_seq_num,
                                                               output_loc,
                                                               motif_idx, 
                                                               gene_name) \
                                                               for motif_idx in tqdm(range(len(motifs_record)),
                                                                                     desc='Matching motifs to sequences',
                                                                                     unit='motifs') \
                                                               for gene_name in tqdm(matching_gene_names,
                                                                                     desc='Motif scanning',
                                                                                     unit='genes'))

    # Concatenate motif matches 
    motif_match_dfs = [] 
    for key, motif_match_df_ in motif_match_dict.items():
        motif_match_df_['gene_name'] = key[0]
        motif_match_df_['motif_name'] = key[1]

        motif_match_dfs.append(motif_match_df_)

    motif_match_df = pd.concat(motif_match_dfs)
    mdata[prog_key].uns['motif_matching'] = motif_match_df

    # Count motif enrichment
    motif_count_df = compute_motif_instances(mdata, motif_match_df, sig=sig, gene_names=gene_names)
    motif_enrich_stat_df, motif_enrich_pval_df = compute_motif_enrichment_(mdata, motif_count_df, prog_key=prog_key, 
                                                                           gene_names=gene_names, weighted=weighted, 
                                                                           num_genes=num_genes, n_jobs=n_jobs)
                         
    # Store motif counts
    mdata[prog_key].uns['motif_counts'] = motif_count_df.loc[gene_names].values
    mdata[prog_key].uns['motif_names'] = motif_count_df.columns.values

    mdata[prog_key].varm['motif_enrich_{}_stat'.format(correlation)] = motif_enrich_stat_df.values
    mdata[prog_key].varm['motif_enrich_{}_pval'.format(correlation)] = motif_enrich_pval_df.values
    mdata[prog_key].uns['motif_names'] = motif_count_df.columns.values

    if not inplace: 

        motif_enrich_stat_df = motif_enrich_stat_df.melt(var_name='motif', value_name='stat')
        motif_enrich_stat_df = motif_enrich_stat_df.reset_index().set_index(['index', 'motif'])

        motif_enrich_pval_df = motif_enrich_pval_df.melt(var_name='motif', value_name='pval')
        motif_enrich_pval_df = motif_enrich_pval_df.reset_index().set_index(['index', 'motif'])

        motif_enrichment_df = motif_enrich_stat_df.merge(motif_enrich_pval_df,
                                                        left_index=True, 
                                                        right_index=True)
        motif_enrichment_df = motif_enrichment_df.reset_index()

        return (motif_match_df,
                motif_count_df.loc[gene_names].values,
                motif_enrichment_df)                
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-mf', '--motif_file', default=None, type=str, required=True) 
    parser.add_argument('-sf', '--seq_file',  default=None, type=str, required=True) 
    parser.add_argument('-cf', '--coords_file', default=None, type=str, required=True) 
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
                             motif_file=args.motif_file, seq_file=arg.seq_file, 
                             coords_file=arg.coords_file, output_loc=args.store_files, 
                             sig=args.significance, num_genes=args.num_genes, 
                             n_jobs=args.n_jobs, correlation=args.correlation, 
                             inplace=args.output)





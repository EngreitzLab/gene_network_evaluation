import os
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymemesuite.common import MotifFile, Sequence, Background, Alphabet
from pymemesuite.fimo import FIMO

import numpy as np
import pandas as pd

from scipy.stats import pearsonr

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import warnings

# Compute background freq based on selected sequences
def compute_background_freq(selected_sequences):

    background = Background.from_sequences(Alphabet.dna(), 
                                           *selected_sequences)
    return background

# Read motif database 
def read_motif_file(motif_file_loc):

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

    # parse sequence file and turn into dictionary
    records = SeqIO.to_dict(SeqIO.parse(open(seq_file_loc), format_))

    return records

# Read coordinates - tab formatted file
# chr, start, end, seq_name, seq_class {promoter, enhancer}, seq_score, gene_name
def read_coords_file(coords_file_loc):

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
def compute_motif_instances(mdata, motif_match_df, sig=0.05, data_key='rna'):

    # Count up significant occurences of motif
    motif_match_df_ = motif_match_df.loc[motif_match_df.qvalue<=sig]
    motif_match_df_ = motif_match_df.value_counts(subset=['gene_name', 'motif_name']).reset_index()
    motif_match_df_ = motif_match_df_.pivot(index='gene_name', columns='motif_name', values='count')

    motif_count_df = pd.DataFrame(index=mdata[data_key].var_names, columns=motif_match_df_.columns)
    motif_count_df.loc[motif_match_df_.index.values] = motif_match_df_ # Gene names should match as this point

    return motif_count_df

# Perform pearson correlation test for motif count enrichment vs gene loadings
# If loadings are dichotomized then this is equivalent to a point biserial correlation test.
def perform_pearsonr(motif_count_df, prog_genes, 
                     motif_enrich_stat_df, motif_enrich_pval_df, 
                     motif_idx, prog_idx):

    loadings = prog_genes.iloc[prog_idx].values.flatten()
    counts = motif_count_df.T.iloc[motif_idx].fillna(0).values.flatten()
    
    stat, pval = pearsonr(loadings, counts)

    motif_enrich_stat_df.iloc[prog_idx, motif_idx]  = stat
    motif_enrich_pval_df.iloc[prog_idx, motif_idx]  = pval

# Count up motif ocurrences and perform diff. test
def compute_motif_enrichment_(mdata, motif_count_df, prog_key='prog', data_key='rna', 
                              weighted=True, num_genes=None, n_jobs=-1):

    # Both weighted and num_genes cannot be set
    if num_genes is not None and weighted:
        raise ValueError('Will not use weighted when num_genes specified.')

    loadings = pd.DataFrame(data=mdata[prog_key].varm['loadings'],
                            index=mdata[prog_key].var_names,
                            columns=mdata[data_key].var_names)
    loadings = loadings.loc[:, motif_count_df.index.values]

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

    # Perform test in parallel across motifs and programs
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(perform_pearsonr)(motif_count_df,
                                                            prog_genes, 
                                                            motif_enrich_stat_df,
                                                            motif_enrich_pval_df,
                                                            motif_idx, 
                                                            prog_idx) \
                                                            for motif_idx in tqdm(range(motif_count_df.columns.values.shape[-1]),
                                                                                     desc='Computing motif enrichment',
                                                                                     unit='motifs') \
                                                            for prog_idx in tqdm(range(mdata[prog_key].var_names.shape[0]),
                                                                                     desc='per program',
                                                                                     unit='programs'))

    return motif_enrich_stat_df, motif_enrich_pval_df

# Compute motif enrichment in enhancers or promoters associated with a gene
def compute_motif_enrichment(mdata, motif_file=None, seq_file=None, coords_file=None, 
                             n_jobs=1, prog_key='prog', data_key='rna', output_loc=None,
                             sig=0.05, num_genes=None, inplace=True, **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

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
                                         mdata[data_key].var_names)

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
                                                                                     unit='motifs'     ) \
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
    mdata[data_key].uns['motif_matching'] = motif_match_df

    # Count motif enrichment
    motif_count_df = compute_motif_instances(mdata, motif_match_df, sig=sig, data_key=data_key)
    motif_enrich_stat_df, motif_enrich_pval_df = compute_motif_enrichment_(mdata, motif_count_df, prog_key=prog_key, 
                                                                           data_key=data_key, weighted=weighted, 
                                                                           num_genes=num_genes, n_jobs=n_jobs) 
                         
    # Store motif counts
    mdata[data_key].varm['motif_counts'] = motif_count_df.loc[mdata[data_key].var_names].values
    mdata[data_key].uns['motif_names'] = motif_count_df.columns.values

    mdata[prog_key].varm['motif_enrich_stat'] = motif_enrich_stat_df.values
    mdata[prog_key].varm['motif_enrich_pval'] = motif_enrich_pval_df.values
    mdata[prog_key].uns['motif_names'] = motif_count_df.columns.values

    return mdata[data_key].uns['motif_matching'], \
           mdata[data_key].varm['motif_counts'], \
           mdata[data_key].uns['motif_names'], \
           mdata[prog_key].varm['motif_enrich_stat'], \
           mdata[prog_key].varm['motif_enrich_pval']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj')
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str)
    parser.add_argument('-mf', '--motif_file', default=None, type=str, required=True) 
    parser.add_argument('-sf', '--seq_file',  default=None, type=str, required=True) 
    parser.add_argument('-cf', '--coords_file', default=None, type=str, required=True) 
    parser.add_argument('--store_files', default=None, type=str)
    parser.add_argument('--significance', default=0.05, choices=range(0,1), 
                                          metavar="(0,1)", type=float)
    parser.add_argument('--num_genes', default=None, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()
    
    compute_motif_enrichment(args.mudataObj, motif_file=args.motif_file, 
                             seq_file=arg.seq_file, coords_file=arg.coords_file, 
                             n_jobs=args.n_jobs, prog_key=args.prog_key, 
                             data_key=args.data_key, output_loc=args.output,
                             sig=args.significance, num_genes=args.num_genes, 
                             inplace=args.output)





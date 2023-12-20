import os
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymemesuite.common import MotifFile, Sequence, Background, Alphabet
from pymemesuite.fimo import FIMO

import numpy as np
import pandas as pd

from joblib import Parallel, delayed
from tqdm.auto import tqdm

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
                        motif_match_dict, motif_idx, gene_name):

    # Subselect coords to gene
    coords_ = coords.loc[coords['gene_name']==gene_name]

    # Get gene associated sequences
    gene_assoc_sequences = get_sequences(sequences_record, coords_)

    # Compute background
    background = compute_background_freq(gene_assoc_sequences)

    # Compute motif matching
    motif = motifs_record[motif_idx]

    pattern_matches = {}
    fimo = FIMO(both_strands=True)
    pattern = fimo.score_motif(motif, gene_assoc_sequences, background)

    motif_match_dict[(gene_name, motif.accession.decode())] = pattern

def compute_motif_enrichment_(mdata, motif_match_df, prog_key='prog'):
    motif_enrich_df=None
    return motif_enrich_df

# Compute motif enrichment in enhancers or promoters associated with a gene
def compute_motif_enrichment(mdata, motif_file=None, seq_file=None, coords_file=None, 
                             n_jobs=1, prog_key='prog', data_key='rna', inplace=True, 
                             **kwargs):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata

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
    Parallel(n_jobs=n_jobs, 
             backend='threading')(delayed(perform_motif_match)(sequences_record, 
                                                               motifs_record, 
                                                               coords,
                                                               motif_match_dict, 
                                                               motif_idx, 
                                                               gene_name) \
                                                               for motif_idx in range(len(motifs_record)) \
                                                               for gene_name in tqdm(matching_gene_names,
                                                                                     desc='Motif scanning',
                                                                                     unit='genes'))

    # # Concatenate motif matches 
    # motif_match_df = pd.DataFrame(columns=['gene_name', 'motif_name', 'seq_name', 'start', 'end', 
    #                                        'strand', 'score', 'pvalue', 'qvalue'])

    # # Compute motif enrichment
    # motif_enrich_df = compute_motif_enrichment_(mdata, prog_key=prog_key, motif_match_df)

    # # Compute & update mudata with eval measure
    # # Store motif_match in varm same as gsea enrich

    # mdata[prog_key].var['eval_measure'] = None

    # if not inplace: return mdata[prog_key].var['eval_measure']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj')
    parser.add_argument('-n', '--n_jobs', default=1, typ=int)
    parser.add_argument('-pk', '--prog_key', default='prog', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str)
    parser.add_argument('-mf', '--motif_file', default=None, typ=str, required=True) 
    parser.add_argument('-sf', '--seq_file',  default=None, typ=str, required=True) 
    parser.add_argument('-cf', '--coords_file', default=None, typ=str, required=True) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()
    
    compute_motif_enrichment(args.mudataObj, motif_file=args.motif_file, 
                             seq_file=arg.seq_file, coords_file=arg.coords_file, 
                             n_jobs=args.n_jobs, prog_key=args.prog_key, 
                             data_key=args.data_key, inplace=args.output)





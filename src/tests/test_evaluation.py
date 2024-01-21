from ..evaluation import *

import logging
logging.basicConfig(level = logging.INFO)
from tqdm.auto import tqdm

from pathlib import Path

# Functional tests        
def test_explained_variance_ratio(mdata):
    try: compute_explained_variance_ratio(mdata)
    except: raise RuntimeError('Explained variance')
    
    try: assert round(mdata['prog'].var['explained_variance_ratio'].sum(),3)==1
    except: raise AssertionError('Explained variance')

# def test_categorical_association(mdata):
#    assert

def test_gene_enrichment(mdata):
    try: compute_geneset_enrichment(mdata)
    except: raise RuntimeError('Gene set enrichment')

def test_motif_enrichment(mdata, coords_file_loc, seq_file_loc):
    
    try: motif_file = read_motif_file((Path(__file__).parent /\
                                       '../../smk/resources/hocomoco_meme.meme').resolve())[:3]
    except: raise RuntimeError('Reading motif file')

    try: coords_file = read_coords_file(coords_file_loc)
    except: raise RuntimeError('Reading coordinate file')

    try: seq_file = read_sequence_file(seq_file_loc)
    except: raise RuntimeError('Reading sequence file')

    try: compute_motif_enrichment(mdata, motif_file, seq_file,
                                  coords_file, output_loc=None)
    except: raise RuntimeError('Motif enrichment')

if __name__=='__main__':

    mdata = create_test_data()

    test_explained_variance_ratio(mdata)
    # test_batch_association(mdata)
    test_gene_enrichment(mdata)
    # test_motif_enrichment(mdata, coords_file, seq_file_loc)
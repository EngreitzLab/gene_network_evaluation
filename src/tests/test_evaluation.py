# Allow this script to be run both as a standalone and as part of the package
if __name__ == "__main__" and __package__ is None:
    import os
    import sys
    current_loc = os.path.dirname(__file__)
    project_src = os.path.abspath(os.path.join(current_loc, ".."))  # .../src
    sys.path.insert(0, project_src)
    __package__ = "tests"

    # Use absolute imports; the block above ensures `evaluation` is on PYTHONPATH when run directly
    from evaluation import (
        compute_explained_variance_ratio,
        compute_categorical_association,
        compute_geneset_enrichment,
        compute_motif_enrichment,
    )

else:
    from ..evaluation import (
        compute_explained_variance_ratio,
        compute_categorical_association,
        compute_geneset_enrichment,
        compute_motif_enrichment,
    )

import mudata

from tqdm.auto import tqdm
from pathlib import Path

import logging
logging.basicConfig(level=logging.INFO)

# Functional tests        
def test_explained_variance_ratio(mdata, prog_key='prog'):
    try: compute_explained_variance_ratio(mdata, prog_key=prog_key, inplace=True)
    except: raise RuntimeError('Explained variance')
    
    # FIXME: See explained_variance computation in evaluations
    try: assert round(mdata[prog_key].var['explained_variance_ratio_X'].sum(),3)==1
    except: raise AssertionError('Explained variance')

    logging.info('Passed variance explained test!')

def test_categorical_association(mdata, prog_key='prog'):

    #TODO: Cover all cases

    try: compute_categorical_association(mdata, prog_key=prog_key, 
                                         categorical_key='celltype', 
                                         pseudobulk_key=None, test='dunn', 
                                         mode='one_vs_one', 
                                         n_jobs=1, inplace=True)
    except: raise RuntimeError('Categorical association - posthoc')

    try: compute_categorical_association(mdata, prog_key=prog_key, 
                                         categorical_key='celltype', 
                                         pseudobulk_key=None, test='pearsonr', 
                                         mode='one_vs_one', 
                                         n_jobs=1, inplace=True)
    except: raise RuntimeError('Categorical association - correlation')

    logging.info('Passed categorical association test!')
    
def test_gene_enrichment(mdata, prog_key='prog'):
    try: compute_geneset_enrichment(mdata, prog_key=prog_key, inplace=True)
    except: raise RuntimeError('Gene set enrichment')

    logging.info('Passed geneset enrichment test!')

def test_motif_enrichment(mdata, seq_file_loc, prog_key='prog'):

    motif_file = (Path(__file__).parent /\
                 './test_data/motifs.meme').resolve()

    loci_file = (Path(__file__).parent /\
                 './test_data/p2g_links.txt').resolve()
    
    try: compute_motif_enrichment(mdata, 
                                  prog_key=prog_key,
                                  motif_file=motif_file, 
                                  seq_file=seq_file_loc,
                                  loci_file=loci_file)
    except: raise RuntimeError('Motif enrichment')

    logging.info('Passed motif enrichment test!')

if __name__=='__main__':

    mdata = mudata.read((Path(__file__).parent /\
                         './test_data/sample_mudata.h5mu').resolve())

    # test_explained_variance_ratio(mdata, prog_key='cNMF')
    test_categorical_association(mdata, prog_key='cNMF')
    test_gene_enrichment(mdata, prog_key='cNMF')
    test_motif_enrichment(mdata, seq_file_loc=sys.argv[1], prog_key='cNMF')

    logging.info('Passed all tests!')
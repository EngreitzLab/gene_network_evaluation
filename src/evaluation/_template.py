import os
import argparse

import mudata

# The template is intended to be a general guide on how to organise your code
# so feel free to write functions as required.
# Extended the parser and function with specific parameters required for
# your benchmark.
# For the parameters already named in the parser, they relate to the 
# standard keys in our mudata input structure and should not be changed.

# Try and organise your code into nested functions
def helper_function(mdata, intermediate=None):
    return

# Rename function to compute_{eval_measure}_
def compute_eval_measure_(some_values, somearg=None):
    return eval_measure

# Rename function to compute_{eval_measure}
def compute_eval_measure(mdata, prog_key='prog', data_key='rna', 
                         somearg=None,  n_jobs=1, inplace=True, 
                         **kwargs):
    """
    Computes eval_measure.

    ARGS
        mdata : MuData
            mudata object containing anndata of program scores and cell-level metadata.
        prog_key: 
            index for the gene program anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        somearg: str (optional)
            description of somearg.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy
    
    RETURNS 
        if not inplace:
            mdata[prog_key].var['eval_measure']
            
    """
    
    # Only copy the anndatas that will be used in the eval
    if not inplace:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})
  
    intermediate_output = helper_function(mdata, 
                                          intermediate=intermediate_input)

    # Compute & update mudata with eval measure
    mdata[prog_key].var['eval_measure'] = compute_eval_measure_(intermediate_output, 
                                                                somearg=somearg)

    # Not all eval measures might be condensed at a program level
    # Store eval measure or additional outputs in
    # var - program indexed series
    # obs - cell indexed series
    # varm - program indexed array
    # obsm - cell indexed array
    # uns - unstructured information
    # obsp - cell x cell structured outputs
    # varp - program x program structured outputs

    # See https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

    # Return only the evaluations if not in inplace mode
    if not inplace: return mdata[prog_key].var['eval_measure']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--somearg', default=None)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-dk', '--data_key', default='rna', type=str) 
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    compute_eval_measure(mdata, prog_key=args.prog_key, data_key=args.data_key, 
                         somearg=args.somearg, n_jobs=args.n_jobs, inplace=args.output)


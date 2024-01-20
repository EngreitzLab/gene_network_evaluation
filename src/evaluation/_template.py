import os
import argparse

# The template is intended to be a general guide on how to organise your code
# so feel free to write functions as required.
# Extended the parser and function with specific parameters required for
# your benchmark.
# For the parameters already named in the parser, they relate to the 
# standard keys in our mudata input structure and should not be changed.

# Try and organise your code into nested functions
def helper_function(mdata, intermediate=None):
    return

def compute_eval_measure_(some_values):
    return eval_measure

# Rename function to compute_{eval_measure}
def compute_eval_measure(mdata, n_jobs=-1, prog_key='prog', 
                         rna_key='rna', atac_key='atac', inplace=True):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata
  
    intermediate_output = helper_function(mdata, 
                                          intermediate=intermediate_input)

    # Compute & update mudata with eval measure
    mdata[prog_key].var['eval_measure'] = compute_eval_measure_(intermediate_output)

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

    if not inplace: return mdata[prog_key].var['eval_measure']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-n', '--n_jobs', default=1, type=int)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-rk', '--rna_key', default='rna', type=str) 
    parser.add_argument('-ak', '--atac_key', default='atac', type=str) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    compute_eval_measure(mdata, n_jobs=args.n_jobs,prog_key=args.prog_key, 
                         rna_key=args.rna_key, atac_key=args.atac_key, inplace=args.output)


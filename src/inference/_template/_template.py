import os
import gin
import argparse

import mudata
import anndata

# NOTE: _template_env.yml contains the essential pacakges that must be present
# in your conda env. Add more packages as neccessary.

# Write a function to run your method

# Most parameters of this function will be supplied by the user 
# using a gin-config file. Generally those should include
# all the user tunable parameters specific for this method.
# Read more here https://github.com/google/gin-config/blob/master/docs/index.md

# Other general parameters such as number of workers etc.
# are supplied separately as they will be common in the 
# pipeline that will use these functions.
# Replace name with run_{program_inference_method}_
@gin.configurable
def run_program_inference_method_(anndata, layer='X',
                                  method_specific_param_1=None, 
                                  method_specific_param_2=None,
                                  method_specific_param_3=None,
                                  method_specific_param_4=None,
                                  n_jobs=-1):
  
    # Write code to run your method
    # You can write additional functions in this script
    # Or create more complex scripts in the scripts folder
    # and import the functions 

    # Place any imports that will only be present
    # in the method env (_template_env.yml) inside
    # this function
    
    return cell_program_scores, program_feature_scores

# Use this function to store method outputs in the relevant
# mudata keys. This function will load the gin config file
# and call the previous function. Add parameters as neccessary
# but use the named parameters in the parser as they relate
# to the keys specified in our mudata input/output specification
# Replace name with run_{program_inference_method}
def run_program_inference_method(mdata, work_dir=None, scratch_dir=None,
                                 prog_key='program_inference_method', 
                                 data_key='rna', layer='X', config_path=None, 
                                 n_jobs=-1, inplace=True):
                                 
    """
    Perform gene program inference using {inference_method}.

    ARGS:
        mdata : MuData
            mudata object containing anndata of data and cell-level metadata.
        work_dir: str
            path to directory to store outputs.
        scratch_dir: str
            path to directory to store intermediate outputs.
        prog_key: 
            index for the anndata object (mdata[prog_key]) in the mudata object.
        data_key: str
            index of the genomic data anndata object (mdata[data_key]) in the mudata object.
        layer: str (default: X)
            anndata layer (mdata[data_key].layers[layer]) where the data is stored.
        config_path: str
            path to gin configurable config file containing method specific parameters.
        n_jobs: int (default: 1)
            number of threads to run processes on.
        inplace: Bool (default: True)
            update the mudata object inplace or return a copy

    RETURNS:
        if not inplace:
            mdata
    
    """

    # Load method specific parameters
    try: gin.parse_config_file(config_path)
    except: raise ValueError('gin config file could not be found')

    if not inplace:
        mdata = mudata.MuData({data_key: mdata[data_key].copy()})
  
    # Compute (structure this as neccessary)
    # Not that method specific params loaded via gin do not
    # have to be passed explicity
    cell_program_scores, program_feature_scores = \
    run_program_inference_method_(mdata[data_key], 
                                  layer = layer, 
                                  n_jobs=n_jobs)

    # Store cell x program and program x feature matrices
    mdata.mod[prog_key] = anndata.AnnData(data=cell_program_scores, 
                                          obs=mdata[data_key].obs)
    mdata[prog_key].varm['loadings'] = program_feature_scores

    # Store additional outputs in
    # var - program indexed series
    # obs - cell indexed series
    # varm - program indexed array
    # obsm - cell indexed array
    # uns - unstructured information
    # obsp - cell x cell structured outputs
    # varp - program x program structured outputs

    # See https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

    if not inplace: return mdata

if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--work_dir', default='./', type=str)
    parser.add_argument('--scratch_dir', default='./', type=str)
    parser.add_argument('-pk', '--prog_key', default='program_inference_method', typ=str) 
    parser.add_argument('-dk', '--data_key', default='rna', typ=str) # could be atac
    parser.add_argument('--layer', default='X', type=str) # layer of data anndata to use for inference
    parser.add_argument('--config_path', default='./_template_config.gin', type=str)
    parser.add_argument('-n', '--n_jobs', default=-1, type=int)
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    mdata = mudata.read(args.mudataObj_path)
    run_program_inference_method(mdata, work_dir=args.work_dir, scratch_dir=args.scratch_dir,
                                 prog_key=args.prog_key, data_key=args.data_key, layer=args.layer, 
                                 config_path=args.config_path, n_jobs=args.n_jobs, inplace=args.output)
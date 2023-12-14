import os
import sys

import mudata as mu

import logging
logging.basicConfig(filename=snakemake.log[0],
                    level = logging.INFO)

# Function to verify data keys
def check_mdata(mdata, config):

    # Check if program scores exist
    try: assert mdata[config['prog_key']]
    except: raise ValueError('Program scores not found.')

    try: assert mdata[config['data_key']]
    except: raise ValueError('Data not found.')

    # Check if expression data exists
    if config['rna_key'] is not None:
        try: assert mdata[config['rna_key']]
        except: raise ValueError('Expression data not found.')

        #TODO: Check if same as used for prog_score

    # Check if accessibility data exists
    if config['atac_key'] is not None:
        try: assert mdata[config['atac_key']]
        except: raise ValueError('Accessiblity data not found.')

    # Check if index of data and prog match
    assert (mdata[config['prog_key']].obs.index==mdata[config['data_key']].obs.index).all()

    # Check if batch key exists in data and prog
    if config['batch_key'] is not None:
        bk = config['batch_key']

        for key in [('prog_key', 'data_key'), 
                    ('data_key', 'prog_key')]:
            key_a, key_b = key[0], key[1]
            if bk in mdata[config[key[0]]].obs.columns:
                if bk in mdata[config[key[1]]].obs.columns:
                    assert (mdata[config[key[0]]].obs[bk]==mdata[config[key[1]]].obs[bk]).all()
                else:
                    mdata[config[key[1]]].obs[bk]=mdata[config[key[0]]].obs[bk]
            else:
                raise ValueError('Batch annotations not found')

    # Check if gene names match organism
    if config['organism'] is not None:
        var_names = mdata[config['data_key']].var_names
        if config['organism']=='human':
            if sum([nam.isupper() for nam in var_names[:10]])<8: # Why not exact?
                logging.warn('Gene names are not in human (ABCD) format')
        elif config['organism']=='mouse':
            if sum([nam.istitle() for nam in var_names[:10]])<8:
                logging.warn('Gene names are not in mouse (Abcd) format')


# Function to load mudata input
def load_mdata(config):

    input_loc = config['input_loc']

    # Check location exists
    try: os.path.exists(input_loc)
    except: raise ValueError('Input does not exist')

    # Check if input is dir or a file
    if os.path.isdir(input_loc):
        fil_nam = [fil for fil in os.listdir(input_loc) if '.h5mu' in fil]

        # Check if unique input
        if len(fil_nam)==0:
            raise ValueError('Input does not exist')
        elif len(fil_nam)>1:
            raise ValueError('Provide unique input') 
        elif len(fil_nam)==1:
            input_loc = os.path.join(input_loc, fil_nam[0])

    if input_loc.endswith('.h5mu'):
        mdata = mu.read(input_loc)
    else:
        raise ValueError('Provide input as a .h5mu file.')

    check_mdata(mdata, config)

    return mdata

# Execution (assumes Snakemake)
with open(snakemake.log[0], "a") as f:
    f.write("\nstderr:\n")
    sys.stderr = sys.stdout = f

    mdata = load_mdata(snakemake.config)
    logging.info('Successfully loaded mudata.')

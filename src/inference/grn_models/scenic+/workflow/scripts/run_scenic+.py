import os
import dill
import argparse
import muon as mu
import numpy as np
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--path_cistopic_obj', required=True)
parser.add_argument('-m','--path_motif_enrichment', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-s', '--path_scenicplus_obj', required=True)
parser.add_argument('-g', '--path_grn', required=True)
parser.add_argument('-r', '--path_r2g', required=True)
parser.add_argument('-t', '--path_tri', required=True)
parser.add_argument('-j','--n_cpu', required=True)
parser.add_argument('-x', '--temp_dir', required=True)
parser.add_argument('-b', '--biomart_host', required=True)
parser.add_argument('-u', '--min_upstream', required=True)
parser.add_argument('-v', '--max_upstream', required=True)
parser.add_argument('-w', '--min_downstream', required=True)
parser.add_argument('-z', '--max_downstream', required=True)

args = vars(parser.parse_args())
path_input = args['path_input']
path_cistopic_obj = args['path_cistopic_obj']
path_motif_enrichment = args['path_motif_enrichment']
organism = args['organism']
path_scenicplus_obj = args['path_scenicplus_obj']
path_grn = args['path_grn']
path_r2g = args['path_r2g']
path_tri = args['path_tri']
n_cpu = int(args['n_cpu'])
temp_dir = args['temp_dir']
biomart_host = args['biomart_host']
min_upstream = int(args['min_upstream'])
max_upstream = int(args['max_upstream'])
min_downstream = int(args['min_downstream'])
max_downstream = int(args['max_downstream'])

# Read rna adata
adata = mu.read(path_input)
del adata.mod['atac']
obs = adata.obs.copy()
adata = adata.mod['rna'].copy()
adata.obs = obs
adata.X = adata.layers["counts"]

# Read cisTopic object
cistopic_obj = dill.load(open(path_cistopic_obj, 'rb'))

# Load motif enrichment object
menr = dill.load(open(path_motif_enrichment, 'rb'))

# Create SCENICPLUS object
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

# Select dbs
if organism == 'human':
    species = "hsapiens"
    assembly = "hg38"
    tf_file = "resources/tf_lists/human.txt"

# Run SCENIC+
path_output = os.path.dirname(path_scenicplus_obj)
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_celltype'],
        species = species,
        assembly = assembly,
        tf_file = tf_file,
        save_path = path_output,
        biomart_host = biomart_host,
        upstream = [min_upstream, max_upstream],
        downstream = [min_downstream, max_downstream],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = False,
        export_to_loom_file = False,
        export_to_UCSC_file = False,
        path_bedToBigBed = 'resources/bin',
        n_cpu = n_cpu,
        _temp_dir=temp_dir
    )
except Exception as e:
    dill.dump(scplus_obj, open(path_scenicplus_obj, 'wb'), protocol=-1)
    print(e, "\nSCENIC+ object saved to", path_scenicplus_obj)

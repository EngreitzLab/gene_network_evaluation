import os
import pandas as pd
import numpy as np
from celloracle import motif_analysis as ma
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser(
    prog="python r2g.py",
    description="Build region to gene links using Cicero coaccessibility scores."
)
parser.add_argument('-d','--path_data', required=True, type=str, help="path to input MuData object, this is not modified throughout pipeline")
parser.add_argument('-a','--all_peaks', required=True, type=str, help="path to scATAC-seq peak list output from r2g.R")
parser.add_argument('-c','--connections', required=True, type=str, help="path to Cicero coaccessibility scores output from r2g.R")
parser.add_argument('-g','--genome', required=True, type=str, help="version of the genome to use, e.g. hg38, mm10")
parser.add_argument('-t','--coaccess_thr', required=True, type=float, help="filter out Cicero connections below this threshold")
parser.add_argument('-o','--path_out', required=True, type=str, help="path to output r2g.csv containing region to gene links")
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_all_peaks = args['all_peaks']
path_connections = args['connections']
genome = args['genome']
coacess_thr = float(args['coaccess_thr'])
path_out = args['path_out']

# Load scATAC-seq peak list
peaks = pd.read_csv(path_all_peaks, index_col=0).x.values.astype('U')
peaks = np.char.replace(peaks, '-', '_')

# Load Cicero coaccessibility scores
cicero_connections = pd.read_csv(path_connections, index_col=0)
cicero_connections['Peak1'] = np.char.replace(cicero_connections['Peak1'].values.astype('U'), '-', '_')
cicero_connections['Peak2'] = np.char.replace(cicero_connections['Peak2'].values.astype('U'), '-', '_')

# Extract tss information
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=genome)

# Integrate
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)

# Process
integrated = integrated[integrated['coaccess'] >= coacess_thr]
integrated['peak_id'] = integrated['peak_id'].str.replace('_', '-')
integrated = integrated.rename(columns={'peak_id': 'cre', 'gene_short_name': 'gene', 'coaccess': 'score'})
integrated = integrated.sort_values(['cre', 'score'], ascending=[True, False])
integrated["pval"] = np.nan

# Remove unexpressed genes
genes = mu.read(os.path.join(path_data, 'rna')).var.index.values.astype('U')
integrated = integrated[integrated['gene'].isin(genes)]

# Write
integrated.to_csv(path_out, index=False)

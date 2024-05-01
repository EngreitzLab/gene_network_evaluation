import pandas as pd
import numpy as np
import celloracle as co
from celloracle import motif_analysis as ma
import muon as mu
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--organism', required=True)

parser.add_argument('-k','--path_peaks', required=True)
parser.add_argument('-c','--path_cicero', required=True)
parser.add_argument('-q','--thr_coaccess', required=True)
parser.add_argument('-j','--path_r2g', required=True)

parser.add_argument('-f','--path_motifs', required=True)
parser.add_argument('-e','--path_TF2peak', required=True)

parser.add_argument('-l','--path_links', required=True)
parser.add_argument('-p','--thr_edge_pval', required=True)
parser.add_argument('-t','--thr_top_edges', required=True)
parser.add_argument('-g','--path_grn', required=True)

parser.add_argument('-b','--path_basegrn', required=True)
parser.add_argument('-r','--path_tri', required=True)

parser.add_argument('-m', '--path_mdata', required=True)
args = vars(parser.parse_args())

# Global
path_input = args['path_input']
organism = args['organism']

# Peak-gene links from Cicero
path_peaks = args['path_peaks']
path_cicero = args['path_cicero']
thr_coaccess = float(args['thr_coaccess'])
path_r2g = args['path_r2g']

# TF-peak links from gimmemotifs scan
path_motifs = args['path_motifs']
path_TF2peak = args['path_TF2peak']

# TF2g links
path_links = args['path_links']
thr_edge_pval = float(args['thr_edge_pval'])
thr_top_edges = int(args['thr_top_edges'])
path_grn = args['path_grn']

# Triplet links
path_basegrn = args['path_basegrn']
path_tri = args['path_tri']

# Mdata
path_mdata = args['path_mdata']

# Determine genome
if organism == 'human':
    genome = 'hg38'
elif organism == 'mouse':
    genome = 'mm10'

# Read mdata to add objects to
data = mu.read(path_input)

# Peak-gene links from Cicero
peaks = pd.read_csv(path_peaks, index_col=0).x.values.astype('U')
peaks = np.char.replace(peaks, '-', '_')
cicero_connections = pd.read_csv(path_cicero, index_col=0)
cicero_connections['Peak1'] = np.char.replace(cicero_connections['Peak1'].values.astype('U'), '-', '_')
cicero_connections['Peak2'] = np.char.replace(cicero_connections['Peak2'].values.astype('U'), '-', '_')
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=genome)
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated, cicero_connections=cicero_connections)
integrated = integrated[integrated.coaccess >= thr_coaccess]
integrated = integrated.rename(columns={"peak_id": "source", "gene_short_name": "target", "coaccess": "weight"})
integrated["pvals"] = 1
integrated = integrated[["source", "target", "weight", "pvals"]]
integrated.to_csv(path_r2g, index=False)
data.uns["r2g"] = integrated

# TF-peak links from gimmemotifs scan
motifs = ma.load_TFinfo(path_motifs)
motifs_df = motifs.scanned_df.sort_values("score", ascending=False)
motifs_df["source"] = motifs_df.apply(lambda x: x["factors_indirect"] if x["factors_direct"] == "" else x["factors_direct"], axis=1)
motifs_df = motifs_df.assign(source=motifs_df.source.str.split(",")).explode("source")
motifs_df = motifs_df.rename(columns={"seqname": "target", "score": "weight"})
motifs_df["pvals"] = 1
motifs_df = motifs_df[["source", "target", "weight", "pvals", "motif_id", "factors_direct", "factors_indirect", "pos", "strand"]]
motifs_df = motifs_df[motifs_df.source != ""]
motifs_df.to_csv(path_TF2peak, index=False)
data.uns["TF2r"] = motifs_df

# Filter links and concatenate
links = co.load_hdf5(file_path=path_links)
links.filter_links(p=thr_edge_pval, weight="coef_abs", threshold_number=thr_top_edges)
grn = []
for celltype in links.links_dict.keys():
    tmp = links.filtered_links[celltype].dropna()
    tmp['celltype'] = celltype
    grn.append(tmp)
grn = pd.concat(grn)[['source', 'target', 'coef_mean', 'p', 'celltype']]
grn = grn.rename(columns={'coef_mean': 'weight', 'p': 'pvals'}).sort_values(['source', 'target']).reset_index(drop=True)
grn.to_csv(path_grn, index=False)
data.uns["grn"] = grn

# Extract triplets and filter by grn
tri = pd.read_csv(path_basegrn, index_col=0)
tri = tri.melt(id_vars=['peak_id', 'gene_short_name'], var_name='TF', value_name='interaction')
tri = tri[tri['interaction'] != 0]
tri = tri[['TF', 'gene_short_name', 'peak_id']].rename(columns={'TF': 'source', 'peak_id': 'region', 'gene_short_name': 'target'})
tri = tri.sort_values(['source', 'target', 'region'])
tri = tri.merge(grn.groupby(['source', 'target']).size().reset_index(), on=['source', 'target'])
tri = tri[['source', 'target', 'region']]
tri['region'] = [r.replace('_', '-') for r in tri['region']]
tri.to_csv(path_tri, index=False)
data.uns["tri"] = tri

# Save data
data.write_h5mu(path_mdata)

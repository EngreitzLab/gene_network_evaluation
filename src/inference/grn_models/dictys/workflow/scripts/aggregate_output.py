import pandas as pd
import numpy as np
import argparse, os, sys
import pickle
import mudata as md
import dictys

"""
5th script for dictys pipeline. 
Assume that there exist a data directory, storing processed mudata, fragment file, mapping file, etc.

This script aggregates inferred GRN for each cell cluster.
"""

parser = argparse.ArgumentParser(description="Aggregate results across individual cell cluster's working directory", usage="")
parser.add_argument('--mudata', type=str, help="Processed mudata, containing clustering results")
parser.add_argument('--input_dir', nargs='*', help="Individual cell cluster's working directory")
parser.add_argument('--output_dir', help="Output directory")

args = parser.parse_args()
mudata = args.mudata
input_dir = args.input_dir
output_dir = args.output_dir
path_TF2peak = os.path.join(output_dir, "TF2r.csv")
path_r2g = os.path.join(output_dir, "r2g.csv")
path_grn = os.path.join(output_dir, "grn.csv")
path_mdata = os.path.join(output_dir, "dictys.h5mu")
                          
print(input_dir)
print(output_dir)

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Read mdata to add objects to
data = md.read(mudata)


# Peak-gene links from distance to closest TSS
peakgene_df = list()
for d in input_dir:
    df = pd.read_csv(os.path.join(d, "tssdist.tsv.gz"), sep="\t")
    df['celltype'] = os.path.basename(d)[6:]
    df['dist'] = - np.abs(df['dist']) / 1e6 * np.log(10)
    df = df.rename(columns={"dist": "weight", 'region': 'source'})
    peakgene_df.append(df)
peakgene_df = pd.concat(peakgene_df)
peakgene_df['pval'] = 1
peakgene_df.to_csv(path_r2g, index=False)
data.uns["r2g"] = peakgene_df.copy(True)

    
# TF-peak links from wellington and homer
motifs_df = list()
for d in input_dir:
    df = pd.read_csv(os.path.join(d, "binding.tsv.gz"), sep="\t")
    df['celltype'] = os.path.basename(d)[6:]
    df = df.rename(columns={"score": "weight", 'TF': 'source', 'loc': 'target'})
    motifs_df.append(df)
motifs_df = pd.concat(motifs_df)
motifs_df['pval'] = 1
motifs_df.to_csv(path_TF2peak, index=False)
data.uns["TF2r"] = motifs_df.copy(True)


# TF-gene links 
tftarget_df = list()
for d in input_dir:
    weights = pd.read_csv(os.path.join(d, "net_nweight.tsv.gz"), sep="\t", index_col=0)
    mask = pd.read_csv(os.path.join(d, "binlinking.tsv.gz"), sep="\t", index_col=0)
    mask = mask[weights.columns].loc[weights.index]

    df = list()
    for i in np.arange(weights.shape[0]):
        tmpdf = [(weights.index[i], weights.columns[j], weights.iloc[i,j]) 
                 for j in np.arange(weights.shape[1]) if mask.iloc[i,j]]
        df.extend(tmpdf)
    df = np.array(df)
    df = pd.DataFrame(df, columns=['source', 'target', 'weight'])
    df['celltype'] = os.path.basename(d)[6:]
    tftarget_df.append(df)
tftarget_df = pd.concat(tftarget_df)
tftarget_df['pval'] = 1
tftarget_df.to_csv(path_grn, index=False)
data.uns["grn"] = tftarget_df.copy(True)

data.write_h5mu(path_mdata)


























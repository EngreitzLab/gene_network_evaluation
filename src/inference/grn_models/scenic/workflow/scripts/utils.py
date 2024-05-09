# -*- coding: utf-8 -*-

"""
Python script with helper functions for runnning and analyzing SCENICprotocol runs
"""


import pandas as pd
import loompy as lp


# Read a single adjecency file output by pyscenic grn
def read_adj(file, tfs=None):
    adj = pd.read_csv(file, sep="\t")
    if tfs != None:
        if isinstance(tfs, str):
            tfs = [tfs]
        adj = adj[adj["TF"].isin(tfs)]
    return adj


# Read multiple adjecency files output by pyscenic grn and rank each target per TF and file.
# Outputs a dataframe with TFs and their targets ranked by feature importance within TF and file groups
def read_multiple_adj(files, tfs=None, rank=True):
    adj = pd.DataFrame()
    for file in sorted(files):
        curr_adj = read_adj(file, tfs=tfs)
        curr_adj["file"] = file
        adj = pd.concat([adj, curr_adj])
    adj["rank"] = adj.groupby(["TF", "file"])["importance"].rank("dense", ascending=False).astype(int)
    return adj


# Takes in output from read_multiple_adj and outputs the average rank and importance of targets within each TF across filetypes
# Outputs a list of TF-target pairs sorted by the mean rank of their link across files
def rank_importance(adj):
    ranks = adj.groupby(["TF", "target"]).agg({"rank": "mean", "importance": ["mean", "std"]})
    ranks = ranks.fillna(0)
    ranks = ranks.reset_index().sort_values(["TF"]).groupby(["TF"], sort=False).apply(lambda x: x.sort_values(("rank", "mean"), ascending=True)).drop("TF", axis=1)
    ranks.index = ranks.index.droplevel(1)
    ranks = ranks.reset_index()
    return ranks


# Combines read_multiple_adj and rank_importance to give both views of the adjacencies across files
def parse_adj(files, tfs=None):
    adj = read_multiple_adj(files, tfs=tfs)
    gene_rank = rank_importance(adj)
    return adj, gene_rank


# Takes in a single or set of loom files with a tf and returns the number of times each gene in the
# dataset was in that tfs regulon along with the per cell activity of the regulon in each file
def parse_loom(files, tf, cells, genes, norm="num_cells"):
    targets = pd.Series(index=genes, data=0, dtype=int)
    activity = pd.DataFrame(index=cells, columns=sorted(files))
    regulon = tf + "(+)"
    for file in sorted(files):
        print(file)
        lf = lp.connect(file, mode='r+', validate=False)

        # Collect regulon targets and activity
        if regulon in lf.ra.Regulons.dtype.names:
            curr_targets = pd.Series(index=lf.ra.Gene, data=lf.ra.Regulons[regulon])
            targets = targets + curr_targets.loc[targets.index]
            curr_activity = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)[regulon]
            activity[file] = curr_activity.loc[activity.index]
        else:
            print("{} regulon not in file".format(regulon))
        # Collect regulon activity
        #curr_activity = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)[regulon]
        #activity[file] = curr_activity.loc[activity.index]

        lf.close()

    targets = targets[targets > 0].sort_values(ascending=False)
    #activity = (activity/activity.sum(axis=0))*len(cells)
    #activity = (activity/activity.sum(axis=0))*activity.sum(axis=0).median()  scale to median
    return targets, activity

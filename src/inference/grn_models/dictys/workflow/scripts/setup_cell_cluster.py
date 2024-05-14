import pandas as pd
import numpy as np
import argparse, os, sys
import pickle
import mudata as md
import dictys

"""
1st script for dictys pipeline. 
Assume that there exist a data directory, storing processed mudata, fragment file, mapping file, etc.
This script will 
    1. create individual folder for each cluster of cells, 
    2. create a file storing all cell barcodes belonging to that cluster,
    3. create a bed file storing the consensus peak list
    4. create a raw expression matrix for cells in that cluster
"""


parser = argparse.ArgumentParser(description="Set up individual cell cluster's working directory", usage="")
parser.add_argument('--mudata', type=str, help="Processed mudata, containing clustering results")
parser.add_argument('--data_dir', type=str, default="./data/static", help="Individual cell cluster's working directory")
parser.add_argument('--bc_mapping', type=str, default=None, help="Path to the file containing mapping between mudata cell barcode to fragment file barcode")


args = parser.parse_args()
mudata = args.mudata
data_dir = args.data_dir
bc_mapping = args.bc_mapping


def barcode_mapping(bc_mapping):
    with open(os.path.join(bc_mapping), "rb") as f:
        mapping = pickle.load(f)
    return mapping


def setup_directory(mudata, data_dir, bc_mapping):
    """
    Set up sub-directories for all cell clusters
    """
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    if bc_mapping is not None:
        mapping = barcode_mapping(bc_mapping)
        
    data = md.read(os.path.join(mudata))
    clus = data.obs["celltype"].unique()
    for count, c in enumerate(clus):
        subset_dir = os.path.join(data_dir, f"Subset{count}")
        if not os.path.exists(subset_dir):
            os.mkdir(subset_dir)
        
        curr_clus = data[data.obs["celltype"]==c]
        curr_rna  = curr_clus.mod["rna"]
        curr_rna_X = pd.DataFrame(np.array(curr_rna.layers['counts'].todense()), columns=curr_rna.var.index, index=curr_rna.obs.index).T
        
        raw_rna_filename = os.path.join(subset_dir, "expression0.tsv.gz")
        proc_rna_filename = os.path.join(subset_dir, "expression.tsv.gz")
        curr_rna_X.to_csv(raw_rna_filename, sep="\t", compression="gzip") 
        dictys.preproc.qc_reads(raw_rna_filename, proc_rna_filename, 50, 10, 0, 0, 0, 0)  # all cell filters are turned off
        
        ref_name = pd.read_csv(proc_rna_filename, sep="\t", index_col=0, nrows=0).columns.tolist()
        if bc_mapping is not None:
            ref_name = [mapping[n] for n in ref_name if n in mapping]# and '_4' in mapping[n]]
        
        atac_fname = os.path.join(subset_dir, "names_atac0.txt")
        with open(atac_fname, "w") as f:
            for n in ref_name:
                f.write(f"{n}\n")
                
        curr_atac = curr_clus.mod["atac"]
        curr_atac_chr = np.unique([n.split(":")[0] for n in curr_atac.var.index]) # assume genomic region follow the pattern chr:start-end
        curr_atac_chr = dict([(n,0) for n in curr_atac_chr if 'chr' in n])

        print(f"Subset {count}: cell type {c}")
        celltypeid = os.path.join(subset_dir, "celltype_label.txt")
        with open(celltypeid, "w") as f:
            f.write(f"Subset {count}: cell type {c}\n")
        
        all_atac = os.path.join(subset_dir, "all_peak.bed")
        with open(all_atac, "w") as f:
            regions = data.mod["atac"].var.index
            for r in regions:
                chrs, loci = r.split(":")
                start, end = loci.split("-")
                if chrs not in curr_atac_chr:
                    continue
                if int(end)-int(start) < 100:
                    continue
                f.write(f"{chrs}\t{start}\t{end}\n")
    
setup_directory(mudata=mudata, data_dir=data_dir, bc_mapping=bc_mapping)


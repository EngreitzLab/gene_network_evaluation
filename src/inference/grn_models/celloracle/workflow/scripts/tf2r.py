import pandas as pd
import numpy as np
import os
from celloracle import motif_analysis as ma
from genomepy import Genome, config
from gimmemotifs.motif import default_motifs
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser(
    prog="python tf2r.py",
    description="Build region to TF links using GimmmeMotifs motif scanning."
)
parser.add_argument('-d','--path_data', required=True, type=str, help="path to input MuData object, this is not modified throughout pipeline")
parser.add_argument('-r','--path_r2g', required=True, type=str, help="path to region to gene links output from r2g.py")
parser.add_argument('-gd','--genome_dir', required=True, type=str, help="path to the genomepy reference genome to use, must be downlaoded first using download_genomes.py")
parser.add_argument('-g','--genome', required=True, type=str, help="version of the genome to use, e.g. hg38, mm10")
parser.add_argument('-f','--fpr', required=True, type=float, help="false positive rate for GimmmeMotifs motif scanning")
parser.add_argument('-b','--blen', required=True, type=int, help="background sequence length for GimmmeMotifs motif scanning")
parser.add_argument('-t','--tfb_thr', required=True, type=float, help="threshold for filtering TF binding predictions")
parser.add_argument('-p','--threads', required=False, type=int, help="number of threads to use for motif scanning")
parser.add_argument('-ti','--path_tfinfo', required=True, type=str, help="path to output TFinfo object")
parser.add_argument('-o','--path_out', required=True, type=str, help="path to output tf2r.csv containing TF to region links")
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_r2g = args['path_r2g']
genome_dir = args['genome_dir']
genome = args['genome']
fpr = float(args['fpr'])
blen = int(args['blen'])
tfb_thr = float(args['tfb_thr'])
threads = args['threads'] if args['threads'] is not None else os.cpu_count()
path_tfinfo = args['path_tfinfo']
path_out = args['path_out']

# Load annotated peak data.
peaks = pd.read_csv(path_r2g)
if peaks.shape[0] == 0:
    tfb = pd.DataFrame(columns=['cre', 'tf', 'score'])
    tfb.to_csv(path_out, index=False)
    exit()
peaks['cre'] = peaks['cre'].str.replace('-', '_')

def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

def check_peak_format(peaks_df, gname, genome_dir):
    """
    Check peak format.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)

    """

    df = peaks_df.copy()
    df = df.rename(columns={'cre': 'peak_id', 'gene': 'gene_short_name'})
    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(int)
    df_decomposed["end"] = df_decomposed["end"].astype(int)

    # Load genome data
    genome_data = Genome(gname, genomes_dir=genome_dir)
    all_chr_list = list(genome_data.keys())

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]

    # Print
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df

# Format and delete peaks
peaks = check_peak_format(peaks, gname=genome, genome_dir=genome_dir)

# Instantiate TFinfo object
tfi = ma.TFinfo(
    peak_data_frame=peaks,
    ref_genome=genome,
    genomes_dir=genome_dir,
)

# Filter motifs based on gene expression
motifs = default_motifs()
dic_motif2TFs = ma.tfinfo_core._get_dic_motif2TFs(
    species='Human',
    motifs=motifs,
    TF_evidence_level='direct_and_indirect',
    formatting='auto'
)
genes = mu.read(os.path.join(path_data, 'rna')).var.index.values.astype('U')
m_keys = []
for m in dic_motif2TFs:
    tfs = np.unique(dic_motif2TFs[m]).astype('U')
    if np.any(np.isin(tfs, genes)):
        m_keys.append(m)
motifs = [m for m in motifs if m.id in m_keys]

# Update config
config.config.config['genomes_dir'] = genome_dir

# Scan
tfi.scan(
    background_length=blen,
    fpr=fpr,
    motifs=motifs,  # Use filtered motifs
    verbose=True,
    n_cpus=threads
)

# Do filtering
tfi.filter_motifs_by_score(threshold=tfb_thr)

# Save TFinfo object
tfi.to_hdf5(file_path=path_tfinfo)

# Extract filtered TF predictions
df = tfi.scanned_filtered[["seqname", "motif_id", "score"]].copy()
df['motif_values'] = df['motif_id'].map(tfi.dic_motif2TFs)
df = df.explode('motif_values')
df = df[df['motif_values'].isin(genes)]  # Filter TFs not in GEX
df = df.groupby(['seqname', 'motif_values'])['score'].max().reset_index()
df = df[['seqname', 'motif_values', 'score']].dropna()
df = df.reset_index(drop=True).rename(columns={'seqname': 'cre', 'motif_values': 'tf'})
df['cre'] = df['cre'].str.replace('_', '-')
df = df.sort_values(['cre', 'score'], ascending=[True, False])
df["pval"] = np.nan

# Write
df.to_csv(path_out, index=False)

print('Done')
os._exit(0)  # Add this else it gets stuck

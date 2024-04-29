import pandas as pd
import argparse, os, sys

"""
2nd script for dictys pipeline. 
Assume that there exist a data directory and individual sub-folder for each cluster.
This script and many of the following scripts need to be ran for each cluster individually.

This script will create a sorted bam file for the cluster of cells, call peaks from the bam file, and intersect the called peaks with the consensus peak list.
"""

parser = argparse.ArgumentParser(description="Splits fragment file by annotated cell clusters and builds .bam file", usage="")
parser.add_argument('--cluster_number', type=int, help="Extract fragments for a cluster and build a binary alignment file")
parser.add_argument('--data_dir', type=str, default="./data/static", help="Individual cell cluster's working directory")

args = parser.parse_args()
clus = args.cluster_number
ddir = args.data_dir

fwflag = 99 # 1 + 2 + 32 + 64
bwflag = 147 # 1 + 2 + 16 + 128
mapq = 60
rnext = '='
lshift = +4
rshift = -5
seqlen = 50
cigar = f'{seqlen}M'
seq = 'N' * seqlen
qual = 'F' * seqlen
valid_chr = [f"chr{i}" for i in range(1,23)] + ['chrX', 'chrY']
valid_chr = dict([(i,0) for i in valid_chr])
    
sam_header_string = """@HD	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr10	LN:133797422
@SQ	SN:chr11	LN:135086622
@SQ	SN:chr12	LN:133275309
@SQ	SN:chr13	LN:114364328
@SQ	SN:chr14	LN:107043718
@SQ	SN:chr15	LN:101991189
@SQ	SN:chr16	LN:90338345
@SQ	SN:chr17	LN:83257441
@SQ	SN:chr18	LN:80373285
@SQ	SN:chr19	LN:58617616
@SQ	SN:chr2	LN:242193529
@SQ	SN:chr20	LN:64444167
@SQ	SN:chr21	LN:46709983
@SQ	SN:chr22	LN:50818468
@SQ	SN:chr3	LN:198295559
@SQ	SN:chr4	LN:190214555
@SQ	SN:chr5	LN:181538259
@SQ	SN:chr6	LN:170805979
@SQ	SN:chr7	LN:159345973
@SQ	SN:chr8	LN:145138636
@SQ	SN:chr9	LN:138394717
@SQ	SN:chrX	LN:156040895
@SQ	SN:chrY	LN:57227415
"""
    
def format_sam(s):
    [chrom, srt, end, bc, rpt] = s.strip().split('\t')
    if chrom.lower() not in valid_chr:
        return    
    qname = f"{chrom}:{srt}:{end}:{bc}"
    fwpos = int(srt) - lshift + 1          # fragment is 0-index, sam is 1-index (bam is 0-index)
    bwpos = int(end) - rshift + 1 - seqlen # reverse strand, left-most position
    tlen  = bwpos + seqlen - fwpos
    for c in range(int(rpt)):
        print(f"{qname}:{c}\t{fwflag}\t{chrom}\t{fwpos}\t{mapq}\t" + 
              f"{cigar}\t{rnext}\t{bwpos}\t{tlen}\t{seq}\t{qual}\tCB:Z:{bc}\n", end="")
        print(f"{qname}:{c}\t{bwflag}\t{chrom}\t{bwpos}\t{mapq}\t" + 
              f"{cigar}\t{rnext}\t{fwpos}\t{tlen*-1}\t{seq}\t{qual}\tCB:Z:{bc}\n", end="")
    
def filter_fragment_file(clus, data_dir):
    atac_fname = os.path.join(data_dir, f"Subset{clus}", "names_atac0.txt")
    ref = pd.read_csv(atac_fname, header=None).to_numpy().flatten()
    ref = dict([(r,0) for r in ref])
    
    print(sam_header_string, end="")
        
    while True:
        ls = sys.stdin.readlines(50000)
        if len(ls) == 0:
            break
        t1 = [l.split('\t') for l in ls]
        wt = [x for x, y in zip(ls, t1) if y[3] in ref and '#' not in ref]
        for w in wt:
            format_sam(w)


filter_fragment_file(clus, data_dir=ddir)













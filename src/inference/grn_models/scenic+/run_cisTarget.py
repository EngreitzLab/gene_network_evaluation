import os
import pickle
import argparse
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path_otsu_bin_topics', required=True)
parser.add_argument('-t', '--path_top_bin_topics', required=True)
parser.add_argument('-m', '--path_markers', required=True)
parser.add_argument('-o', '--organism', required=True)
parser.add_argument('-d', '--path_motif_enrichment', required=True)

args = vars(parser.parse_args())
path_otsu_bin_topics = args['path_otsu_bin_topics']
path_top_bin_topics = args['path_top_bin_topics']
path_markers = args['path_markers']
organism = args['organism']
path_motif_enrichment = args['path_motif_enrichment']

# Load region sets
region_bin_topics_otsu = pickle.load(open(path_otsu_bin_topics, 'rb'))
region_bin_topics_top = pickle.load(open(path_top_bin_topics, 'rb'))
markers_dict = pickle.load(open(path_markers, 'rb'))
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top.keys():
    regions = region_bin_topics_top[topic].index[region_bin_topics_top[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    if len(regions) > 0:
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

# Select dbs
if organism == 'human':
    rankings_db = "resources/rankings/human/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"  # TODO: move to resources
    scores_db = "resources/scores/human/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"  # TODO: move to resources
    motif_annotation = "resources/annotations/human/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"  # TODO: move to resources 
    species = "homo_sapiens"

# Run pycisTarget
path_output = os.path.dirname(path_motif_enrichment)
run_pycistarget(
    region_sets = region_sets,
    species = species,
    save_path =  path_output,
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = "/cellar/users/aklie/tmp",  # TODO: no hardcoding
    annotation_version = 'v10nr_clust',
)

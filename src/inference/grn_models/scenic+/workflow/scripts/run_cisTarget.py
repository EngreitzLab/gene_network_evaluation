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
parser.add_argument('-a', '--path_annotation', required=True)
parser.add_argument('-s', '--path_scores', required=True)
parser.add_argument('-r', '--path_rankings', required=True)
parser.add_argument('-j','--n_cpu', required=True)
parser.add_argument('-x', '--temp_dir', required=True)
parser.add_argument('-w', '--run_without_promoters', required=False)

args = vars(parser.parse_args())
path_otsu_bin_topics = args['path_otsu_bin_topics']
path_top_bin_topics = args['path_top_bin_topics']
path_markers = args['path_markers']
organism = args['organism']
path_motif_enrichment = args['path_motif_enrichment']
path_annotation = args['path_annotation']
path_scores = args['path_scores']
path_rankings = args['path_rankings']
n_cpu = int(args['n_cpu'])
temp_dir = args['temp_dir']
run_without_promoters = bool(args['run_without_promoters'])

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
if len(region_sets['topics_otsu']) == 0:
    print('No regions found in topics_otsu')
    del region_sets['topics_otsu']
for topic in region_bin_topics_top.keys():
    regions = region_bin_topics_top[topic].index[region_bin_topics_top[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
if len(region_sets['topics_otsu']) == 0:
    print('No regions found in topics_otsu')
    del region_sets['topics_otsu']
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    if len(regions) > 0:
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
if len(region_sets['DARs']) == 0:
    print('No DARs found')
    del region_sets['DARs']

# Select dbs
if organism == 'human':
    species = "homo_sapiens"

# Run pycisTarget
path_output = os.path.dirname(path_motif_enrichment)
run_pycistarget(
    region_sets = region_sets,
    species = species,
    save_path = path_output,
    ctx_db_path = path_rankings,
    dem_db_path = path_scores,
    path_to_motif_annotations = path_annotation,
    run_without_promoters = run_without_promoters,
    n_cpu = n_cpu,
    _temp_dir = temp_dir,
    annotation_version = 'v10nr_clust',
)

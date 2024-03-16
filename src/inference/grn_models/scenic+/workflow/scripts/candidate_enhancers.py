import pickle
import argparse
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import impute_accessibility, normalize_scores, find_highly_variable_features, find_diff_features


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--path_cistopic_obj', required=True)
parser.add_argument('-o', '--path_otsu_bin_topics', required=True)
parser.add_argument('-n', '--ntop', required=True)
parser.add_argument('-t', '--path_top_bin_topics', required=True)
parser.add_argument('-s', '--scale_factor_impute', required=True)
parser.add_argument('-r', '--scale_factor_normalize', required=True)
parser.add_argument('-m', '--path_markers', required=True)

args = vars(parser.parse_args())
path_cistopic_obj = args['path_cistopic_obj']
path_otsu_bin_topics = args['path_otsu_bin_topics']
ntop = int(args['ntop'])
path_top_bin_topics = args['path_top_bin_topics']
scale_factor_impute = float(args['scale_factor_impute'])
scale_factor_normalize = float(args['scale_factor_normalize'])
path_markers = args['path_markers']

# Read cisTopic object
cistopic_obj = pickle.load(open(path_cistopic_obj, 'rb'))

# Binarize topics using otsu and top
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, 
    method='otsu',
)
region_bin_topics_top = binarize_topics(
    cistopic_obj, 
    method='ntop', 
    ntop=ntop
)

# Save regions as pickle
pickle.dump(region_bin_topics_otsu, open(path_otsu_bin_topics, 'wb'))
pickle.dump(region_bin_topics_top, open(path_top_bin_topics, 'wb'))

# Impute accessibility by matrix multiplication
imputed_acc_obj = impute_accessibility(
        cistopic_obj,  # cisTopic object
        selected_cells=None,  # certain cells to use for imputation
        selected_regions=None,  # certain regions to use for imputation
        scale_factor=scale_factor_impute,  # A number to multiply the imputed values for.
        project=cistopic_obj.project + "impute"  # A string to add to the name of the imputed object.
    )

# Normalize accessibility scores
normalized_imputed_acc_obj = normalize_scores(
        imputed_acc_obj, 
        scale_factor=scale_factor_normalize  # Scale factor for cell-level normalization
    )

# Find highly variable regions
variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj, 
        plot=False
    )
del normalized_imputed_acc_obj

# Find differentially accessible regions
markers = find_diff_features(
        cistopic_obj, 
        imputed_acc_obj, 
        variable="celltype",
        var_features=variable_regions, 
        split_pattern="_"
    )

# Save markers as pickle
pickle.dump(markers, open(path_markers, 'wb'))

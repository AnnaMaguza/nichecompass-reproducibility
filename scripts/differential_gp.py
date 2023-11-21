# Imports =================
print('imports')
import argparse
import glob
import os
import pandas as pd
from nichecompass.models import NicheCompass

# Parameters =================
print('define parameters')

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model_type')
parser.add_argument('-s', '--suffix')
parser.add_argument('-r', '--leiden_resolution', type=float)

args = parser.parse_args()
model_label = args.model_type
model_suffix = args.suffix
latent_leiden_resolution = args.leiden_resolution
latent_cluster_key = f"latent_leiden_{latent_leiden_resolution}"

# Define paths
dataset = 'nanostring_cosmx_human_nsclc'
model_folder_path = glob.glob(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/models/{model_label}/*_{model_suffix}')[0]
result_folder_path = model_folder_path.replace('models','results')
load_timestamp = model_folder_path.split('/')[-1]

## AnnData keys
gp_names_key = "nichecompass_gp_names"
structure_dict_key = 'cluster_groups'
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"

# Create required directories
os.makedirs(result_folder_path, exist_ok=True)

# Load =================
print('load model')
model = NicheCompass.load(dir_path=model_folder_path,
                          adata=None,
                          adata_file_name=f"{dataset}_{model_label}.h5ad",
                          gp_names_key=gp_names_key)
model.adata.uns['nichecompass_active_gp_names'] = model.get_active_gps()

# Run differential gp testing =================
print('differential testing')
#compare = {
#    'tumor_clusters': ['2','4','5','9','10','11','12'],
#    'stroma_clusters': ['0','1','13','14'],
#    'neutrophil_clusters': ['6','8'],
#    'macrophage_clusters': ['7','15','13']
#}

compare = model.adata.uns[structure_dict_key][latent_cluster_key]
res = {}

# Run differential gp testing
log_bayes_factor_thresh = 2.3 # 2.3 strong threshold; 4.6 decisive threshold (https://en.wikipedia.org/wiki/Bayes_factor)

# test path before execution
a = pd.DataFrame({'test':[1,2], 'train':[3,4]})
a.to_csv(f'{result_folder_path}/gpTest_test.csv')
for structure, clusters in compare.items():
    print(structure)
    if len(clusters) > 1:
        '''
        cl = clusters[0]
        print(cl)
        enriched_gps = model.run_differential_gp_tests(
            cat_key=latent_cluster_key,
            selected_cats = [cl],
            comparison_cats=[x for x in clusters if x != cl],
            log_bayes_factor_thresh=log_bayes_factor_thresh)
        # res[f'{structure}_{cl}'] = model.adata.uns['nichecompass_differential_gp_test_results'][model.adata.uns['nichecompass_differential_gp_test_results']['p_h1']<0.05]
        res[f'{structure}_{cl}'] = model.adata.uns[differential_gp_test_results_key]
        res[f'{structure}_{cl}'].to_csv(f'{result_folder_path}/gpTest_{structure}_{cl}_r{latent_leiden_resolution}.csv')
        
        '''
        for cl in clusters:
            if not os.path.exists(f'{result_folder_path}/gpTest_{structure}_{cl}_r{latent_leiden_resolution}.csv'):
                print(cl)
                enriched_gps = model.run_differential_gp_tests(
                    cat_key=latent_cluster_key,
                    selected_cats = [cl],
                    comparison_cats=[x for x in clusters if x != cl],
                    log_bayes_factor_thresh=log_bayes_factor_thresh)
                res[f'{structure}_{cl}'] = model.adata.uns[differential_gp_test_results_key]
                res[f'{structure}_{cl}'].to_csv(f'{result_folder_path}/gpTest_{structure}_{cl}_r{latent_leiden_resolution}.csv')
        





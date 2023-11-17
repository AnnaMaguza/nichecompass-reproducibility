import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import gc
import glob
from nichecompass.models import NicheCompass
from nichecompass.utils import add_gps_from_gp_dict_to_adata
from nichecompass.utils import aggregate_obsp_matrix_per_cell_type
from nichecompass.utils import compute_communication_gp_network

def load_adata(suffix, model='reference_query_mapping', dataset='nanostring_cosmx_human_nsclc'):    
    model_folder = glob.glob(f'{base_path}/artifacts/{dataset}/models/{model}/*_{suffix}')[0]
    adata_path = f"{model_folder}/{dataset}_{model}.h5ad"
    adata = sc.read_h5ad(adata_path)
    
    from matplotlib.colors import to_hex
    batch_colors = np.apply_along_axis(to_hex, 1, np.array(plt.get_cmap('tab20b').colors)[[0,1,3,5,8,11,13,17],:])
    batch_colors = {b: c for b, c in zip(['lung5_rep1','lung5_rep2','lung5_rep3','lung6','lung9_rep1','lung9_rep2','lung12','lung13'], batch_colors)}
    adata.uns['batch_colors'] = [batch_colors[b] for b in adata.obs.batch.unique()]
    
    return adata, adata_path

base_path = '/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility'


ref =  NicheCompass.load(dir_path=f'{base_path}/artifacts/nanostring_cosmx_human_nsclc/models/reference/19102023_172844_43/',
                  adata=None,
                  adata_file_name='nanostring_cosmx_human_nsclc_reference.h5ad',
                  gp_names_key='nichecompass_gp_names')

mmodel =  NicheCompass.load(dir_path=f'{base_path}/artifacts/nanostring_cosmx_human_nsclc/models/reference_query_mapping/19102023_172844_43_8/',
                  adata=None,
                  adata_file_name='nanostring_cosmx_human_nsclc_reference_query_mapping.h5ad',
                  gp_names_key='nichecompass_gp_names')

mmodel.adata.uns['nichecompass_source_genes_idx'] = ref.adata.uns['nichecompass_source_genes_idx']
mmodel.adata.uns['nichecompass_sources_categories_label_encoder'] = ref.adata.uns['nichecompass_sources_categories_label_encoder']
mmodel.adata.uns['nichecompass_target_genes_idx'] = ref.adata.uns['nichecompass_target_genes_idx']
mmodel.adata.uns['nichecompass_targets_categories_label_encoder'] = ref.adata.uns['nichecompass_targets_categories_label_encoder']
mmodel.adata.uns['nichecompass_source_genes_idx'] = ref.adata.uns['nichecompass_source_genes_idx']
mmodel.adata.varm = ref.adata.varm
del ref
gc.collect()


gp_names_key='nichecompass_gp_names'
active_gp_names_key='nichecompass_active_gp_names'
mmodel.adata.uns[active_gp_names_key] = mmodel.get_active_gps()


#nx_s13 = compute_communication_gp_network(gp_list=['Spp1_ligand_receptor_target_gene_GP'],  model=mmodel, group_key="latent_leiden_0.7", filter_key='batch', filter_cat='lung13')
#nx_s13.to_csv(f'{base_path}/artifacts/nanostring_cosmx_human_nsclc/results/reference_query_mapping/19102023_172844_43_8/Spp1_ligand_receptor_target_gene_GP_lung13.csv')

nx_s6 = compute_communication_gp_network(gp_list=['Spp1_ligand_receptor_target_gene_GP'],  model=mmodel, group_key="latent_leiden_0.7", filter_key='batch', filter_cat='lung6')
nx_s6.to_csv(f'{base_path}/artifacts/nanostring_cosmx_human_nsclc/results/reference_query_mapping/19102023_172844_43_8/Spp1_ligand_receptor_target_gene_GP_lung6.csv')
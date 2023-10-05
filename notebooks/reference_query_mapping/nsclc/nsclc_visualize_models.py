import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import glob
import os
from matplotlib.colors import to_hex

dataset = 'nanostring_cosmx_human_nsclc'
model_type = 'reference_query_mapping'
available_models = glob.glob(f'artifacts/{dataset}/models/{model_type}/*')
for model in available_models:
    figure_path = model.replace('models','figures')
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)
    if os.path.exists(f"{model}/{dataset}_{model_type}.h5ad") and not os.path.exists(f"{figure_path}/umap_all.png"):
        adata = sc.read_h5ad(f"{model}/{dataset}_{model_type}.h5ad")
        adata.obs_names = adata.obs['cell_ID'].values
        adata.obs['cell type'] = adata.obs['cell_type_original'].map({
            'tumor 9': 'tumor', 'tumor 6': 'tumor', 'tumor 5': 'tumor', 'tumor 13': 'tumor', 'tumor 12': 'tumor',
            'fibroblast': 'fibroblast', 'neutrophil': 'neutrophil', 'plasmablast': 'plasmablast', 'B-cell': 'B-cell', 'endothelial': 'endothelial', 'epithelial': 'epithelial',
            'T CD4 naive': 'T cell/NK', 'T CD4 memory': 'T cell/NK', 'T CD8 naive': 'T cell/NK', 'T CD8 memory': 'T cell/NK', 'Treg': 'T cell/NK', 'NK': 'T cell/NK',
            'mDC': 'DC/monocyte', 'pDC': 'DC/monocyte', 'monocyte': 'DC/monocyte', 'macrophage': 'macrophage', 'mast': 'mast'
        }).astype("category")
        adata.uns['cell type_colors'] = np.apply_along_axis(to_hex, 1, np.array(plt.get_cmap('tab20').colors))[[6,10,2,14,4,5,0,1,8,12,18]]
        adata.uns['cell type_colors'][3] = '#ffff33'
        adata.uns['mapping_entity_colors'] = ['#00AFBB', '#999999']
        adata.obs['cell type'] = adata.obs['cell type'].cat.reorder_categories(['tumor','fibroblast','endothelial','epithelial','macrophage','DC/monocyte','B-cell','plasmablast', 'T cell/NK', 'mast', 'neutrophil'], ordered=False)
        sc.pl.umap(adata, color=['batch','cell type','niche','mapping_entity'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_all.png", bbox_inches="tight")
        sc.pl.umap(adata[adata.obs.mapping_entity == 'reference'], color=['batch','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_reference.png", bbox_inches="tight")
        new_batches = adata.obs.loc[adata.obs.mapping_entity == 'query','batch'].unique()
        if len(new_batches) == 2:
            sc.pl.umap(adata[adata.obs.batch != new_batches[1]], color=['mapping_entity','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
            plt.savefig(f"{figure_path}/umap_reference_query{new_batches[0]}.png", bbox_inches="tight")
            sc.pl.umap(adata[adata.obs.batch != new_batches[0]], color=['mapping_entity','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
            plt.savefig(f"{figure_path}/umap_reference_query{new_batches[1]}.png", bbox_inches="tight")


model_type = 'reference'
available_models = glob.glob(f'artifacts/{dataset}/models/{model_type}/*')

for model in available_models:
    figure_path = model.replace('models','figures')
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)
    if os.path.exists(f"{model}/{dataset}_{model_type}.h5ad") and not os.path.exists(f"{figure_path}/umap_all.png"):
        adata = sc.read_h5ad(f"{model}/{dataset}_{model_type}.h5ad")
        adata.obs_names = adata.obs['cell_ID'].values
        adata.obs['cell type'] = adata.obs['cell_type_original'].map({
            'tumor 9': 'tumor', 'tumor 6': 'tumor', 'tumor 5': 'tumor', 'tumor 13': 'tumor', 'tumor 12': 'tumor',
            'fibroblast': 'fibroblast', 'neutrophil': 'neutrophil', 'plasmablast': 'plasmablast', 'B-cell': 'B-cell', 'endothelial': 'endothelial', 'epithelial': 'epithelial',
            'T CD4 naive': 'T cell/NK', 'T CD4 memory': 'T cell/NK', 'T CD8 naive': 'T cell/NK', 'T CD8 memory': 'T cell/NK', 'Treg': 'T cell/NK', 'NK': 'T cell/NK',
            'mDC': 'DC/monocyte', 'pDC': 'DC/monocyte', 'monocyte': 'DC/monocyte', 'macrophage': 'macrophage', 'mast': 'mast'
        }).astype("category")
        adata.uns['cell type_colors'] = np.apply_along_axis(to_hex, 1, np.array(plt.get_cmap('tab20').colors))[[6,10,2,14,4,5,0,1,8,12,18]]
        adata.uns['cell type_colors'][3] = '#ffff33'
        adata.uns['mapping_entity_colors'] = ['#00AFBB', '#999999']
        adata.obs['cell type'] = adata.obs['cell type'].cat.reorder_categories(['tumor','fibroblast','endothelial','epithelial','macrophage','DC/monocyte','B-cell','plasmablast', 'T cell/NK', 'mast', 'neutrophil'], ordered=False)
        sc.pl.umap(adata, color=['batch','cell type','niche','mapping_entity'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_all.png", bbox_inches="tight")


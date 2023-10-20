import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import glob
import gc
import argparse
from nichecompass.models import NicheCompass
from matplotlib.colors import to_hex

def load_model(suffix, cluster=False, latent_leiden_resolution=0.5):
    model_folder = glob.glob(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/models/{model}/*_{suffix}')[0]
    figure_path = model_folder.replace('models','figures')
    adata = sc.read_h5ad(f"{model_folder}/{dataset}_{model}.h5ad")
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

    if model == 'reference_query_mapping':
        new_batches = adata.obs.loc[adata.obs.mapping_entity == 'query','batch'].unique()
        sc.pl.umap(adata[adata.obs.mapping_entity == 'reference'], color=['batch','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_reference.png", bbox_inches="tight")
        sc.pl.umap(adata[adata.obs.batch != new_batches[1]], color=['mapping_entity','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_reference_query_{new_batches[0]}.png", bbox_inches="tight")
        sc.pl.umap(adata[adata.obs.batch != new_batches[0]], color=['mapping_entity','cell type','niche'], ncols=4, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_reference_query_{new_batches[1]}.png", bbox_inches="tight")

    if cluster:
        if not f"latent_leiden_{latent_leiden_resolution}" in adata.obs.columns:
            sc.tl.leiden(adata, resolution=latent_leiden_resolution, key_added=f"latent_leiden_{latent_leiden_resolution}", neighbors_key="nichecompass_latent")
            adata.write_h5ad(f"{model_folder}/{dataset}_{model}.h5ad") # re-write file with new metadata  

        sc.pl.umap(adata, color=['batch', f'latent_leiden_{latent_leiden_resolution}','cell type','niche','mapping_entity'], ncols=5, wspace=0.5, size=0.5, show=False)
        plt.savefig(f"{figure_path}/umap_all.png", bbox_inches="tight")        
        

        # Plot sample 5 for fov effects
        for b in adata.obs.loc[adata.obs['patient']=='Lung5','batch'].unique():
            sc.pl.embedding(adata[adata.obs.batch==b], basis="spatial", color=['cell type', f"latent_leiden_{latent_leiden_resolution}"], ncols=2, wspace=0.3, size=3, show=False)
            plt.savefig(f"{figure_path}/spatial_{b}.png", bbox_inches="tight")    
    
    adata = NicheCompass.load(dir_path=model_folder,
                          adata=adata,
                          gp_names_key='nichecompass_gp_names')
    
    # check number of active gps
    active_gps = adata.get_active_gps()
    print(f"{len(active_gps)} active GPs out of {len(adata.adata.uns['nichecompass_gp_names'])} total GPs.")
    
    gc.collect()
    return



parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model_type')
parser.add_argument('-s', '--suffix')
parser.add_argument('-c', '--cluster', action='store_true')
parser.add_argument('-r', '--leiden_resolution', type=float)

args = parser.parse_args()

dataset = 'nanostring_cosmx_human_nsclc'
model = args.model_type
load_model(args.suffix, cluster=args.cluster, latent_leiden_resolution=args.leiden_resolution)

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model')

model_label = parser.parse_args().model
dataset = 'nanostring_cosmx_human_nsclc'
model_type = 'reference'
cell_type_key='cell_type'
batch_key='batch'
latent_key='nichecompass_latent'

# load model
adata = sc.read_h5ad(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/models/{model_type}/{model_label}/{dataset}_{model_type}.h5ad')
adata.obs_names = adata.obs['cell_ID'].values
adata.obs[cell_type_key] = adata.obs['cell_type_original'].map({
    'tumor 9': 'tumor', 'tumor 6': 'tumor', 'tumor 5': 'tumor', 'tumor 13': 'tumor', 'tumor 12': 'tumor',
    'fibroblast': 'fibroblast', 'neutrophil': 'neutrophil', 'plasmablast': 'plasmablast', 'B-cell': 'B-cell', 'endothelial': 'endothelial', 'epithelial': 'epithelial',
    'T CD4 naive': 'T cell/NK', 'T CD4 memory': 'T cell/NK', 'T CD8 naive': 'T cell/NK', 'T CD8 memory': 'T cell/NK', 'Treg': 'T cell/NK', 'NK': 'T cell/NK',
    'mDC': 'DC/monocyte', 'pDC': 'DC/monocyte', 'monocyte': 'DC/monocyte', 'macrophage': 'macrophage', 'mast': 'mast'
}).astype("category")

# Run scBI benchmark
bm = Benchmarker(
    adata,
    batch_key=batch_key,
    label_key=cell_type_key,
    embedding_obsm_keys=[latent_key],
    bio_conservation_metrics=BioConservation(isolated_labels=False, nmi_ari_cluster_labels_leiden=False, nmi_ari_cluster_labels_kmeans=False, silhouette_label=True, clisi_knn=True),
    batch_correction_metrics=BatchCorrection(silhouette_batch=True, ilisi_knn=True, kbet_per_label=True, graph_connectivity=False, pcr_comparison=True),
    n_jobs=20
)
bm.benchmark()
df = bm.get_results(min_max_scale=False)
df.to_csv(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/results/{model_type}/{model_label}/scbi_metrics_m.csv')

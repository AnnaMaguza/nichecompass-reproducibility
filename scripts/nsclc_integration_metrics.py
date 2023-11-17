from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model')
parser.add_argument('-b', '--batch')
parser.add_argument('-l', '--latent')
parser.add_argument('-wS', '--withinSample', action='store_true')

model_label = parser.parse_args().model
batch_key = parser.parse_args().batch
latent_key = parser.parse_args().latent
within_sample = parser.parse_args().withinSample

dataset = 'nanostring_cosmx_human_nsclc'
model_type = 'reference_query_mapping'
cell_type_key='cell_type'
# batch_key='mapping_entity'
# latent_keys=['nichecompass_latent','X_pca']

# load model
adata = sc.read_h5ad(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/models/{model_type}/{model_label}/{dataset}_{model_type}.h5ad')
adata.obs_names = adata.obs['cell_ID'].values
adata.obs[cell_type_key] = adata.obs['cell_type_original'].map({
    'tumor 9': 'tumor', 'tumor 6': 'tumor', 'tumor 5': 'tumor', 'tumor 13': 'tumor', 'tumor 12': 'tumor',
    'fibroblast': 'fibroblast', 'neutrophil': 'neutrophil', 'plasmablast': 'plasmablast', 'B-cell': 'B-cell', 'endothelial': 'endothelial', 'epithelial': 'epithelial',
    'T CD4 naive': 'T cell/NK', 'T CD4 memory': 'T cell/NK', 'T CD8 naive': 'T cell/NK', 'T CD8 memory': 'T cell/NK', 'Treg': 'T cell/NK', 'NK': 'T cell/NK',
    'mDC': 'DC/monocyte', 'pDC': 'DC/monocyte', 'monocyte': 'DC/monocyte', 'macrophage': 'macrophage', 'mast': 'mast'
}).astype("category")


# keep only those from the same replicate as the query
if within_sample:
    sample_name = adata.obs.loc[adata.obs.mapping_entity == 'query','patient'].unique()[0]
    batches = adata.obs.loc[adata.obs.patient == sample_name,'batch'].unique()
    if len(batches) > 1:
        adata = adata[adata.obs.patient == sample_name]
    

# Run scBI benchmark
if latent_key == 'X_pca':
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
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
df.to_csv(f'/lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/artifacts/{dataset}/results/{model_type}/{model_label}/scbi_metrics_{latent_key}_{batch_key}.csv')

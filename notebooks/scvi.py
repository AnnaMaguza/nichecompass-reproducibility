import argparse
import scanpy as sc
import scvi
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-cvCat', '--categorical_covariate', action='append', default=None)
args = parser.parse_args()

# File name
covariates = args.categorical_covariate 
covariates = '-'.join(covariates)

# Create anndata
adata = sc.read_h5ad("datasets/srt_data/gold/nanostring_cosmx_human_nsclc_all_raw.h5ad")

# Train model
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=covariates)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
latent=vae.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.write("datasets/srt_data/gold/nanostring_cosmx_human_nsclc_all_raw_scVI.h5ad")
latent = pd.DataFrame(latent, index=adata.obs_names, columns=[f'scvi_{i+1}' for i in range(latent.shape[1])])
latent.to_csv("datasets/srt_data/irene_temp/nanostring_cosmx_human_nsclc/scVI_latent.csv")


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import rcParams
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_class_weight

def preprocess_subsample(adata, samples, n_neighbours=10, n_pcs=30, exclude_cells={}, plot_covariates=[]):
    # Subsample
    subadata = adata[adata.obs['sample'].isin(samples)].copy()
    for key, value in exclude_cells.items():
        subadata = subadata[~subadata.obs[key].isin(value)].copy()
    subadata.raw = subadata
    print(subadata.shape)
    
    # PCA, UMAP, Cluster
    sc.pp.scale(subadata)
    sc.tl.pca(subadata, svd_solver='arpack')
    rcParams['figure.figsize'] = (5, 3)
    sc.pl.pca_variance_ratio(subadata, log=True)
    sc.pp.neighbors(subadata, n_neighbors=n_neighbours, n_pcs=n_pcs)
    sc.tl.umap(subadata)
    sc.tl.leiden(subadata)
    
    # Some UMAP plots
    plot_covariates.append('leiden')
    rcParams['figure.figsize'] = (5, 3)
    sc.pl.umap(subadata, color=plot_covariates, ncols=2)
    
    return(subadata)

def cell_signatures(adata, cell_label_column, ratio_threshold=2):
    rcParams['figure.figsize'] = (5, 3)
    sc.tl.rank_genes_groups(adata, cell_label_column, method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    
    result = adata.uns['rank_genes_groups']
    deg_all = pd.DataFrame({})
    for cell in adata.obs[cell_label_column].unique():
        deg = pd.DataFrame({'pvals_adj': result['pvals_adj'][cell], 'logFC': result['logfoldchanges'][cell]},
                          index=result['names'][cell])
        deg.loc[(deg.pvals_adj > 0.05) | (deg.logFC < 0), 'pvals_adj'] = 1
        deg.sort_values(by='pvals_adj', inplace=True)
        deg.rename(columns={'pvals_adj': f'{cell}_adjPval'}, inplace=True)
        deg.drop(columns=['logFC'], inplace=True)
        deg_all = pd.merge(deg_all, deg, how="outer", left_index=True, right_index=True)

    total_signif = (deg_all != 1.0).sum(axis=1)
    total_signif.hist(bins=50)
    plt.show()
    
    # 0 pvalues to min pvalue
    deg_all[deg_all == 0] = deg_all[deg_all != 0].min().min()
    
    # if a gene is differential in more than 1 cell type, keep only those that are more relevant
    # logFC(pval) > min(logFC(pval))/2
    for gene in deg_all.index.values:
        if total_signif[gene] > 1:
            logpvals = -np.log10(deg_all.loc[gene,:])
            threshold = max(logpvals)/ratio_threshold
            deg_all.loc[gene,logpvals < threshold] = 1.0

    (deg_all != 1.0).sum(axis=1).hist(bins=50)
    plt.show()
    print((deg_all != 1.0).sum(axis=0))
    
    markers = {}
    for col in deg_all.columns:
        cell = col.replace('_adjPval','')
        markers[cell] = list(deg_all.index.values[deg_all[col] < 0.05])

    return(adata, markers)

def score_cells(adata, markers, cell_label_column):
    score_cols = []
    for cell in markers:
        if len(markers[cell]) > 0:
            sc.tl.score_genes(adata, markers[cell], score_name=f'{cell}_score')
            score_cols.append(f'{cell}_score')

    rcParams['figure.figsize'] = (5, 4)
    sc.pl.umap(adata, color=cell_label_column, ncols=1)
    rcParams['figure.figsize'] = (3, 2)
    sc.pl.umap(adata, color=score_cols, ncols=5)
    
    adata.obs['cell_type_predicted'] = adata.obs.loc[:,adata.obs.columns.str.contains('_score')].idxmax(axis=1).str.replace('_score','').values
    adata.obs['cell_type_predicted_ok'] = adata.obs['cell_type_predicted'] == adata.obs[cell_label_column]
    print(adata.obs['cell_type_predicted_ok'].value_counts())
    print(adata.obs.groupby(cell_label_column)['cell_type_predicted_ok'].value_counts())
    return(adata)

def plot_cluster_proportions(cluster_props, 
                             cluster_palette=None,
                             xlabel_rotation=0): 
    fig, ax = plt.subplots(dpi=150)
    fig.patch.set_facecolor("white")
    
    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title="Cluster")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout()
    
    return fig

def predict(adata_test, models, scaler):
    
    X_test = np.array(adata_test.raw.X.todense())
    X_test = scaler.transform(X_test)
    
    for model_name, clf in models.items():
        y_hat_test = clf.predict(X_test)
        adata_test.obs[model_name] = y_hat_test
        sc.pl.umap(adata_test, color=model_name, ncols=1)
    
    return adata_test
    
def train_classifier(adata_train, y_col, model_config, train_pct=0.8):
    model_store = {}
    
    # Exclude unknown
    adata_train = adata_train[adata_train.obs[y_col] != 'Unknown',:]
    
    # train-val , X-y split
    val_id = int(np.trunc(adata_train.shape[0] * train_pct))
    X_train = np.array(adata_train.raw.X[0:val_id,:].todense())
    y_train = adata_train.obs[y_col].astype(str).values[0:val_id]
    X_val = np.array(adata_train.raw.X[val_id:,:].todense())
    y_val = adata_train.obs[y_col].astype(str).values[val_id:]
    
    # Train models
    for model in model_config:
        print(model['name'])

        # Scale
        scaler = StandardScaler().fit(X_train)
        X_train = scaler.transform(X_train)
        X_val = scaler.transform(X_val)    
        
        # Fit model
        if model['weight_samples']:
            weight_dict = compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
            weight_dict = {cell_type: weight for cell_type, weight in zip(np.unique(y_train), weight_dict)}
            sample_weight = np.array([weight_dict[cell_type] for cell_type in y_train])
            clf = model['model'](**model['parameters']).fit(X_train, y_train, sample_weight)
            
        else:
            clf = model['model'](**model['parameters']).fit(X_train, y_train)
        
        model_store[f'{y_col}_{model["name"]}'] = clf
        
        # Predict
        y_hat_train = clf.predict(X_train)
        y_hat_val = clf.predict(X_val)    

        # Metrics
        print(f"Train F1 score: {metrics.f1_score(y_train, y_hat_train, average='weighted')}")  
        print(f"Validation F1 score: {metrics.f1_score(y_val, y_hat_val, average='weighted')}")      
        print("Train confuison matrix")    
        contingency = pd.DataFrame(metrics.confusion_matrix(y_train, y_hat_train), index=np.unique(y_val), columns=np.unique(y_val))
        contingency = contingency.div(contingency.sum(axis=1), axis=0)
        sns.clustermap(contingency, cmap='Blues', row_cluster=False, col_cluster=False, figsize=(7, 7))
        plt.show()
        print("Validation confusion matrix")
        contingency = pd.DataFrame(metrics.confusion_matrix(y_val, y_hat_val), index=np.unique(y_val), columns=np.unique(y_val))
        contingency = contingency.div(contingency.sum(axis=1), axis=0)
        sns.clustermap(contingency, cmap='Blues', row_cluster=False, col_cluster=False, figsize=(7, 7))
        plt.show()
        
        # UMAPs
        rcParams['figure.figsize'] = (5, 4)
        adata_train.obs[f'{y_col}_{model["name"]}'] = np.concatenate((y_hat_train, y_hat_val))
        sc.pl.umap(adata_train, color=f'{y_col}_{model["name"]}', ncols=1)
        
    return adata_train, model_store, scaler
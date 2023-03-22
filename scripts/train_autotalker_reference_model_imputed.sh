python train_autotalker_reference_model.py \
--dataset starmap_plus_mouse_cns \
--reference_batches batch1_imputed batch2_imputed batch3_imputed \
--n_neighbors 12 \
--filter_genes \
--n_hvg 4000 \
--nichenet_max_n_target_genes_per_gp 100 \
--include_mebocost_gps \
--mebocost_species mouse \
--counts_key None \
--n_cond_embed 60 \
--no-log_variational \
--n_epochs 40 \
--n_epochs_all_gps 20 \
--lambda_group_lasso 0. \
--lambda_l1_masked 0. \
--edge_batch_size 128 \
--node_batch_size 16

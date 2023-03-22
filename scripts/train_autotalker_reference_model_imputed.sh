python train_autotalker_reference_model.py \
--dataset starmap_plus_mouse_cns \
--reference_batches batch1_imputed \
--n_neighbors 4 \
--filter_genes \
--n_hvg 4000 \
--nichenet_max_n_target_genes_per_gp 100 \
--include_mebocost_gps \
--mebocost_species mouse \
--counts_key None \
--no-log_variational \
--n_epochs 1 \
--n_epochs_all_gps 1 \
--lambda_group_lasso 0.01 \
--lambda_l1_masked 0.01 \
--edge_batch_size 256 \
--node_batch_size 32

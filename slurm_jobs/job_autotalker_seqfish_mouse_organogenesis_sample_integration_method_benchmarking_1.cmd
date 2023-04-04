#!/bin/bash
#SBATCH -J autotalker_seqfish_mouse_organogenesis_sample_integration_method_benchmarking_1
#SBATCH -o /home/sbirk/projects/autotalker-reproducibility/slurm_jobs/logs/out_autotalker_seqfish_mouse_organogenesis_sample_integration_method_benchmarking_1.txt
#SBATCH -e /home/sbirk/projects/autotalker-reproducibility/slurm_jobs/logs/err_autotalker_seqfish_mouse_organogenesis_sample_integration_method_benchmarking_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate autotalker
cd /
cd /home/sbirk/projects/autotalker-reproducibility/scripts
python train_autotalker_benchmarking_models.py --adata_new_name None  --n_neighbors_list 4 4 8 8 12 12 16 16 20 20  --edge_batch_size_list 512 512 256 256 256 256 128 128 128 128  --node_batch_size_list 64 64 32 32 32 32 16 16 16 16  --seeds 0 1 2 3 4 5 6 7 8 9  --run_index 1 2 3 4 5 6 7 8 9 10  --cell_type_key celltype_mapped_refined  --nichenet_keep_target_genes_ratio 0.01  --nichenet_max_n_target_genes_per_gp 25344  --include_mebocost_gps  --mebocost_species mouse  --gp_filter_mode subset  --combine_overlap_gps  --overlap_thresh_source_genes 0.9  --overlap_thresh_target_genes 0.9  --overlap_thresh_genes 0.9  --dataset seqfish_mouse_organogenesis  --reference_batches batch1 batch2 batch3 batch4 batch5 batch6  --counts_key counts  --condition_key batch  --spatial_key spatial  --adj_key spatial_connectivities  --mapping_entity_key mapping_entity  --no-filter_genes  --gp_targets_mask_key autotalker_gp_targets  --gp_sources_mask_key autotalker_gp_sources  --gp_names_key autotalker_gp_names  --model_label sample_integration_method_benchmarking  --active_gp_names_key autotalker_active_gp_names  --latent_key autotalker_latent  --active_gp_thresh_ratio 0.03  --gene_expr_recon_dist nb  --cond_embed_injection gene_expr_decoder  --log_variational  --n_layers_encoder 1  --conv_layer_encoder gcnconv  --n_epochs 40  --n_epochs_all_gps 20  --lr 0.001  --lambda_edge_recon 10.  --lambda_gene_expr_recon 0.01  --lambda_cond_contrastive 10.  --contrastive_logits_ratio 0.1  --lambda_group_lasso 0.  --lambda_l1_masked 0. 

#!/bin/bash
#SBATCH -J staci_single_sample_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_1pct_batch5_metrics_computation_1
#SBATCH -o ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/out_staci_single_sample_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_1pct_batch5_metrics_computation_1.txt
#SBATCH -e ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/err_staci_single_sample_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_1pct_batch5_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=160G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset nanostring_cosmx_human_nsclc_subsample_1pct_batch5 --task single_sample_method_benchmarking --file_name nanostring_cosmx_human_nsclc_subsample_1pct_batch5_staci.h5ad --cell_type_key cell_type --batch_key None --latent_key staci_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi

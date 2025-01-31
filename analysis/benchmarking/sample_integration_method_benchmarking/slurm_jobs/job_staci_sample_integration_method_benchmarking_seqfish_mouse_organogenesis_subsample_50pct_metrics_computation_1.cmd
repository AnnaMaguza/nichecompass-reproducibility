#!/bin/bash
#SBATCH -J staci_sample_integration_method_benchmarking_seqfish_mouse_organogenesis_subsample_50pct_metrics_computation_1
#SBATCH -o ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/out_staci_sample_integration_method_benchmarking_seqfish_mouse_organogenesis_subsample_50pct_metrics_computation_1.txt
#SBATCH -e ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/err_staci_sample_integration_method_benchmarking_seqfish_mouse_organogenesis_subsample_50pct_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=128G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/sample_integration_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset seqfish_mouse_organogenesis_subsample_50pct --task sample_integration_method_benchmarking --file_name seqfish_mouse_organogenesis_subsample_50pct_staci.h5ad --cell_type_key cell_type --batch_key batch --batches batch1 batch2 batch3 batch4 batch5 batch6 --latent_key staci_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi basw blisi kbet pcr

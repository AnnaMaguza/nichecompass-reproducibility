#!/bin/bash
#SBATCH -J spatial_atac_rna_seq_mouse_brain_nichecompass_one-hop-norm_reference_model
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference/slurm_jobs/logs/out_spatial_atac_rna_seq_mouse_brain_nichecompass_one-hop-norm_reference_model.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference/slurm_jobs/logs/err_spatial_atac_rna_seq_mouse_brain_nichecompass_one-hop-norm_reference_model.txt
#SBATCH -t 12:00:00
#SBATCH -p interactive_gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive_gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference
bash train_spatial_atac_rna_seq_mouse_brain_nichecompass_one-hop-norm_reference_model.sh

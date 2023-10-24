#!/bin/bash
#SBATCH -J 49_3_cluster_latent
#SBATCH -o logs/out_49_3_cluster_latent.txt
#SBATCH -e logs/err_49_3_cluster_latent.txt
#SBATCH -t 24:00:00
#SBATCH -p gpu_p
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=64G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference_query/slurm_jobs
python ../../nsclc_visualize_and_cluster_models.py --model_type reference_query_mapping --suffix 49_3 --cluster --leiden_resolution 0.5

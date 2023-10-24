#!/bin/bash
#SBATCH -J 27_3_cluster_latent
#SBATCH -o logs/out_27_3_cluster_latent.txt
#SBATCH -e logs/err_27_3_cluster_latent.txt
#SBATCH -t 24:00:00
#SBATCH -p gpu_p
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --qos=gpu_normal
#SBATCH --mem=64G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference_query
python ../nsclc_visualize_and_cluster_models.py --model_type reference_query_mapping --suffix 27_3 --cluster --leiden_resolution 0.5

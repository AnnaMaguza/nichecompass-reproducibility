#!/bin/bash
#SBATCH -J 25_cluster_latent
#SBATCH -o logs/out_25_cluster_latent.txt
#SBATCH -e logs/err_25_cluster_latent.txt
#SBATCH -t 24:00:00
#SBATCH -p cpu_p
#SBATCH -c 8
#SBATCH --qos=cpu_normal
#SBATCH --mem=64G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference
python ../nsclc_visualize_and_cluster_models.py --model_type reference --suffix 25 --cluster --leiden_resolution 0.5

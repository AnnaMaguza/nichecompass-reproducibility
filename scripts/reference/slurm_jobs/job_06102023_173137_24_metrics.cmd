#!/bin/bash
#SBATCH -J 06102023_173137_24_metrics
#SBATCH -o logs/out_06102023_173137_24_metrics.txt
#SBATCH -e logs/err_06102023_173137_24_metrics.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 20
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=160G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference
python ../nsclc_integration_metrics.py --model 06102023_173137_24

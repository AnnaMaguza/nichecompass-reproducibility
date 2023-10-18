#!/bin/bash
#SBATCH -J 05102023_213300_20_metrics
#SBATCH -o logs/out_05102023_213300_20_metrics.txt
#SBATCH -e logs/err_05102023_213300_20_metrics.txt
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
python ../nsclc_integration_metrics.py --model 05102023_213300_20

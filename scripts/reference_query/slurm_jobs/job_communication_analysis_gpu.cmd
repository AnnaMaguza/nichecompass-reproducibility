#!/bin/bash
#SBATCH -J gcommunication
#SBATCH -o out_gcommunication.txt
#SBATCH -e err_gcommunication.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 20
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_long
#SBATCH --mem=240G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference_query
python ../communication.py
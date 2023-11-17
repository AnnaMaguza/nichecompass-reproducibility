#!/bin/bash
#SBATCH -J communication
#SBATCH -o out_communication.txt
#SBATCH -e err_communication.txt
#SBATCH -t 48:00:00
#SBATCH -p cpu_p
#SBATCH -c 20
#SBATCH --qos=cpu_long
#SBATCH --mem=400G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference_query
python ../communication.py

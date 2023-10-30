#!/bin/bash
#SBATCH -J 43_8_gp_test
#SBATCH -o logs/out_43_8_gp_test.txt
#SBATCH -e logs/err_43_8_gp_test.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=64G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference_query/slurm_jobs
python ../../differential_gp.py --model_type reference_query_mapping --suffix 43_8 --leiden_resolution 0.5

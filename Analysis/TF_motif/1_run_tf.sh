#!/bin/bash 
#SBATCH --job-name=DL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --nodelist=gpu02


export PATH=/public/home/guogjgroup/ggj/anaconda3/envs/pytorch/bin/:$PATH
python TF_motif_insert.py  \
    --data Humanbrain.sample.h5 \
    --model_weight $model_dir/best_model.pth.tar \
    --tf_file tf_name_id.csv \

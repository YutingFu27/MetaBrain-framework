#!/bin/bash 
#SBATCH --job-name=MBRS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
export PATH=~/anaconda3/envs/pytorch/bin/:$PATH
python MBRS.py \
    --data Humanbrain.sample.h5 \
    --model_weight $model_dir/best_model.pth.tar \
    --caqtl_file caqtl.csv \
    --genome_file Ref/hg38.fa \
    --anno_file Humanbrain.sample.anno.csv



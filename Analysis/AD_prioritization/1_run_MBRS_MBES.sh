#!/bin/bash 
#SBATCH --job-name=DL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
python gwas_MBRS.py \
    --data Humanbrain.sample.h5 \
    --model_weight $model_dir/best_model.pth.tar \
    --gwas_file gwas.csv \
    --genome_file Ref/hg19.fa \
    --anno_file brain_anno.csv 


python gwas_MBES.py \
    --data Humanbrain.sample.h5 \
    --model_weight $model_dir/best_model.pth.tar \
    --vcf_file predicted.vcf \
    --genome_file Ref/hg19.fa \
    --modellist_file modelList.txt


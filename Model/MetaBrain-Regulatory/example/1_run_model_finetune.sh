#!/bin/bash 
#SBATCH --job-name=DL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --nodelist=gpu01
export PATH= ~/anaconda3/envs/pytorch/bin/:$PATH
pretrained_dir= ~/FYT/AI/CodeTest/20240411_beluga/train/20240825_train/model_output/
python Finetune.py \
    --data Humanbrain.new.h5 \
    --lr 1e-4 \
    --patience 10 \
    --batch_size 250 \
    --pretrained_weight $pretrained_dir/best_model.pth.tar

model_dir=./model_output/
python eva.py \
    --data Humanbrain.new.h5 \
    --batch_size 250 \
    --model_weight $model_dir/best_model.pth.tar
    
python cellembed.py --model_weight $model_dir/best_model.pth.tar



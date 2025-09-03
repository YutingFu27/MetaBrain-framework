#!/bin/bash 
#SBATCH --job-name=DL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --nodelist=gpu01

#python Check_gpu.py
export PATH=/public/home/guogjgroup/ggj/anaconda3/envs/pytorch/bin/:$PATH
model_dir=./model_output/
python chromatin.py \
    --data Humanbrain.sample.h5 \
    --model_weight $model_dir/best_model.pth.tar
python merge.py

if [ -f "Xreducedall_all.npy" ]; then
  rm *_chromatin.npy
fi



import h5py, os, argparse, logging, time
import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

# args
parser = argparse.ArgumentParser()
parser.add_argument("--model_weight")
args = parser.parse_args()
model_weight=args.model_weight
load_params = torch.load(model_weight,map_location='cpu')['state_dict']
print(load_params.keys())
cell_embed=load_params['classifier.4.weight']
cell_embed.shape

w=cell_embed.data.cpu().numpy()
w.shape

np.save('cellembed.npy',w)

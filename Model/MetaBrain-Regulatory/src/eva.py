import os
import sys
import numpy as np
import pandas as pd
import glob
import argparse
import logging

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
import torch.optim as optim
import time
from time import perf_counter as pc
from tqdm import tqdm

from utils import *
sys.path.append('./arch/')
from arch import deepsea,MetaBrain,ResNet

# args
parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--model_weight")
parser.add_argument("--batch_size", dest="batch_size", default=250, type=int)
args = parser.parse_args()

#input
import h5py
h5file = h5py.File(args.data, 'r')
X = h5file["pmat"]["X"][:].swapaxes(-1,1).astype(np.float32)
peak_idx = h5file['pmat']['pmat_sc']['i'][:]
cell_idx = h5file['pmat']['pmat_sc']['j'][:]
x = h5file['pmat']['pmat_sc']['x'][:]
dim = h5file['pmat']['pmat_sc']['dim'][:]
y = np.zeros((dim[0], dim[1]), dtype = np.float32)
y[peak_idx, cell_idx] = x
features = h5file["pmat"]["peak"][:]
h5file.close()


# unpack anno
n_tasks = y.shape[-1]
mask = features[:,-1].astype(str)
test_idx = mask=='test'
x_test = X[test_idx]
y_test = y[test_idx]


# define data loader
batch_size = args.batch_size
test_loader = DataLoader(list(zip(x_test, y_test)), batch_size=batch_size, 
                            shuffle=False, num_workers=0, drop_last=False, pin_memory=True)


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
from MetaBrain import *
model = MetaBrain(num_target=n_tasks, load_base=False)
model.to(device)

model_weight=args.model_weight
load_params = torch.load(model_weight)['state_dict']
model_params = model.state_dict()

# extract the corresponding params
co_params = {k: v for k, v in load_params.items() if k in model_params.keys()}

model_params.update(co_params)
model.load_state_dict(model_params)

model.eval()
model.to(device)
target_list=[]
pred_list=[]
for idx, data in enumerate(tqdm(test_loader, 0)):        
    inputs, targets = data
    inputs, targets= inputs.to(device), targets.to(device)
            
        # testing
    with torch.no_grad():
        outputs= model.forward(inputs)        
        # append target and output
    target_list.append(targets.cpu().detach())
    pred_list.append(outputs.cpu().detach())
    
test_targets = torch.cat(target_list, dim=0) 
test_predictions = torch.cat(pred_list, dim=0)  



roc_auc_cell = calculate_roc_cell(test_targets, test_predictions)
auroc_cell = [roc_auc_cell[k] for k in roc_auc_cell.keys()]

roc_auc_peak = calculate_roc_peak(test_targets, test_predictions)
auroc_peak = [roc_auc_peak[k] for k in roc_auc_peak.keys()]




pd.DataFrame({"auroc_cell":auroc_cell}).T.to_csv("./Metric_cell.csv")

pd.DataFrame({"auroc_peak":auroc_peak}).T.to_csv("./Metric_peak.csv")

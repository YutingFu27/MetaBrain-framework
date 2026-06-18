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
parser.add_argument("--pretrained_weight")
parser.add_argument("--lr", dest="lr", default=1e-4, type=float)
parser.add_argument("--batch_size", dest="batch_size", default=250, type=int)
parser.add_argument("--patience", dest="patience", default=10, type=int)
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
train_idx, val_idx, test_idx = mask=='train', mask=='val', mask=='test'
x_train, x_val, x_test = X[train_idx], X[val_idx], X[test_idx]
y_train, y_val, y_test = y[train_idx], y[val_idx], y[test_idx]
del X,y



# define data loader
batch_size = args.batch_size
train_loader = DataLoader(list(zip(x_train, y_train)), batch_size=batch_size,
                            shuffle=True, num_workers=0, drop_last=False, pin_memory=True)
val_loader = DataLoader(list(zip(x_val, y_val)), batch_size=batch_size, 
                            shuffle=False, num_workers=0, drop_last=False, pin_memory=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

pretrained_weight=args.pretrained_weight
from MetaBrain import *
model = MetaBrain(num_target=n_tasks, load_base=False)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
load_params = torch.load("./best_model.pth.tar",map_location='cpu')['state_dict']
model_params = model.state_dict()

# extract the corresponding params
co_params = {k: v for k, v in load_params.items() if k in model_params.keys()}

model_params.update(co_params)
model.load_state_dict(model_params)

model.to(device)
os.makedirs("Log", exist_ok=True)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=time.strftime('./Log/log.%m%d.%H:%M:%S.txt'),
                    filemode='w')
logger = logging.getLogger(__name__)
logger.info('finish loading model')

learning_rate = args.lr
optimizer = optim.Adam(model.parameters(), lr=learning_rate)

warmup_epoches = 10
max_iter_epoches = 100
scheduler = CosineWarmupScheduler(optimizer=optimizer, warmup=warmup_epoches, max_iters=max_iter_epoches)

path2outprfx = './model_output/'
os.makedirs(path2outprfx,exist_ok=True)

# Training hyperparameters
num_epoches = 100
criterion = nn.BCELoss()
since = time.time()
logs = {"train_loss_list":[], 
            "val_loss_list":[], 
            "val_auc_cell":[], 
        "val_auc_peak":[]
            }
best_loss = 100.0
best_auc_cell = 0.1
best_auc_peak = 0.1

model.zero_grad() 
patience=args.patience
stale=0

########################
# train + valid
for epoch in range(num_epoches):
    logger.info('-' * 10)
    logger.info("Target Learning Rate: {}".format(str(learning_rate)))
    logger.info('Epoch {}/{}'.format(epoch + 1, num_epoches))

    
    epoch_target_list = []
    epoch_pred_list = []
    
    val_target_list = []
    val_pred_list = []
    
    running_loss = 0.0
    val_running_loss = 0.0
    
    model.train()
    for idx, data in enumerate(tqdm(train_loader,position=0)):
        
        
        # get the inputs; data is a list of [inputs, labels]
        inputs, targets = data
        inputs, targets= inputs.to(device), targets.to(device)

       
        
        outs = model(inputs)
        loss = criterion(outs, targets)
        
        loss = loss.nanmean()
        running_loss += loss.item()
        
        L2_regularization=1e-6
        for name,param in model.named_parameters():
            if 'classifier' in name:
                loss += L2_regularization*torch.norm(param)**2/2
        
        loss.backward() 
    

        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5) # Clip to avoid exploding gradient issues
        torch.cuda.empty_cache()   
        optimizer.step()
        optimizer.zero_grad()
        
    
  
        
    epoch_loss = running_loss / len(train_loader)
    logs['train_loss_list'].append(epoch_loss)
    
    model.eval()
    for idx, data in enumerate(tqdm(val_loader,position=0)):        
        inputs, targets = data
        inputs, targets= inputs.to(device), targets.to(device)
            
        with torch.no_grad():
            outputs= model.forward(inputs)
        val_loss = criterion(outputs, targets)
        
            
        val_loss = val_loss.nanmean()
        val_running_loss += val_loss.item()
        
        # append target and output
        val_target_list.append(targets.cpu().detach())
        val_pred_list.append(outputs.cpu().detach())
        
    # calculate auc statistics
    val_target_list = torch.cat(val_target_list, dim=0) 
    val_pred_list = torch.cat(val_pred_list, dim=0) 
    
    #auroc
    auc_peak_list = calculate_roc_peak(val_target_list, val_pred_list) 
   
    auc_peak_list=[auc_peak_list[k] for k in auc_peak_list.keys()] 
  
    val_auc_peak = float(np.nanmedian(auc_peak_list))  

    
    auc_cell_list = calculate_roc_cell(val_target_list, val_pred_list) 
  
    auc_cell_list=[auc_cell_list[k] for k in auc_cell_list.keys()] 
  
    val_auc_cell = float(np.nanmedian(auc_cell_list))  
   
    
    val_loss = val_running_loss / len(val_loader)
    logs['val_loss_list'].append(val_loss)
    logs['val_auc_peak'].append(val_auc_peak)
    logs['val_auc_cell'].append(val_auc_cell)
  
    
    show_train_log(loss_val=logs['val_loss_list'], output_fname='current_loss_val.pdf')
    show_train_log(loss_train=logs['train_loss_list'], output_fname='current_loss_train.pdf')
    show_train_log(acc_val=logs['val_auc_peak'], output_fname='current_auc_peak.pdf')
    show_train_log(acc_val=logs['val_auc_cell'], output_fname='current_auc_cell.pdf')
  
    show_train_log(loss_train=logs['train_loss_list'], 
                    loss_val=logs['val_loss_list'], 
                    output_fname='current_logs_loss.pdf')    


    scheduler.step() 
    curr_lr = scheduler.get_lr() 
    
    logger.info("Current Learning Rate: {}".format(str(curr_lr)))
    logger.info('Model_Loss: {:.4f}  Valid_Loss: {:.4f}  Valid_auc_peak: {:.4f}  Valid_auc_cell: {:.4f}' .format(
          epoch_loss, val_loss, val_auc_peak,val_auc_cell))
    torch.cuda.empty_cache()
    if(val_auc_cell >= best_auc_cell and  val_auc_peak >=best_auc_peak):
        best_auc_peak = val_auc_peak
        best_auc_cell = val_auc_cell
        best_loss = val_loss
        checkpoint = {'state_dict' : model.state_dict(), 
                      'current_model_loss': epoch_loss,
                      'current_valid_loss': val_loss,
                      'current_valid_auc_peak': val_auc_peak,
                      'current_valid_auc_cell': val_auc_cell
                     }
        save_checkpoint(checkpoint, path2outprfx +  'best_model.pth.tar')
        stale = 0
    else:
        stale += 1
        if stale > patience:
            with open('log.txt', 'a') as f:
                print(f"No improvment {patience} consecutive epochs, early stopping",file=f)
            break


    if(epoch + 1 == num_epoches):
        save_checkpoint(checkpoint, path2outprfx + '.'  + '.last_epoch.pth.tar')

time_elapsed = time.time() - since
logger.info('-' * 10)
logger.info('Training complete in {:.0f}m {:.0f}s'.format(time_elapsed // 60, time_elapsed % 60))

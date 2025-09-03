from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.metrics import auc, roc_curve, precision_recall_curve, average_precision_score
import matplotlib
import os
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
import torch.optim as optim
import numpy as np
import pandas as pd 
from tqdm import tqdm
def calculate_roc_cell(target, prediction):
    # assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    fpr, tpr, roc_auc = {}, {}, {} # orderedDict after python3.8
    n_classes = target.shape[-1]
    for index in tqdm(range(n_classes)):
        feature_targets = target[:, index]
        feature_preds = prediction[:, index]
        if len(np.unique(feature_targets)) > 1:
            fpr[index], tpr[index], _ = roc_curve(feature_targets, feature_preds)
            roc_auc[index] = auc(fpr[index], tpr[index])
        else:
          #  logging.warning("roc value was underestimated!")
            roc_auc[index] = 0
   
    return roc_auc

def calculate_roc_peak(target, prediction):
    # assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    fpr, tpr, roc_auc = {}, {}, {} # orderedDict after python3.8
    n_peak = target.shape[0]
    for index in tqdm(range(n_peak)):
        feature_targets = target[index]
        feature_preds = prediction[index]
        if len(np.unique(feature_targets)) > 1:
            fpr[index], tpr[index], _ = roc_curve(feature_targets, feature_preds)
            roc_auc[index] = auc(fpr[index], tpr[index])
        else:
          #  logging.warning("roc value was underestimated!")
            roc_auc[index] = 0
    
    return roc_auc


def save_checkpoint(state, filename):
    print("=> saving checkpoint")
    torch.save(state, filename)
    
def load_checkpoint(checkpoint):
    print("=> loading checkpoint")
    model.load_state_dict(checkpoint['state_dict'])
    


    
def show_train_log(loss_train=None, loss_val=None, acc_val=None, 
                   loss_test=None, acc_test=None,
                   lrs=None,
                    fig_size=(12,8),
                    save=True,
                    output_dir='Figures',
                    output_fname="Training_loss_log.pdf",
                    style="seaborn-colorblind",
                    fig_title="Training Log",
                    dpi=500):
    """function show train log.

    Parameters
    ----------
    loss_train : list
        traing loss
    loss_val : list
        validation loss
    kernel_size : int, optional
        Size of the convolving kernel

    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("PDF")
    
    plt.style.use(style)
    plt.figure()

    if loss_train:
        plt.plot(range(1, len(loss_train)+1), loss_train, 'b', label='Training Loss')
    if loss_val:
        plt.plot(range(1, len(loss_val)+1), loss_val, 'r', label='Validation Loss')
    if loss_test:
        plt.plot(range(1, len(loss_test)+1), loss_test, 'black', label='Test Loss')
    rate = 1
    if acc_val:
        plt.plot(range(1, len(acc_val)+1), list(map(lambda x: x*rate, acc_val)), 'g', label=str(rate)+'X Validation Accuracy')
    if acc_test:
        plt.plot(range(1, len(acc_test)+1), list(map(lambda x: x*rate, acc_test)), 'purple', label=str(rate)+'X Test Accuracy')
    if lrs:
        rate = int(1/(np.median(lrs)+1e-6))
        plt.plot(range(1, len(lrs)+1), list(map(lambda x: x*rate, lrs)), 'y--', alpha=0.2, label=str(rate)+'X Learning Rates')
    plt.title(fig_title)
    plt.legend()
    if save:
        plt.savefig(os.path.join(output_dir, output_fname),
                    format="pdf", dpi=dpi)
    else:
        plt.show()
    plt.close()
    

class CosineWarmupScheduler(optim.lr_scheduler._LRScheduler):

    def __init__(self, optimizer, warmup, max_iters):
        self.warmup = warmup
        self.max_num_iters = max_iters
        super().__init__(optimizer)

    def get_lr(self):
        lr_factor = self.get_lr_factor(epoch=self.last_epoch)
        return [base_lr * lr_factor for base_lr in self.base_lrs]

    def get_lr_factor(self, epoch):
        lr_factor = 0.5 * (1 + np.cos(np.pi * epoch / self.max_num_iters))
        if epoch <= self.warmup:
            lr_factor *= epoch * 1.0 / self.warmup
        return lr_factor

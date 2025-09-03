import h5py
import pysam
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix,csc_matrix
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import os
import anndata
import torch
sys.path.append('./arch/')
from arch import deepsea,MetaBrain,ResNet

# args
parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--model_weight")
parser.add_argument("--caqtl_file")
parser.add_argument("--genome_file")
parser.add_argument("--anno_file")
args = parser.parse_args()

snp_df=pd.read_csv(args.caqtl_file)
snp_df

peaks=snp_df.loc[:,['Chr','Start','End']]
peaks

onehot_nuc = {'A':[1,0,0,0],
            'C':[0,1,0,0],
            'G':[0,0,1,0],
            'T':[0,0,0,1],
            'N':[0,0,0,0]}
            
def _onehot_seq(seq):
    return np.array([onehot_nuc[nuc] for nuc in str(seq).upper()])

# genome
genome = pysam.Fastafile(args.genome_file)
genome

from tqdm import tqdm

seq_onehot = []

for peak in tqdm(peaks.values):
    seqnames, start, end = peak[:3]
    #chrom = str(seqnames.replace("chr",""))
    start, end = int(start), int(end)
    chrom = seqnames
    # catch overflowed error
    chrom_size = genome.get_reference_length(chrom)
    if end > chrom_size:
        print(peak[-1])
        pad = 'N' * (end - chrom_size) # pad N
        end = chrom_size
    # fetch sequence 
    seq = genome.fetch(reference=chrom, start=start, end=end)
    # pad N
    if start + 2000 > chrom_size:
        seq += pad
    # onehot    
    seq = _onehot_seq(seq)
    seq_onehot.append(seq)

seq_onehot = np.array(seq_onehot, dtype=np.bool)
genome.close()

seq_onehot.shape
seq_onehot_rc = seq_onehot[:, ::-1, ::-1]

#load model
h5file = h5py.File(args.data, 'r')
dim = h5file['pmat']['pmat_sc']['dim'][:]
y = np.zeros((dim[0], dim[1]), dtype = np.float32)
h5file.close()
n_tasks = y.shape[-1]
from MetaBrain import *
model = MetaBrain(num_target=n_tasks, load_base=False)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model_weight=args.model_weight
load_params = torch.load(model_weight,map_location='cpu')['state_dict']
model_params = model.state_dict()

# extract the corresponding params
co_params = {k: v for k, v in load_params.items() if k in model_params.keys()}

model_params.update(co_params)
model.load_state_dict(model_params)

model.eval()
model.to(device)

anno=pd.read_csv(args.anno_file)
anno

#dif
a= ['A','C','G','T']
b= [[ True, False, False, False],[ False,True, False, False],[ False, False,True,False],[False, False, False,True]]
dic = dict(zip(a, b))
all_ct_pred=pd.DataFrame()
for i in tqdm(range(seq_onehot.shape[0])):
    seq_ref_1hot = seq_onehot[i]
    seq_ref_1hot_ref=seq_ref_1hot.copy()
    seq_ref_1hot_ref[999]=dic[snp_df['Ref_Allele'][i]]
    output_ref = model(torch.from_numpy(seq_ref_1hot_ref.swapaxes(-1,0).astype(np.float32)).unsqueeze(0).to(device))
    seq_ref_1hot_alt=seq_ref_1hot.copy()
    seq_ref_1hot_alt[999]=dic[snp_df['Alt_Allele'][i]]
    output_alt = model(torch.from_numpy(seq_ref_1hot_alt.swapaxes(-1,0).astype(np.float32)).unsqueeze(0).to(device))
    ref_pred = output_ref.cpu().detach().numpy()
    alt_pred = output_alt.cpu().detach().numpy()
    x=ref_pred
    y=alt_pred
    dif=np.log2(y/(1-y))-np.log2(x/(1-x))
    pred=pd.DataFrame(dif.T)
        
    ct_pred=pd.DataFrame(columns=['ct','pred','base','ref','alt'])
    ct_pred['ct']=anno['cell type']
    ct_pred['pred']=pred[0]
    stat=ct_pred.groupby(['ct']).mean()
    stat=stat.reset_index()
    stat['base']=snp_df['SnpId'][i]
    stat['ref']=snp_df['Ref_Allele'][i]
    stat['alt']=snp_df['Alt_Allele'][i]
    all_ct_pred=pd.concat([all_ct_pred,stat])
        
all_ct_pred.to_csv('MBRS.csv')


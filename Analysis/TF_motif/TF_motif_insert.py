import sys
import h5py, os, argparse, logging, time
import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
import anndata
from Bio import SeqIO

from utils import *
sys.path.append('./arch/')
from arch import deepsea,MetaBrain,ResNet

# args
parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--model_weight")
parser.add_argument("--tf_file")
args = parser.parse_args()


#load model
h5file = h5py.File(args.data, 'r')
dim = h5file['pmat']['pmat_sc']['dim'][:]
y = np.zeros((dim[0], dim[1]), dtype = np.float32)
h5file.close()
n_tasks = y.shape[-1]


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


TF=pd.read_csv(args.tf_file)
TF

TF_id=list(TF['tf_id'])
TF_name=list(TF['tf_name'])


motif_fasta_folder = 'jaspar/Homo_sapiens_motif_fasta/'
tf_ad_list=[]
fasta_bg = "%s/shuffled_peaks.fasta" % motif_fasta_folder
pred_bg = pred_on_fasta(fasta_bg, model, batch_size=500)

for tf in tqdm(TF_id):
    fasta_motif = "%s/shuffled_peaks_motifs/%s.fasta" % (motif_fasta_folder, tf)
  
    pred_motif = pred_on_fasta(fasta_motif, model, batch_size=500)
 
    tf_score = pred_motif.mean(axis=0) - pred_bg.mean(axis=0)
    tf_score = (tf_score - tf_score.mean()) / tf_score.std()
    tf_ad=anndata.AnnData(np.expand_dims(tf_score,axis=0))
    tf_ad_list.append(tf_ad)
    
tf_ad_total=anndata.concat(tf_ad_list,join='outer')


mat=pd.DataFrame(tf_ad_total.X.T,index=ad.obs_names,columns=list(TF_name))
mat.to_csv('./mat_insertion_tf.csv')

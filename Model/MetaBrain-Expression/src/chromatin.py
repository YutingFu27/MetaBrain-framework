import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

from scipy.sparse import coo_matrix
from scipy.sparse import load_npz,save_npz
import h5py
import numpy as np
import pandas as pd
import sys
sys.path.append('./arch/')
from arch import deepsea,MetaBrain,ResNet
# args
parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--model_weight")
args = parser.parse_args()

#input
import h5py
h5file = h5py.File(args.data, 'r')
dim = h5file['pmat']['pmat_sc']['dim'][:]
y = np.zeros((dim[0], dim[1]), dtype = np.float32)
h5file.close()
n_tasks = y.shape[-1]

model_weight=args.model_weight
load_params = torch.load(model_weight,map_location='cpu')['state_dict']
model_params = model.state_dict()

model = MetaBrain(num_target=n_tasks, load_base=False)

# extract the corresponding params
co_params = {k: v for k, v in load_params.items() if k in model_params.keys()}

model_params.update(co_params)
model.load_state_dict(model_params)

model.eval()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

import math

def encodeSeqs(seqs, inputsize=2000):
    """Convert sequences to 0-1 encoding and truncate to the input size.
    The output concatenates the forward and reverse complement sequence
    encodings.
    Args:
        seqs: list of sequences (e.g. produced by fetchSeqs)
        inputsize: the number of basepairs to encode in the output
    Returns:
        numpy array of dimension: (2 x number of sequence) x 4 x inputsize
    2 x number of sequence because of the concatenation of forward and reverse
    complement sequences.
    """
    seqsnp = np.zeros((len(seqs), 4, inputsize), np.bool_)

    mydict = {'A': np.asarray([1, 0, 0, 0]), 'C': np.asarray([0, 1, 0, 0]),
            'G': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
            'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
            'a': np.asarray([1, 0, 0, 0]), 'c': np.asarray([0, 1, 0, 0]),
            'g': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
            'n': np.asarray([0, 0, 0, 0]), '-': np.asarray([0, 0, 0, 0])}

    n = 0
    for line in seqs:
        cline = line[int(math.floor(((len(line) - inputsize) / 2.0))):int(math.floor(len(line) - (len(line) - inputsize) / 2.0))]
        for i, c in enumerate(cline):
            seqsnp[n, :, i] = mydict[c]
        n = n + 1

    # get the complementary sequences
    #dataflip = seqsnp[:, ::-1, ::-1]
    #seqsnp = np.concatenate([seqsnp, dataflip], axis=0)
    return seqsnp.astype("float32")


import pysam
hg19_path = "/public/home/guogjgroup/ggj/FYT/AI/Ref/hg19.fa"
genome =pysam.Fastafile(hg19_path)

#now read in the gene file
import pandas as pd
gene_file = "geneanno.csv"
geneanno=pd.read_csv('./geneanno.csv')
chrom_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 
               'chrX', 'chrY']
df=geneanno
df['seqnames'] = pd.Categorical(df['seqnames'], categories=chrom_order, ordered=True)
df_sorted = df.sort_values('seqnames')
geneanno=df_sorted
geneanno.to_csv('geneanno_sort.csv')

windowsize = 2000
maxshift=20000
shifts = np.array(list(range(-20000,20000,200)))
print("shifts:\n",shifts)
assert len(shifts)==200
from tqdm import tqdm
for chr_id in chrom_order:
    print(chr_id)
    gene_df=geneanno[geneanno['seqnames']==chr_id]
    predictions_withrc = []
    for gene in tqdm(gene_df.values):
        gene_id,symbol,chrom,strand,TSS,CAGE_TSS,gene_type=gene[:7]
        tss=int(CAGE_TSS)
        strand=(1 if strand=="+" else -1)
        seqs_to_predict = []
        for shift in shifts:
            seq=genome.fetch(reference=chrom, start=tss + shift*strand -
                               int(0.5*windowsize), end=tss + shift*strand + int(0.5*windowsize))
            seqs_to_predict.append(seq)
        seqsnp = encodeSeqs(seqs_to_predict)

        model_input = torch.from_numpy(np.array(seqsnp)).unsqueeze(3)
        rc_model_input = torch.from_numpy(np.array(seqsnp[:,::-1,::-1])).unsqueeze(3)
        model_input = model_input.to(device)
        rc_model_input = rc_model_input.to(device)
        print(model_input.shape)
        prediction = model(model_input.squeeze()).cpu().detach().numpy().copy()
        rc_prediction = model(rc_model_input.squeeze()).detach().cpu().numpy().copy()
    #predictions_fwdonly.append(prediction)
        predictions_withrc.append(np.array(0.5*(prediction+rc_prediction)))

#predictions_fwdonly=np.array(predictions_fwdonly)
    predictions_withrc=np.array(predictions_withrc)

    pos_weight_shifts = shifts
    pos_weights = np.vstack([
        np.exp(-0.01*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.02*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.05*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.1*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.2*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.01*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.02*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.05*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.1*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.2*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0)])
#reconstructed_expecto_fwdonly = np.sum(pos_weights[None,:,:,None]*predictions_fwdonly[:,None,:,:],axis=2)
    reconstructed_expecto_withrc = np.sum(pos_weights[None,:,:,None]*predictions_withrc[:,None,:,:],axis=2)
#Xreducedall=0.5*(reconstructed_expecto_fwdonly+reconstructed_expecto_withrc)
    print(reconstructed_expecto_withrc.shape)
    np.save(chr_id+'_chromatin.npy',reconstructed_expecto_withrc)

genome.close()

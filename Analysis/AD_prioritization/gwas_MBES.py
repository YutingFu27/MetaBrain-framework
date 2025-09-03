import argparse
import math
import pysam
import torch
from torch import nn
import numpy as np
import pandas as pd
import h5py

sys.path.append('./arch/')
from arch import deepsea,MetaBrain,ResNet

from utils import *

# args
parser = argparse.ArgumentParser()
parser.add_argument("--data")
parser.add_argument("--model_weight")
parser.add_argument("--vcf_file")
parser.add_argument("--genome_file")
parser.add_argument("--modellist_file")
args = parser.parse_args()

genome = pysam.Fastafile(args.genome_file)


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

model=model.eval()
model.to(device)
print('Model loaded')

CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']


inputfile = args.vcf_file
maxshift = 800
inputsize = 2000
batchSize = 32
windowsize = inputsize + 100


vcf = pd.read_csv(inputfile, sep='\t', header=None, comment='#')

# standardize
vcf.iloc[:, 0] = 'chr' + vcf.iloc[:, 0].map(str).str.replace('chr', '')
vcf = vcf[vcf.iloc[:, 0].isin(CHRS)]
vcf

for shift in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200)):
    refseqs = []
    altseqs = []
    ref_matched_bools = []
    for i in range(vcf.shape[0]):
        refseq, altseq, ref_matched_bool = fetchSeqs(
            vcf[0][i], vcf[1][i], vcf[2][i], vcf[3][i], shift=shift, inputsize=inputsize)
        refseqs.append(refseq)
        altseqs.append(altseq)
        ref_matched_bools.append(ref_matched_bool)

    if shift == 0:
        # only need to be checked once
        print("Number of variants with reference allele matched with reference genome:")
        print(np.sum(ref_matched_bools))
        print("Number of input variants:")
        print(len(ref_matched_bools))

    ref_encoded = encodeSeqs(refseqs, inputsize=inputsize).astype(np.float32)
    alt_encoded = encodeSeqs(altseqs, inputsize=inputsize).astype(np.float32)

    ref_preds = []
    for i in range(int(1 + (ref_encoded.shape[0]-1) / batchSize)):
        input = torch.from_numpy(ref_encoded[int(i*batchSize):int((i+1)*batchSize),:,:])
        input = input.to(device)
        ref_preds.append(model.forward(input).cpu().detach().numpy().copy())
    ref_preds = np.vstack(ref_preds)

    alt_preds = []
    for i in range(int(1 + (alt_encoded.shape[0]-1) / batchSize)):
        input = torch.from_numpy(alt_encoded[int(i*batchSize):int((i+1)*batchSize),:,:])
        input = input.to(device)
        alt_preds.append(model.forward(input).cpu().detach().numpy().copy())
    alt_preds = np.vstack(alt_preds)

    diff = alt_preds - ref_preds
    f = h5py.File(inputfile + '.shift_' + str(shift) + '.diff.h5', 'w')
    f.create_dataset('pred', data=diff)
    f.close()
    



import h5py
from six.moves import reduce
import pickle

#load resources
modelList = pd.read_table(args.modellist_file,sep='\t',header=0)
models = []
for file in modelList['ModelName']:
        with open(file,'rb') as f:
            model=pickle.load(f)
        models.append(model)
        



#load input data
maxshift = int(800)
snpEffectFilePattern='predicted.vcf.shift_SHIFT.diff.h5'
nfeatures=n_tasks
batchSize=500


snpEffects = []
for shift in [str(n) for n in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200))]:
    h5f = h5py.File(snpEffectFilePattern.replace(
        'SHIFT', shift), 'r')['/pred']

    
    index_start = 0
    index_end = int(h5f.shape[0] / 2)

    snp_temp = (np.asarray(h5f[index_start:index_end,:])+ np.asarray(h5f[index_start+int(h5f.shape[0]/2):index_end+int(h5f.shape[0]/2),:]))/2.0
    snpEffects.append(snp_temp)


coor = pd.read_csv(args.vcf_file,sep='\t',header=None,comment='#')
coor = coor.iloc[index_start:index_end,:]
coor

n_snps=coor.shape[0]
snpExpEffects = compute_effects(snpEffects, \
                               n_snps,
                                models, maxshift=maxshift, nfeatures=nfeatures,
                                batchSize = batchSize)
#write output
snpExpEffects_df = coor
snpExpEffects_df=pd.concat([snpExpEffects_df.reset_index(),pd.DataFrame(snpExpEffects, columns = modelList.iloc[:,1])],axis=1,ignore_index =False)


snpExpEffects_df.to_csv('./MBES.csv')
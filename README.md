# MetaBrain

**MetaBrain** is a dual-module deep learning framework, trained on multi-omics human brain data at single-cell resolution, jointly predicts single-cell chromatin accessibility and gene expression solely from sequence. 

<img src="https://github.com/YutingFu27/MetaBrain-framework/blob/main/Image/MetaBrain_arch.jpg" width=100% height=100%>



## Requirements

- Python packages

```python
anndata==0.12.2
Bio==1.8.0
biopython==1.81
h5py==3.8.0
matplotlib==3.6.2
numpy==1.21.5
pandas==1.3.5
pysam==0.21.0
scanpy==1.11.4
scikit_learn==1.2.2
scipy==1.8.1
seaborn==0.13.2
six==1.16.0
torch==2.0.1
tqdm==4.65.0
```



- R packages

```R
library(ggplot2)
library(dplyr)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)
library(MotifDb)
library(motifbreakR)
library(Seurat)
library(SeuratDisk)
library(ArchR)
library(chromVAR)
library(Signac)
```



- Other tools

```shell
> plink
> ldsc
```



## Tutorial: Training MetaBrain on single-cell multi-omics datasets (human frontal cortex)  



#### MetaBrain-Regulatory

###### (1) Prepare scATAC-seq datasets

[0_dataset_pre.ipynb](https://github.com/YutingFu27/MetaBrain-framework/blob/main/Model/MetaBrain-Regulatory/example/0_dataset_pre.ipynb)   

Example data can be downloaded [here](https://bis.zju.edu.cn/BrainCis/download.html)   



###### (2) Train 

```shell
sbatch 1_run_model_train.sh
```



###### (3) Evaluate

[2_Evaluate.ipynb](https://github.com/YutingFu27/MetaBrain-framework/blob/main/Model/MetaBrain-Regulatory/example/2_evaluate.ipynb)   







#### MetaBrain-Expression

###### (1) generate chromatin feature 

```shell
sbatch 0_chromatin.sh
```



###### (2) Prepare scRNA-seq datasets

 [1_exp.ipynb](https://github.com/YutingFu27/MetaBrain-framework/blob/main/Model/MetaBrain-Expression/example/1_exp.ipynb)   

Example data can be downloaded [here](https://bis.zju.edu.cn/BrainCis/download.html)   



###### (3) Regression

[2_ridge.ipynb](https://github.com/YutingFu27/MetaBrain-framework/blob/main/Model/MetaBrain-Expression/example/2_ridge.ipynb)   

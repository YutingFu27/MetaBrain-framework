#pre statictis
conda activate sldp

#preprocess
preprocesspheno --sumstats-stem caqtl \
    --refpanel-name KG3.95 \
    --svd-stem data/svds_95percent/ \
    --print-snps data/1000G_hm3_noMHC.rsid \
    --ldscores-chr data/LDscore/LDscore. \
    --ld-blocks data/pickrell_ldblocks.hg19.eur.bed \
    --bfile-chr data/plink_files/1000G.EUR.QC.


###pre signed anno
preprocessannot --sannot-chr Caqtl/neuron/signedanno/celltype/  \
    --bfile-chr data/plink_files/1000G.EUR.QC. \
    --print-snps data/1000G_hm3_noMHC.rsid \
    --ld-blocks data/pickrell_ldblocks.hg19.eur.bed


###run sldp
sldp --pss-chr Caqtl/neuron/sumstats/caqtl.KG3.95/ \
    --sannot-chr Caqtl/neuron/signedanno/celltype/ \
    --background-sannot-chr data/maf5/ \
    --outfile-stem caqtl_celltype \
    --ld-blocks data/pickrell_ldblocks.hg19.eur.bed \
    --svd-stem data/svds_95percent/ \
    --bfile-reg-chr data/plink_files/1000G.EUR.QC.hm3_noMHC. \
    --seed 0
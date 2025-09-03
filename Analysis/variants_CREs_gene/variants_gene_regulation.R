setwd('V2G/')
library(data.table)
library(dplyr)
library(GenomicRanges)

use.ct=use.ct
caqtl=fread(paste0(use.ct,'_caqtl.csv'))
head(caqtl)
caqtl=caqtl[order(caqtl$quantile,decreasing = T),]
caqtl=caqtl[caqtl$quantile>=0.75,]

caqtl$label1=paste(caqtl$SnpId,caqtl$Chr,caqtl$Pos,caqtl$Ref,caqtl$Alt,sep = '_')
caqtl$label2=paste(caqtl$SnpId,caqtl$Chr,caqtl$Pos,caqtl$Alt,caqtl$Ref,sep = '_')
head(caqtl)
caqtl.bak=caqtl

#overlap with cell-type-eqtl
caqtl=caqtl.bak
load('significant_eqtl/ct8_eqtl.rda')
table(ct8_eqtl$ct)
ct8_eqtl$label=paste(ct8_eqtl$SNP,ct8_eqtl$chr,ct8_eqtl$pos,ct8_eqtl$other_allele,ct8_eqtl$effect_allele,sep = '_')
ct_eqtl=ct8_eqtl[ct8_eqtl$ct==use.ct,]
caqtl_ov_ct_eqtl=merge(caqtl,ct_eqtl,by.x='SnpId',by.y='SNP')
caqtl_ov_ct_eqtl$symbol=reshape2::colsplit(caqtl_ov_ct_eqtl$gene,'_',names = c('c1','c2'))$c1


caqtl$ov_snp=0
caqtl[caqtl$SnpId%in%unique(caqtl_ov_ct_eqtl$SnpId),]$ov_snp=1

caqtl$ov_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  
  if(snp%in%unique(caqtl_ov_ct$SnpId)){
    print(i)
    use.gene=caqtl_ov_ct_eqtl[caqtl_ov_ct_eqtl$SnpId==snp,]$symbol%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    caqtl$ov_gene[i]=use.gene
  }
}

caqtl_ov_ct_eqtl=caqtl
save(caqtl_ov_ct_eqtl,file = 'caqtl_ov_ct_eqtl.rda')  

#overlap with gtex-eqtl
caqtl=caqtl.bak
load('/media/ggj/Guo-4T-AI/FYT/MultiOmic/CodeTest/20240822_eqtl/GTEx_Analysis_v8_eQTL_updated/gtex_eqtl.rda')
head(gtex_eqtl)
gtex_eqtl$label=paste(gtex_eqtl$rs_id_dbSNP155_GRCh38p13,gtex_eqtl$chr,gtex_eqtl$variant_pos,gtex_eqtl$ref,gtex_eqtl$alt,sep = '_')
head(gtex_eqtl)
caqtl_ovgtxeqtl=merge(caqtl,gtex_eqtl,by.x='SnpId',by.y='rs_id_dbSNP155_GRCh38p13')
head(caqtl_ovgtxeqtl)



caqtl$ov_snp=0
caqtl[caqtl$SnpId%in%unique(caqtl_ovgtxeqtl$SnpId),]$ov_snp=1
table(caqtl$ov_snp)

caqtl$ov_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  
  if(snp%in%unique(caqtl_ovgtxeqtl$SnpId)){
    print(i)
    use.gene=caqtl_ovgtxeqtl[caqtl_ovgtxeqtl$SnpId==snp,]$gene_name%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    caqtl$ov_gene[i]=use.gene
  }
}

caqtl_ovgtxeqtl=caqtl
save(caqtl_ovgtxeqtl,file = 'caqtl_ovgtxqtl.rda')  

#LD
caqtl=caqtl.bak
load('/media/ggj/Guo-4T-AI/FYT/MultiOmic/CodeTest/20250124_brain_eqtl/LDblock/ct_ld.rda')

LD$symbol=reshape2::colsplit(LD$gene,'_',names = c('c1','c2'))$c1
caqtl_ld=merge(caqtl,LD,by.x='SnpId',by.y='SNP_A')

caqtl$ld_snp=''
caqtl$ld_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  if(snp%in%unique(caqtl_ld$SnpId)){
    print(i)
    use.snp=caqtl_ld[caqtl_ld$SnpId==snp,]$SNP_B%>%unique()
    use.gene=caqtl_ld[caqtl_ld$SnpId==snp,]$symbol%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    if (length(use.snp)>1) {
      use.snp=paste(use.snp,collapse = ',') 
    }
    caqtl$ld_snp[i]=use.snp
    caqtl$ld_gene[i]=use.gene
  }
}
caqtl_ovld=caqtl
save(caqtl_ovld,file = 'caqtl_ovld.rda')

#cicero
load('cicero.rda')
caqtl=caqtl.bak
caqtl$cicero_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  use.caqtl.gr=caqtl.gr[caqtl.gr$SNP==snp]
  ov=cicero.gr[queryHits(findOverlaps(cicero.gr,use.caqtl.gr))]%>%as.data.frame()
  ov
  if(dim(ov)[1]>=1){
    print(i)
    use.gene=ov$gene%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    caqtl$cicero_gene[i]=use.gene
    
  }
}

caqtl_cicero=caqtl
save(caqtl_cicero,file = './caqtl_cicero.rda')

#ABC
load('./ABC.rda')
caqtl=caqtl.bak
caqtl$ABC_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  use.caqtl.gr=caqtl.gr[caqtl.gr$SNP==snp]
  ov=ABC_result.gr[queryHits(findOverlaps(ABC_result.gr,use.caqtl.gr))]%>%as.data.frame()
  ov
  if(dim(ov)[1]>=1){
    print(i)
    use.gene=ov$gene%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    caqtl$ABC_gene[i]=use.gene
    
  }
}
caqtl_ABC=caqtl
save(caqtl_ABC,file = './OGC_caqtl_ABC.rda')

#Plac-seq
load('PLAC.rda')
caqtl=caqtl.bak
caqtl$PLAC_gene=''
for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  use.caqtl.gr=caqtl.gr[caqtl.gr$SNP==snp]
  ov=PLAC.gr[queryHits(findOverlaps(PLAC.gr,use.caqtl.gr))]%>%as.data.frame()
  ov
  if(dim(ov)[1]>=1){
    print(i)
    use.gene=ov$gene%>%unique()
    if (length(use.gene)>1) {
      use.gene=paste(use.gene,collapse = ',') 
    }
    caqtl$PLAC_gene[i]=use.gene
    
  }
}
caqtl_PLAC=caqtl
save(caqtl_PLAC,file = './caqtl_PLAC.rda')

#TFBS
caqtl=caqtl.bak
TFBS=fread('snp_motifbreadR_output.tsv')
TFBS
TFBS$variant=paste(TFBS$seqnames,TFBS$end,TFBS$REF,TFBS$ALT,sep = '_')
TFBS$variant[1:5]

caqtl$variant=paste(caqtl$Chr,caqtl$Pos,caqtl$Ref,caqtl$Alt,sep = '_')


intersect(caqtl$SnpId,TFBS$SNP_id)%>%length()


TFBS=TFBS[TFBS$variant%in%intersect(caqtl$variant,TFBS$variant),]


caqtl$TFBS_motif1=''
caqtl$TFBS_motif2=''

for (i in 1:dim(caqtl)[1]) {
  snp=caqtl$SnpId[i]
  
  if(snp%in%unique(TFBS$SNP_id)){
    print(i)
    use.TFBS=TFBS[TFBS$SNP_id==snp,]
    use.gene1=use.TFBS[use.TFBS$alleleDiff>0,]$geneSymbol%>%unique()
    use.gene2=use.TFBS[use.TFBS$alleleDiff<0,]$geneSymbol%>%unique()
    if (length(use.gene1)>1) {
      use.gene1=paste(use.gene1,collapse = ',') 
    }
    if (length(use.gene2)>1) {
      use.gene2=paste(use.gene2,collapse = ',') 
    }
    if (length(use.gene1)>=1) {
      caqtl$TFBS_motif1[i]=paste0("Gain:",use.gene1)
    }
    
    if (length(use.gene2)>=1) {
      caqtl$TFBS_motif2[i]=paste0("Loss:",use.gene2)
    }
    
  }
}



caqtl_TFBS=caqtl
save(caqtl_TFBS,file = './caqtl_TFBS.rda')

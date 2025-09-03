setwd('AD_Wightman_2021/')

library(data.table)
library(dplyr)
gwas=fread('Wightman_gwas.txt.gz')
head(gwas)
#filter significant gwas snps
gwas=gwas[as.numeric(gwas$P)<5e-8,]

gwas=gwas[gwas$A2%in%c('A','T','C','G','N'),]
gwas=gwas[gwas$A1%in%c('A','T','C','G','N'),]
fwrite(as.data.frame(gwas$SNP),file = 'snp.txt',quote = F,col.names = F,row.names = F,sep = '\t')

gwas$label1=paste(gwas$SNP,gwas$A1,gwas$A2,sep = '_')
gwas$label2=paste(gwas$SNP,gwas$A2,gwas$A1,sep = '_')

snp_anno=fread('./snp_anno.txt') #from ensembl vep
head(snp_anno)
snp_anno$label=paste(snp_anno$`#Uploaded_variation`,snp_anno$Allele,snp_anno$REF_ALLELE,sep='_')
snp_anno=snp_anno[snp_anno$`#Uploaded_variation`%in%gwas$SNP,]

table(snp_anno$Consequence)
snp_anno=snp_anno[!grepl('missense_variant',snp_anno$Consequence),]
snp_anno=snp_anno[!grepl('synonymous_variant',snp_anno$Consequence),]
snp_anno=snp_anno[!grepl('coding_sequence_variant',snp_anno$Consequence),]
unique(snp_anno$`#Uploaded_variation`)%>%length()

gwas=gwas[gwas$SNP%in%unique(snp_anno$`#Uploaded_variation`),]


#output for MetaBrain-Regulatory
head(gwas)
gwas_df=gwas[,c(8,1,2,4,3,6)]
head(gwas_df)
colnames(gwas_df)[c(1:5)]=c("SnpId","Chr","Pos","Ref_Allele","Alt_Allele")
snp=gwas_df
snp$Start=snp$Pos-1000
snp$End=snp$Pos+1000
snp$Chr=paste0('chr',snp$Chr)
head(snp)
fwrite(snp,file = 'gwas.csv')

#output for MetaBrain-Expression
aa<-snp[,2:5]
aa <- data.frame(aa)
head(aa)
write.table(aa, file = 'predicted.vcf',row.names = F,col.names = F,sep = "\t",quote = F)

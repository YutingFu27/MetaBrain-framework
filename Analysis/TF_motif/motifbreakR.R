library(ggplot2)
library(dplyr)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)
library(MotifDb)
library(motifbreakR)

setwd('motifbreakR/')
snp=read.csv('caqtl.csv')
head(snp)

variants <- snps.from.rsid(rsid = snp$SnpId,
                           dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                           search.genome = BSgenome.Hsapiens.UCSC.hg38)

results <- motifbreakR(snpList = variants, filterp = TRUE,
                       pwmList = subset(MotifDb, 
                                        dataSource %in% c("jaspar2022")),
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::SerialParam())

table(results$effect)
exportMBtable(results, file = "snp_motifbreakR_output.tsv", format = "tsv")
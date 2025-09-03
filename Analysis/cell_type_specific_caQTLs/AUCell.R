setwd('aucell/')

#atac mat
atac_mat=fread('2E_Accessibility_score.tsv.gz')
all.ct=colnames(atac_mat)
all.ct
dim(atac_mat)

atac_mat=atac_mat[,c('V1','GABA 1','GABA 2','Glutamatergic 1','Glutamatergic 2','Oligo Precursor','Oligodendrocyte','Microglia','Astrocyte 1','Astrocyte 2')]

peak=atac_mat$V1
peak=gsub(':','_',peak)
peak=gsub('-','_',peak)
peak
atac_mat=as.data.frame(atac_mat)
atac_mat=atac_mat[,-1]
rownames(atac_mat)=peak

colnames(atac_mat)=c('GABA1','GABA2','GLUT1','GLUT2','OPC','OGC','MGC','ASC1','ASC2')

chr=reshape2::colsplit(peak,'_',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'_',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'_',names = c('c1','c2','c3'))$c3
peak.gr=GRanges(seqnames = chr,IRanges(start,end),peak=peak)
peak.gr



cts=c('ASC','MGC','OGC','OPC','GABA','GLUT')
all_auc_df=data.frame()
for (use.ct in cts) {
  message(use.ct)
  snp_df=fread(paste0(use.ct,'_caqtl.csv'))
  snp_df=snp_df[snp_df$quantile>=0.75,] #top25% 
  snp.gr=GRanges(seqnames = snp_df$Chr,IRanges(snp_df$Pos,snp_df$Pos))
  ov.peak=peak.gr[queryHits(findOverlaps(peak.gr,snp.gr))]$peak

  
  ov.peak_sets=list(overlap=ov.peak)
  message('Start to perform AUCell')
  rank=AUCell_buildRankings(as.matrix(atac_mat),plotStats = T)
  auc_scores=AUCell_calcAUC(ov.peak_sets,rank)
  auc_df=as.data.frame(getAUC(auc_scores))
  auc_df=melt(auc_df)
  auc_df$group=use.ct
  all_auc_df=rbind(all_auc_df,auc_df)
}

auc_df_mat=dcast(all_auc_df,variable~group)


p=pheatmap::pheatmap(auc_df_mat,
                     color =rev(colorRampPalette(brewer.pal(n=7,name = 'RdBu'))(1000)),
                     scale = 'column',
                     cluster_rows = F,
                     cluster_cols = F,
                     angle_col = 315,
                     fontsize = 6,
                     width = 3,
                     height = 3,border_color = 'white',treeheight_row = 0,treeheight_col = 0)
ggsave(p,filename = './aucell.pdf',width = 6,height = 6)

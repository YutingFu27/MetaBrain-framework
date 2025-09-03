#Fisher
setwd('fisher/')

Q_score <- function(x){
  rank(x,ties.method = "min")/length(x)
}

library(RColorBrewer)
library(data.table)
library(ggplot2)
library(ArchR)
library(Seurat)
library(dplyr)
library(ggsci)
library(epitools)

all.ct=c('ASC','GABA','GLUT','MGC','OGC','OPC')
mycolor=pal_npg()(6)
names(mycolor)=all.ct


#for each ct
for (use.ct in all.ct) {
  message(use.ct)
  
  caqtl_use=fread(paste0(use.ct,'_caqtl.csv'))
  caqtl_re=fread(paste0(use.ct,'_caQTL_MBRS.csv'))
  caqtl_re=caqtl_re[caqtl_re$ct==use.ct,]
  caqtl_use$MBRS=abs(caqtl_use$pred)
  
  snp_df=caqtl_use
  snp_df$quantile=Q_score(abs(snp_df$MBRS))
  snp_df=snp_df[order(snp_df$quantile,decreasing = T),]
  
  peak_df=fread(paste0('ct_specific_peak/',use.ct,'_use.bed.gz'))
  peak.gr=GRanges(seqnames = peak_df$V1,IRanges(peak_df$V2,peak_df$V3))
  
  #fisher 
  message('Start to peform Fisher extact test')
  fisher_df=data.frame()
  for (cutoff in seq(25,95,5)) {
    message(cutoff)
    snp_df_lower=snp_df[snp_df$quantile<cutoff/100,]
    snp_df_higher=snp_df[snp_df$quantile>=cutoff/100,]
    
    snp_df_lower.gr=GRanges(seqnames = snp_df_lower$Chr,IRanges(snp_df_lower$Pos,snp_df_lower$Pos))
    snp_df_higher.gr=GRanges(seqnames = snp_df_higher$Chr,IRanges(snp_df_higher$Pos,snp_df_higher$Pos))
    
    queryHits(findOverlaps(snp_df_lower.gr,peak.gr))%>%unique()%>%length()/length(snp_df_lower.gr)
    queryHits(findOverlaps(snp_df_higher.gr,peak.gr))%>%unique()%>%length()/length(snp_df_higher.gr)
    
    a=queryHits(findOverlaps(snp_df_lower.gr,peak.gr))%>%unique()%>%length()
    b=length(snp_df_lower.gr)-a
    c=queryHits(findOverlaps(snp_df_higher.gr,peak.gr))%>%unique()%>%length()
    d=length(snp_df_higher.gr)-c
    data <- matrix(c(a,c,b,d),ncol = 2,dimnames = list(Group=c('lower','higher'),Result=c('Positive','Negative')))
    data=data[c(2:1),]
    fisher_result <- fisher.test(data)
    print(fisher_result$estimate)
    fishe_tmp=data.frame('OR'=as.numeric(fisher_result$estimate),'Lower_CI'=fisher_result$conf.int[1],'Upper_CI'=fisher_result$conf.int[2],'Quantile'=cutoff)
    fisher_df=rbind(fisher_df,fishe_tmp)
  }
  
  fisher_df$log2_OR=log2(fisher_df$OR)
  fisher_df$log2_lowCI=log2(fisher_df$Lower_CI)
  fisher_df$log2_upperCI=log2(fisher_df$Upper_CI)
  fisher_df$ct=use.ct
  
  message('Start to draw pictures')
  p=ggplot(fisher_df,aes(x=Quantile,y=log2_OR,fill=use.ct))+geom_ribbon(aes(ymin=log2_lowCI,ymax=log2_upperCI),alpha=0.3,color=NA)+geom_line(aes(color=ct),size=1)+scale_color_manual(values =as.character(mycolor[use.ct]))+scale_fill_manual(values = as.character(mycolor[use.ct]))+scale_x_continuous(breaks=seq(25,95,10))+theme_ArchR()+NoLegend()+ylab('log2OR(95%CI)')+xlab('MetaBrain score(%)')+geom_hline(yintercept = 0,linetype='dashed',color='darkred')

  ggsave(p,filename = paste0(use.ct,'_fisher.pdf'),height = 6,width = 6)
  message('Done')
}







setwd('japsar/')
library(data.table)
library(dplyr)
library(Biostrings)
library(rtracklayer)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(universalmotif)

source("utils.R")
dir.create("Homo_sapiens_motif_fasta")




bed <- import.bed("peaks.bed")
bed
set.seed(10)
examples <- sample(bed, 2000)
seqs <- get.seqs(BSgenome.Hsapiens.UCSC.hg38, examples, no.cores = 1)
seqs
writeXStringSet(seqs, "Homo_sapiens_motif_fasta/example_peaks.fasta", format="fasta", width=2000) 

cmd <- "fasta_ushuffle -k 2 < Homo_sapiens_motif_fasta/example_peaks.fasta > Homo_sapiens_motif_fasta/shuffled_peaks.fasta"
system(cmd)

# read motifs 
motifs <- read_meme("JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt")
shuffled_pks <- readDNAStringSet("Homo_sapiens_motif_fasta/shuffled_peaks.fasta", format="fasta")


# write fasta + motif 
dir.create("Homo_sapiens_motif_fasta/shuffled_peaks_motifs/")

all_tf_df=data.frame()
for (i in 1:length(motifs)) {
  tf <- motifs[[i]]@name
  tf_name<- motifs[[i]]@altname
  tf_df=data.frame('tf_id'=tf,'tf_name'=tf_name)
  all_tf_df=rbind(all_tf_df,tf_df)
  pwm <- motifs[[i]]@motif
  set.seed(10)
  out <- apply(pwm, 2, function(x) {
    return(sample(rownames(pwm), 2000, replace=T, prob=x))
  })
  motif_seqs <- apply(out, 1, function(x) paste(x, collapse=""))
  

  left_coord <- width(shuffled_pks)[1]/2 - floor(ncol(pwm)/2)
  left <- as.character(subseq(shuffled_pks, start=1, end=left_coord))
  right <- as.character(subseq(shuffled_pks, start=left_coord+ncol(pwm)+1, end=width(shuffled_pks)[1]))
  
  shuffled_pks_motifs <- DNAStringSet(paste0(left, motif_seqs, right))
  writeXStringSet(shuffled_pks_motifs, paste0("Homo_sapiens_motif_fasta/shuffled_peaks_motifs/", tf, ".fasta"), format="fasta")

  print(i)
  print(tf_name)
}
fwrite(all_tf_df,file = './tf_name_id.csv')

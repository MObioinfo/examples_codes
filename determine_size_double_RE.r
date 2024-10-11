
setwd("/net/fs-1/home01/maribaud/")
library(Biostrings)

setwd("/net/fs-1/home01/maribaud/salmon_salar_sequence")
read.table("styI_all.txt",header=TRUE)->styI
read.table("hindIII_all.txt",header=TRUE)->hindIII

rbind(styI,hindIII)->double_digest
double_digest[order(double_digest[,1]),]->double_digest2
double_digest2[order(double_digest2[,2]),]->double_digest3
 inter<-c()
 for (i in 1:length(double_digest3[,1]))
 {
 double_digest3[i+1,2]-double_digest3[i,2]+1->inter[i]
 }
 print("summary_size_frag_sty_hind")
 print(summary(inter))
 write.table(inter,"size_frag_sty_hind.txt",sep="\t",row.names=FALSE,col.names=FALSE)
 
 
  unique(c(as.character(styI[,4]),as.character(styI[,5])))->uniqseq
 unique(c(as.character(hindIII[,4]),as.character(hindIII[,5])))->uniqseq2
 
 nb_each_uniq_seq<-c()
 for (g in 1:length(uniqseq))
 {
 length(which(styI[,4]==uniqseq[g]))->nb_each_uniq_seq[g]
 }
 print("nb_uniq_seq_F_sty")
 print(length(subset(nb_each_uniq_seq,nb_each_uniq_seq==1)))
 
 nb_each_uniq_seq2<-c()
for (g in 1:length(uniqseq))
 {
 length(which(styI[,5]==uniqseq[g]))->nb_each_uniq_seq2[g]
 }
 print("nb_uniq_seq_R_sty")
 print(length(subset(nb_each_uniq_seq2,nb_each_uniq_seq2==1)))

 nb_each_uniq_seq3<-c()
 for (g in 1:length(uniqseq2))
 {
 length(which(hindIII[,4]==uniqseq2[g]))->nb_each_uniq_seq3[g]
 }
  print("nb_uniq_seq_F_hindIII")
 print(length(subset(nb_each_uniq_seq3,nb_each_uniq_seq3==1)))
 
 nb_each_uniq_seq4<-c()
for (g in 1:length(uniqseq2))
 {
 length(which(hindIII[,5]==uniqseq2[g]))->nb_each_uniq_seq4[g]
 }
  print("nb_uniq_seq_R_hindIII")
 print(length(subset(nb_each_uniq_seq4,nb_each_uniq_seq4==1)))
 
 
 
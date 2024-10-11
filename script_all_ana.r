setwd("C:/Users/maribaud/Desktop/RNA_seq_analysis_July2019/")
read.table("all_data_RNA_seq.txt",header=TRUE)->TPM
TPM[,-c(12:13)]->TPM
apply(TPM[,8:20],1,mean)->TPM$mean
apply(TPM[,8:20],1,sd)->TPM$sd
subset(TPM,TPM$sd!=0&TPM$mean!=0)->TPM2
(TPM2[,8]-TPM2$mean)/TPM2$sd->TPM2$eye_Zscore
(TPM2[,9]-TPM2$mean)/TPM2$sd->TPM2$head_kidney_Zscore
(TPM2[,10]-TPM2$mean)/TPM2$sd->TPM2$nose_Zscore
(TPM2[,11]-TPM2$mean)/TPM2$sd->TPM2$pyloric_caeca_Zscore
(TPM2[,12]-TPM2$mean)/TPM2$sd->TPM2$skin_Zscore
(TPM2[,13]-TPM2$mean)/TPM2$sd->TPM2$liver_Zscore
(TPM2[,14]-TPM2$mean)/TPM2$sd->TPM2$brain_Zscore
(TPM2[,15]-TPM2$mean)/TPM2$sd->TPM2$gill_Zscore
(TPM2[,16]-TPM2$mean)/TPM2$sd->TPM2$gut_Zscore
(TPM2[,17]-TPM2$mean)/TPM2$sd->TPM2$heart_Zscore
(TPM2[,18]-TPM2$mean)/TPM2$sd->TPM2$kidney_Zscore
(TPM2[,19]-TPM2$mean)/TPM2$sd->TPM2$muscle_Zscore
(TPM2[,20]-TPM2$mean)/TPM2$sd->TPM2$spleen_Zscore

write.table(TPM2,"data_genes_expression_without_sexual_tissues.txt", sep="\t",row.names=FALSE)
setwd("C:/Users/maribaud/Desktop/RNA_seq_analysis_July2019/")
read.table("data_genes_expression_without_sexual_tissues.txt", header=TRUE)->TPM2

library(fpc)
library(cluster)

data.matrix(TPM2[,23:35])->x
km_level10 <- kmeans(t(scale(t(x))), 10)
km_level11 <- kmeans(t(scale(t(x))), 11)
km_level12 <- kmeans(t(scale(t(x))), 12)
km_level13 <- kmeans(t(scale(t(x))), 13)
km_level14 <- kmeans(t(scale(t(x))), 14)

km_level10$cluster->TPM2$cluster_level10
km_level11$cluster->TPM2$cluster_level11
km_level12$cluster->TPM2$cluster_level12
km_level13$cluster->TPM2$cluster_level13
km_level14$cluster->TPM2$cluster_level14

TPM2[,36:40]->allu
library(alluvial)

aggregate(allu,by=allu[,1:5],length)->C
C[,1:6]->C
colnames(C)<-c("level10","level11","level12","level13","level14","freq")

png("alluvial_plot.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", ifelse(C[,1]==2,"orange",ifelse(C[,1]==3,"pink",
ifelse(C[,1]==4,"red",ifelse(C[,1]==5,"purple",ifelse(C[,1]==6,"blue",ifelse(C[,1]==7,"darkblue",ifelse(C[,1]==8,"green",ifelse(C[,1]==9,"darkgreen","grey"))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()


png("alluvial1.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 2, "orange", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial3.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 3, "pink", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial4.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 4, "red", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial5.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 5, "purple", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial6.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 6, "blue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial7.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 7, "darkblue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial8.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 8, "green", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial9.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 9, "darkgreen", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial10.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 10, "grey", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

subset(C,C$freq>50)->C

png("alluvial_plot2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", ifelse(C[,1]==2,"orange",ifelse(C[,1]==3,"pink",
ifelse(C[,1]==4,"red",ifelse(C[,1]==5,"purple",ifelse(C[,1]==6,"blue",ifelse(C[,1]==7,"darkblue",ifelse(C[,1]==8,"green",ifelse(C[,1]==9,"darkgreen","grey"))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()


png("alluvial1.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial2.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 2, "orange", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial3.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 3, "pink", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial4.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 4, "red", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial5.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 5, "purple", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial6.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 6, "blue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial7.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 7, "darkblue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial8.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 8, "green", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial9.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 9, "darkgreen", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial10.2.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 10, "grey", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

TPM2$final<-"NA"
 TPM2[which(TPM2$cluster_level10==1&TPM2$cluster_level11==5&TPM2$cluster_level12==10&TPM2$cluster_level13==12&TPM2$cluster_level14==3),41]<-1
 TPM2[which(TPM2$cluster_level10==1&TPM2$cluster_level11==5&TPM2$cluster_level12==10&TPM2$cluster_level13==12&TPM2$cluster_level14==1),41]<-2
  TPM2[which(TPM2$cluster_level10==2&TPM2$cluster_level11==11&TPM2$cluster_level12==12&TPM2$cluster_level13==13&TPM2$cluster_level14==2),41]<-3
   TPM2[which(TPM2$cluster_level10==2&TPM2$cluster_level11==11&TPM2$cluster_level12==3&TPM2$cluster_level13==13&TPM2$cluster_level14==10),41]<-4
    TPM2[which(TPM2$cluster_level10==3&TPM2$cluster_level11==4&TPM2$cluster_level12==1&TPM2$cluster_level13==3&TPM2$cluster_level14==14),41]<-5
	 TPM2[which(TPM2$cluster_level10==4&TPM2$cluster_level11==10&TPM2$cluster_level12==7&TPM2$cluster_level13==6&TPM2$cluster_level14==5),41]<-6
	  TPM2[which(TPM2$cluster_level10==5&TPM2$cluster_level11==7&TPM2$cluster_level12==4&TPM2$cluster_level13==11&TPM2$cluster_level14==9),41]<-7
	   TPM2[which(TPM2$cluster_level10==6&TPM2$cluster_level11==1&TPM2$cluster_level12==1&TPM2$cluster_level13==2&TPM2$cluster_level14==6),41]<-8
	   TPM2[which(TPM2$cluster_level10==7&TPM2$cluster_level11==8&TPM2$cluster_level12==8&TPM2$cluster_level13==8&TPM2$cluster_level14==11),41]<-9
	   TPM2[which(TPM2$cluster_level10==8&TPM2$cluster_level11==9&TPM2$cluster_level12==6&TPM2$cluster_level13==7&TPM2$cluster_level14==4),41]<-10
		TPM2[which(TPM2$cluster_level10==8&TPM2$cluster_level11==6&TPM2$cluster_level12==11&TPM2$cluster_level13==9&TPM2$cluster_level14==12),41]<-11
		 TPM2[which(TPM2$cluster_level10==8&TPM2$cluster_level11==6&TPM2$cluster_level12==5&TPM2$cluster_level13==5&TPM2$cluster_level14==8),41]<-12
		  TPM2[which(TPM2$cluster_level10==9&TPM2$cluster_level11==3&TPM2$cluster_level12==2&TPM2$cluster_level13==1&TPM2$cluster_level14==13),41]<-13
		   TPM2[which(TPM2$cluster_level10==10&TPM2$cluster_level11==2&TPM2$cluster_level12==9&TPM2$cluster_level13==10&TPM2$cluster_level14==7),41]<-14
		   TPM2[which(TPM2$cluster_level10==10&TPM2$cluster_level11==2&TPM2$cluster_level12==9&TPM2$cluster_level13==10&TPM2$cluster_level14==1),41]<-15
			
## creation of a supplemental cluster null with all the genes excluded from the clusters defined just previously

 ifelse(TPM2$final=="NA",16,TPM2$final)->TPM2$final

 length(subset(TPM2$final,TPM2$final==1))
 
 
subset(TPM2,TPM2$final==1)->C1
subset(TPM2,TPM2$final==2)->C2
subset(TPM2,TPM2$final==3)->C3
subset(TPM2,TPM2$final==4)->C4
subset(TPM2,TPM2$final==5)->C5
subset(TPM2,TPM2$final==6)->C6
subset(TPM2,TPM2$final==7)->C7
subset(TPM2,TPM2$final==8)->C8
subset(TPM2,TPM2$final==9)->C9
subset(TPM2,TPM2$final==10)->C10
subset(TPM2,TPM2$final==11)->C11
subset(TPM2,TPM2$final==12)->C12
subset(TPM2,TPM2$final==13)->C13
subset(TPM2,TPM2$final==14)->C14
subset(TPM2,TPM2$final==15)->C15
subset(TPM2,TPM2$final==16)->C16

library(ggplot2)

C1[,23:35]->c1bis
C2[,23:35]->c2bis
C3[,23:35]->c3bis
C4[,23:35]->c4bis
C5[,23:35]->c5bis
C6[,23:35]->c6bis
C7[,23:35]->c7bis
C8[,23:35]->c8bis
C9[,23:35]->c9bis
C10[,23:35]->c10bis
C11[,23:35]->c11bis
C12[,23:35]->c12bis
C13[,23:35]->c13bis
C14[,23:35]->c14bis
C15[,23:35]->c15bis
C16[,23:35]->c16bis



cbind(c1bis[,1],"eye")->s1
cbind(c1bis[,2],"head_kidney")->s2
cbind(c1bis[,3],"nose")->s3
cbind(c1bis[,4],"pyloric_caeca")->s4
cbind(c1bis[,5],"skin")->s5
cbind(c1bis[,6],"liver")->s6
cbind(c1bis[,7],"brain")->s7
cbind(c1bis[,8],"gill")->s8
cbind(c1bis[,9],"gut")->s9
cbind(c1bis[,10],"heart")->s10
cbind(c1bis[,11],"kidney")->s11
cbind(c1bis[,12],"muscle")->s12
cbind(c1bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster1_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster1_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster1_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 1","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster1_k16.png",500,800)
heatmap(data.matrix(c1bis))
dev.off()
png("heatmap_cluster1_2_k16.png",500,500)
heatmap(data.matrix(c1bis),Colv=NA)
dev.off()
png("heatmap_cluster1_3_k16.png",500,800)
heatmap(data.matrix(C1[,8:19]))
dev.off()
png("heatmap_cluster1_4_k16.png",500,500)
heatmap(data.matrix(C1[,8:19]),Colv=NA)
dev.off()

## all ready for presentation

cbind(c2bis[,1],"eye")->s1
cbind(c2bis[,2],"head_kidney")->s2
cbind(c2bis[,3],"nose")->s3
cbind(c2bis[,4],"pyloric_caeca")->s4
cbind(c2bis[,5],"skin")->s5
cbind(c2bis[,6],"liver")->s6
cbind(c2bis[,7],"brain")->s7
cbind(c2bis[,8],"gill")->s8
cbind(c2bis[,9],"gut")->s9
cbind(c2bis[,10],"heart")->s10
cbind(c2bis[,11],"kidney")->s11
cbind(c2bis[,12],"muscle")->s12
cbind(c2bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster2_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster2_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster2_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 2","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster2_k16.png",500,800)
heatmap(data.matrix(c2bis))
dev.off()
png("heatmap_cluster2_2_k16.png",500,500)
heatmap(data.matrix(c2bis),Colv=NA)
dev.off()
png("heatmap_cluster2_3_k16.png",500,800)
heatmap(data.matrix(C2[,8:22]))
dev.off()
png("heatmap_cluster2_4_k16.png",500,500)
heatmap(data.matrix(C2[,8:22]),Colv=NA)
dev.off()


cbind(c3bis[,1],"eye")->s1
cbind(c3bis[,2],"head_kidney")->s2
cbind(c3bis[,3],"nose")->s3
cbind(c3bis[,4],"pyloric_caeca")->s4
cbind(c3bis[,5],"skin")->s5
cbind(c3bis[,6],"liver")->s6
cbind(c3bis[,7],"brain")->s7
cbind(c3bis[,8],"gill")->s8
cbind(c3bis[,9],"gut")->s9
cbind(c3bis[,10],"heart")->s10
cbind(c3bis[,11],"kidney")->s11
cbind(c3bis[,12],"muscle")->s12
cbind(c3bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster3_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster3_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster3_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 3","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster3_k16.png",500,800)
heatmap(data.matrix(c3bis))
dev.off()
png("heatmap_cluster3_2_k16.png",500,500)
heatmap(data.matrix(c3bis),Colv=NA)
dev.off()
png("heatmap_cluster3_3_k16.png",500,800)
heatmap(data.matrix(C3[,8:22]))
dev.off()
png("heatmap_cluster3_4_k16.png",500,500)
heatmap(data.matrix(C3[,8:22]),Colv=NA)
dev.off()


cbind(c4bis[,1],"eye")->s1
cbind(c4bis[,2],"head_kidney")->s2
cbind(c4bis[,3],"nose")->s3
cbind(c4bis[,4],"pyloric_caeca")->s4
cbind(c4bis[,5],"skin")->s5
cbind(c4bis[,6],"liver")->s6
cbind(c4bis[,7],"brain")->s7
cbind(c4bis[,8],"gill")->s8
cbind(c4bis[,9],"gut")->s9
cbind(c4bis[,10],"heart")->s10
cbind(c4bis[,11],"kidney")->s11
cbind(c4bis[,12],"muscle")->s12
cbind(c4bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster4_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster4_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster4_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 4","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster4_k16.png",500,800)
heatmap(data.matrix(c4bis))
dev.off()
png("heatmap_cluster4_2_k16.png",500,500)
heatmap(data.matrix(c4bis),Colv=NA)
dev.off()
png("heatmap_cluster4_3_k16.png",500,800)
heatmap(data.matrix(C4[,8:22]))
dev.off()
png("heatmap_cluster4_4_k16.png",500,500)
heatmap(data.matrix(C4[,8:22]),Colv=NA)
dev.off()


cbind(c5bis[,1],"eye")->s1
cbind(c5bis[,2],"head_kidney")->s2
cbind(c5bis[,3],"nose")->s3
cbind(c5bis[,4],"pyloric_caeca")->s4
cbind(c5bis[,5],"skin")->s5
cbind(c5bis[,6],"liver")->s6
cbind(c5bis[,7],"brain")->s7
cbind(c5bis[,8],"gill")->s8
cbind(c5bis[,9],"gut")->s9
cbind(c5bis[,10],"heart")->s10
cbind(c5bis[,11],"kidney")->s11
cbind(c5bis[,12],"muscle")->s12
cbind(c5bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster5_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster5_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster5_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 5","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster5_k16.png",500,800)
heatmap(data.matrix(c5bis))
dev.off()
png("heatmap_cluster5_2_k16.png",500,500)
heatmap(data.matrix(c5bis),Colv=NA)
dev.off()
png("heatmap_cluster5_3_k16.png",500,800)
heatmap(data.matrix(C5[,8:22]))
dev.off()
png("heatmap_cluster5_4_k16.png",500,500)
heatmap(data.matrix(C5[,8:22]),Colv=NA)
dev.off()


cbind(c6bis[,1],"eye")->s1
cbind(c6bis[,2],"head_kidney")->s2
cbind(c6bis[,3],"nose")->s3
cbind(c6bis[,4],"pyloric_caeca")->s4
cbind(c6bis[,5],"skin")->s5
cbind(c6bis[,6],"liver")->s6
cbind(c6bis[,7],"brain")->s7
cbind(c6bis[,8],"gill")->s8
cbind(c6bis[,9],"gut")->s9
cbind(c6bis[,10],"heart")->s10
cbind(c6bis[,11],"kidney")->s11
cbind(c6bis[,12],"muscle")->s12
cbind(c6bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster6_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster6_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster6_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 6","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster6_k16.png",500,800)
heatmap(data.matrix(c6bis))
dev.off()
png("heatmap_cluster6_2_k16.png",500,500)
heatmap(data.matrix(c6bis),Colv=NA)
dev.off()
png("heatmap_cluster6_3_k16.png",500,800)
heatmap(data.matrix(C6[,8:22]))
dev.off()
png("heatmap_cluster6_4_k16.png",500,500)
heatmap(data.matrix(C6[,8:22]),Colv=NA)
dev.off()



cbind(c7bis[,1],"eye")->s1
cbind(c7bis[,2],"head_kidney")->s2
cbind(c7bis[,3],"nose")->s3
cbind(c7bis[,4],"pyloric_caeca")->s4
cbind(c7bis[,5],"skin")->s5
cbind(c7bis[,6],"liver")->s6
cbind(c7bis[,7],"brain")->s7
cbind(c7bis[,8],"gill")->s8
cbind(c7bis[,9],"gut")->s9
cbind(c7bis[,10],"heart")->s10
cbind(c7bis[,11],"kidney")->s11
cbind(c7bis[,12],"muscle")->s12
cbind(c7bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster7_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster7_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster7_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 7","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster7_k16.png",500,800)
heatmap(data.matrix(c7bis))
dev.off()
png("heatmap_cluster7_2_k16.png",500,500)
heatmap(data.matrix(c7bis),Colv=NA)
dev.off()
png("heatmap_cluster7_3_k16.png",500,800)
heatmap(data.matrix(C7[,8:22]))
dev.off()
png("heatmap_cluster7_4_k16.png",500,500)
heatmap(data.matrix(C7[,8:22]),Colv=NA)
dev.off()


cbind(c8bis[,1],"eye")->s1
cbind(c8bis[,2],"head_kidney")->s2
cbind(c8bis[,3],"nose")->s3
cbind(c8bis[,4],"pyloric_caeca")->s4
cbind(c8bis[,5],"skin")->s5
cbind(c8bis[,6],"liver")->s6
cbind(c8bis[,7],"brain")->s7
cbind(c8bis[,8],"gill")->s8
cbind(c8bis[,9],"gut")->s9
cbind(c8bis[,10],"heart")->s10
cbind(c8bis[,11],"kidney")->s11
cbind(c8bis[,12],"muscle")->s12
cbind(c8bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster8_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster8_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster8_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 8","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster8_k16.png",500,800)
heatmap(data.matrix(c8bis))
dev.off()
png("heatmap_cluster8_2_k16.png",500,500)
heatmap(data.matrix(c8bis),Colv=NA)
dev.off()
png("heatmap_cluster8_3_k16.png",500,800)
heatmap(data.matrix(C8[,8:22]))
dev.off()
png("heatmap_cluster8_4_k16.png",500,500)
heatmap(data.matrix(C8[,8:22]),Colv=NA)
dev.off()

cbind(c9bis[,1],"eye")->s1
cbind(c9bis[,2],"head_kidney")->s2
cbind(c9bis[,3],"nose")->s3
cbind(c9bis[,4],"pyloric_caeca")->s4
cbind(c9bis[,5],"skin")->s5
cbind(c9bis[,6],"liver")->s6
cbind(c9bis[,7],"brain")->s7
cbind(c9bis[,8],"gill")->s8
cbind(c9bis[,9],"gut")->s9
cbind(c9bis[,10],"heart")->s10
cbind(c9bis[,11],"kidney")->s11
cbind(c9bis[,12],"muscle")->s12
cbind(c9bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster9_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster9_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster9_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 9","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster9_k16.png",500,800)
heatmap(data.matrix(c9bis))
dev.off()
png("heatmap_cluster9_2_k16.png",500,500)
heatmap(data.matrix(c9bis),Colv=NA)
dev.off()
png("heatmap_cluster9_3_k16.png",500,800)
heatmap(data.matrix(C9[,8:22]))
dev.off()
png("heatmap_cluster9_4_k16.png",500,500)
heatmap(data.matrix(C9[,8:22]),Colv=NA)
dev.off()


cbind(c10bis[,1],"eye")->s1
cbind(c10bis[,2],"head_kidney")->s2
cbind(c10bis[,3],"nose")->s3
cbind(c10bis[,4],"pyloric_caeca")->s4
cbind(c10bis[,5],"skin")->s5
cbind(c10bis[,6],"liver")->s6
cbind(c10bis[,7],"brain")->s7
cbind(c10bis[,8],"gill")->s8
cbind(c10bis[,9],"gut")->s9
cbind(c10bis[,10],"heart")->s10
cbind(c10bis[,11],"kidney")->s11
cbind(c10bis[,12],"muscle")->s12
cbind(c10bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster10_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster10_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster10_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 10","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster10_k16.png",500,800)
heatmap(data.matrix(c10bis))
dev.off()
png("heatmap_cluster10_2_k16.png",500,500)
heatmap(data.matrix(c10bis),Colv=NA)
dev.off()
png("heatmap_cluster10_3_k16.png",500,800)
heatmap(data.matrix(C10[,8:22]))
dev.off()
png("heatmap_cluster10_4_k16.png",500,500)
heatmap(data.matrix(C10[,8:22]),Colv=NA)
dev.off()



cbind(c11bis[,1],"eye")->s1
cbind(c11bis[,2],"head_kidney")->s2
cbind(c11bis[,3],"nose")->s3
cbind(c11bis[,4],"pyloric_caeca")->s4
cbind(c11bis[,5],"skin")->s5
cbind(c11bis[,6],"liver")->s6
cbind(c11bis[,7],"brain")->s7
cbind(c11bis[,8],"gill")->s8
cbind(c11bis[,9],"gut")->s9
cbind(c11bis[,10],"heart")->s10
cbind(c11bis[,11],"kidney")->s11
cbind(c11bis[,12],"muscle")->s12
cbind(c11bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster11_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster11_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster11_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 11","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster11_k16.png",500,800)
heatmap(data.matrix(c11bis))
dev.off()
png("heatmap_cluster11_2_k16.png",500,500)
heatmap(data.matrix(c11bis),Colv=NA)
dev.off()
png("heatmap_cluster11_3_k16.png",500,800)
heatmap(data.matrix(C11[,8:22]))
dev.off()
png("heatmap_cluster11_4_k16.png",500,500)
heatmap(data.matrix(C11[,8:22]),Colv=NA)
dev.off()


cbind(c12bis[,1],"eye")->s1
cbind(c12bis[,2],"head_kidney")->s2
cbind(c12bis[,3],"nose")->s3
cbind(c12bis[,4],"pyloric_caeca")->s4
cbind(c12bis[,5],"skin")->s5
cbind(c12bis[,6],"liver")->s6
cbind(c12bis[,7],"brain")->s7
cbind(c12bis[,8],"gill")->s8
cbind(c12bis[,9],"gut")->s9
cbind(c12bis[,10],"heart")->s10
cbind(c12bis[,11],"kidney")->s11
cbind(c12bis[,12],"muscle")->s12
cbind(c12bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster12_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster12_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster12_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 12","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster12_k16.png",500,800)
heatmap(data.matrix(c12bis))
dev.off()
png("heatmap_cluster12_2_k16.png",500,500)
heatmap(data.matrix(c12bis),Colv=NA)
dev.off()
png("heatmap_cluster12_3_k16.png",500,800)
heatmap(data.matrix(C12[,8:22]))
dev.off()
png("heatmap_cluster12_4_k16.png",500,500)
heatmap(data.matrix(C12[,8:22]),Colv=NA)
dev.off()


cbind(c13bis[,1],"eye")->s1
cbind(c13bis[,2],"head_kidney")->s2
cbind(c13bis[,3],"nose")->s3
cbind(c13bis[,4],"pyloric_caeca")->s4
cbind(c13bis[,5],"skin")->s5
cbind(c13bis[,6],"liver")->s6
cbind(c13bis[,7],"brain")->s7
cbind(c13bis[,8],"gill")->s8
cbind(c13bis[,9],"gut")->s9
cbind(c13bis[,10],"heart")->s10
cbind(c13bis[,11],"kidney")->s11
cbind(c13bis[,12],"muscle")->s12
cbind(c13bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster13_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster13_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster13_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 13","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster13_k16.png",500,800)
heatmap(data.matrix(c13bis))
dev.off()
png("heatmap_cluster13_2_k16.png",500,500)
heatmap(data.matrix(c13bis),Colv=NA)
dev.off()
png("heatmap_cluster13_3_k16.png",500,800)
heatmap(data.matrix(C13[,8:22]))
dev.off()
png("heatmap_cluster13_4_k16.png",500,500)
heatmap(data.matrix(C13[,8:22]),Colv=NA)
dev.off()


cbind(c14bis[,1],"eye")->s1
cbind(c14bis[,2],"head_kidney")->s2
cbind(c14bis[,3],"nose")->s3
cbind(c14bis[,4],"pyloric_caeca")->s4
cbind(c14bis[,5],"skin")->s5
cbind(c14bis[,6],"liver")->s6
cbind(c14bis[,7],"brain")->s7
cbind(c14bis[,8],"gill")->s8
cbind(c14bis[,9],"gut")->s9
cbind(c14bis[,10],"heart")->s10
cbind(c14bis[,11],"kidney")->s11
cbind(c14bis[,12],"muscle")->s12
cbind(c14bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster14_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster14_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster14_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 14","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster14_k16.png",500,800)
heatmap(data.matrix(c14bis))
dev.off()
png("heatmap_cluster14_2_k16.png",500,500)
heatmap(data.matrix(c14bis),Colv=NA)
dev.off()
png("heatmap_cluster14_3_k16.png",500,800)
heatmap(data.matrix(C14[,8:22]))
dev.off()
png("heatmap_cluster14_4_k16.png",500,500)
heatmap(data.matrix(C14[,8:22]),Colv=NA)
dev.off()


cbind(c15bis[,1],"eye")->s1
cbind(c15bis[,2],"head_kidney")->s2
cbind(c15bis[,3],"nose")->s3
cbind(c15bis[,4],"pyloric_caeca")->s4
cbind(c15bis[,5],"skin")->s5
cbind(c15bis[,6],"liver")->s6
cbind(c15bis[,7],"brain")->s7
cbind(c15bis[,8],"gill")->s8
cbind(c15bis[,9],"gut")->s9
cbind(c15bis[,10],"heart")->s10
cbind(c15bis[,11],"kidney")->s11
cbind(c15bis[,12],"muscle")->s12
cbind(c15bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster15_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster15_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster15_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 15","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster15_k16.png",500,800)
heatmap(data.matrix(c15bis))
dev.off()
png("heatmap_cluster15_2_k16.png",500,500)
heatmap(data.matrix(c15bis),Colv=NA)
dev.off()
png("heatmap_cluster15_3_k16.png",500,800)
heatmap(data.matrix(C15[,8:22]))
dev.off()
png("heatmap_cluster15_4_k16.png",500,500)
heatmap(data.matrix(C15[,8:22]),Colv=NA)
dev.off()



cbind(c16bis[,1],"eye")->s1
cbind(c16bis[,2],"head_kidney")->s2
cbind(c16bis[,3],"nose")->s3
cbind(c16bis[,4],"pyloric_caeca")->s4
cbind(c16bis[,5],"skin")->s5
cbind(c16bis[,6],"liver")->s6
cbind(c16bis[,7],"brain")->s7
cbind(c16bis[,8],"gill")->s8
cbind(c16bis[,9],"gut")->s9
cbind(c16bis[,10],"heart")->s10
cbind(c16bis[,11],"kidney")->s11
cbind(c16bis[,12],"muscle")->s12
cbind(c16bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster16_k16_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster16_k16_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster16_k16.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 16","\n","(16 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster16_k16.png",500,800)
heatmap(data.matrix(c16bis))
dev.off()
png("heatmap_cluster16_2_k16.png",500,500)
heatmap(data.matrix(c16bis),Colv=NA)
dev.off()
png("heatmap_cluster16_3_k16.png",500,800)
heatmap(data.matrix(C16[,8:22]))
dev.off()
png("heatmap_cluster16_4_k16.png",500,500)
heatmap(data.matrix(C16[,8:22]),Colv=NA)
dev.off()

write.table(C1,"cluster1_method1.txt",sep="\t",row.names=FALSE)
write.table(C2,"cluster2_method1.txt",sep="\t",row.names=FALSE)
write.table(C3,"cluster3_method1.txt",sep="\t",row.names=FALSE)
write.table(C4,"cluster4_method1.txt",sep="\t",row.names=FALSE)
write.table(C5,"cluster5_method1.txt",sep="\t",row.names=FALSE)
write.table(C6,"cluster6_method1.txt",sep="\t",row.names=FALSE)
write.table(C7,"cluster7_method1.txt",sep="\t",row.names=FALSE)
write.table(C8,"cluster8_method1.txt",sep="\t",row.names=FALSE)
write.table(C9,"cluster9_method1.txt",sep="\t",row.names=FALSE)
write.table(C10,"cluster10_method1.txt",sep="\t",row.names=FALSE)
write.table(C11,"cluster11_method1.txt",sep="\t",row.names=FALSE)
write.table(C12,"cluster12_method1.txt",sep="\t",row.names=FALSE)
write.table(C13,"cluster13_method1.txt",sep="\t",row.names=FALSE)
write.table(C14,"cluster14_method1.txt",sep="\t",row.names=FALSE)
write.table(C15,"cluster15_method1.txt",sep="\t",row.names=FALSE)
write.table(C16,"cluster16_method1.txt",sep="\t",row.names=FALSE)

## without brain

setwd("C:/Users/maribaud/Desktop/RNA_seq_analysis_July2019/")
read.table("all_data_RNA_seq.txt",header=TRUE)->TPM
TPM[,-c(12,13,16)]->TPM
apply(TPM[,8:19],1,mean)->TPM$mean
apply(TPM[,8:19],1,sd)->TPM$sd
subset(TPM,TPM$sd!=0&TPM$mean!=0)->TPM2
(TPM2[,8]-TPM2$mean)/TPM2$sd->TPM2$eye_Zscore
(TPM2[,9]-TPM2$mean)/TPM2$sd->TPM2$head_kidney_Zscore
(TPM2[,10]-TPM2$mean)/TPM2$sd->TPM2$nose_Zscore
(TPM2[,11]-TPM2$mean)/TPM2$sd->TPM2$pyloric_caeca_Zscore
(TPM2[,12]-TPM2$mean)/TPM2$sd->TPM2$skin_Zscore
(TPM2[,13]-TPM2$mean)/TPM2$sd->TPM2$liver_Zscore
(TPM2[,14]-TPM2$mean)/TPM2$sd->TPM2$gill_Zscore
(TPM2[,15]-TPM2$mean)/TPM2$sd->TPM2$gut_Zscore
(TPM2[,16]-TPM2$mean)/TPM2$sd->TPM2$heart_Zscore
(TPM2[,17]-TPM2$mean)/TPM2$sd->TPM2$kidney_Zscore
(TPM2[,18]-TPM2$mean)/TPM2$sd->TPM2$muscle_Zscore
(TPM2[,19]-TPM2$mean)/TPM2$sd->TPM2$spleen_Zscore

write.table(TPM2,"data_genes_expression_without_sexual_tissues_and_without_brain.txt", sep="\t",row.names=FALSE)
setwd("C:/Users/maribaud/Desktop/RNA_seq_analysis_July2019/")
read.table("data_genes_expression_without_sexual_tissues_and_without_brain.txt", header=TRUE)->TPM2

library(fpc)
library(cluster)

data.matrix(TPM2[,22:33])->x
km_level10 <- kmeans(t(scale(t(x))), 10)
km_level11 <- kmeans(t(scale(t(x))), 11)
km_level12 <- kmeans(t(scale(t(x))), 12)
km_level13 <- kmeans(t(scale(t(x))), 13)
km_level14 <- kmeans(t(scale(t(x))), 14)

km_level10$cluster->TPM2$cluster_level10
km_level11$cluster->TPM2$cluster_level11
km_level12$cluster->TPM2$cluster_level12
km_level13$cluster->TPM2$cluster_level13
km_level14$cluster->TPM2$cluster_level14

TPM2[,34:38]->allu
library(alluvial)

aggregate(allu,by=allu[,1:5],length)->C
C[,1:6]->C
colnames(C)<-c("level10","level11","level12","level13","level14","freq")

png("alluvial_plotwobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", ifelse(C[,1]==2,"orange",ifelse(C[,1]==3,"pink",
ifelse(C[,1]==4,"red",ifelse(C[,1]==5,"purple",ifelse(C[,1]==6,"blue",ifelse(C[,1]==7,"darkblue",ifelse(C[,1]==8,"green",ifelse(C[,1]==9,"darkgreen","grey"))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()


png("alluvial1wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 2, "orange", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial3wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 3, "pink", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial4wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 4, "red", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial5wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 5, "purple", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial6wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 6, "blue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial7wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 7, "darkblue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial8wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 8, "green", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial9wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 9, "darkgreen", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial10wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 10, "grey", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

subset(C,C$freq>50)->C

png("alluvial_plot2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", ifelse(C[,1]==2,"orange",ifelse(C[,1]==3,"pink",
ifelse(C[,1]==4,"red",ifelse(C[,1]==5,"purple",ifelse(C[,1]==6,"blue",ifelse(C[,1]==7,"darkblue",ifelse(C[,1]==8,"green",ifelse(C[,1]==9,"darkgreen","grey"))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()


png("alluvial1.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 1, "yellow", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial2.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 2, "orange", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial3.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 3, "pink", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial4.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 4, "red", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial5.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 5, "purple", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial6.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 6, "blue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial7.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 7, "darkblue", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial8.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 8, "green", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial9.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 9, "darkgreen", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

png("alluvial10.2wobrain.png",800,500)
alluvial(C[,1:5],freq=C[,6],
col = ifelse(C[,1] == 10, "grey", "white"),
                 cex = 0.7,
				 border=NA
)
dev.off()

TPM2$final<-"NA"
 TPM2[which(TPM2$cluster_level10==1&TPM2$cluster_level11==7&TPM2$cluster_level12==6&TPM2$cluster_level13==13&TPM2$cluster_level14==12),39]<-1
 TPM2[which(TPM2$cluster_level10==1&TPM2$cluster_level11==7&TPM2$cluster_level12==6&TPM2$cluster_level13==10&TPM2$cluster_level14==12),39]<-2
  TPM2[which(TPM2$cluster_level10==2&TPM2$cluster_level11==3&TPM2$cluster_level12==9&TPM2$cluster_level13==9&TPM2$cluster_level14==2),39]<-3
   TPM2[which(TPM2$cluster_level10==3&TPM2$cluster_level11==4&TPM2$cluster_level12==12&TPM2$cluster_level13==5&TPM2$cluster_level14==14),39]<-4
    TPM2[which(TPM2$cluster_level10==4&TPM2$cluster_level11==1&TPM2$cluster_level12==10&TPM2$cluster_level13==11&TPM2$cluster_level14==9),39]<-5
	 TPM2[which(TPM2$cluster_level10==4&TPM2$cluster_level11==1&TPM2$cluster_level12==8&TPM2$cluster_level13==11&TPM2$cluster_level14==10),39]<-6
	  TPM2[which(TPM2$cluster_level10==4&TPM2$cluster_level11==1&TPM2$cluster_level12==8&TPM2$cluster_level13==11&TPM2$cluster_level14==4),39]<-7
	   TPM2[which(TPM2$cluster_level10==5&TPM2$cluster_level11==6&TPM2$cluster_level12==11&TPM2$cluster_level13==1&TPM2$cluster_level14==8),39]<-8
	   TPM2[which(TPM2$cluster_level10==5&TPM2$cluster_level11==6&TPM2$cluster_level12==11&TPM2$cluster_level13==4&TPM2$cluster_level14==6),39]<-9
	   TPM2[which(TPM2$cluster_level10==6&TPM2$cluster_level11==8&TPM2$cluster_level12==2&TPM2$cluster_level13==6&TPM2$cluster_level14==13),39]<-10
		TPM2[which(TPM2$cluster_level10==6&TPM2$cluster_level11==5&TPM2$cluster_level12==5&TPM2$cluster_level13==2&TPM2$cluster_level14==1),39]<-11
		 TPM2[which(TPM2$cluster_level10==7&TPM2$cluster_level11==9&TPM2$cluster_level12==4&TPM2$cluster_level13==12&TPM2$cluster_level14==5),39]<-12
		  TPM2[which(TPM2$cluster_level10==7&TPM2$cluster_level11==9&TPM2$cluster_level12==7&TPM2$cluster_level13==12&TPM2$cluster_level14==7),39]<-13
		   TPM2[which(TPM2$cluster_level10==8&TPM2$cluster_level11==11&TPM2$cluster_level12==12&TPM2$cluster_level13==7&TPM2$cluster_level14==14),39]<-14
		   TPM2[which(TPM2$cluster_level10==8&TPM2$cluster_level11==11&TPM2$cluster_level12==9&TPM2$cluster_level13==7&TPM2$cluster_level14==2),39]<-15
		    TPM2[which(TPM2$cluster_level10==9&TPM2$cluster_level11==10&TPM2$cluster_level12==1&TPM2$cluster_level13==3&TPM2$cluster_level14==11),39]<-16
			 TPM2[which(TPM2$cluster_level10==10&TPM2$cluster_level11==2&TPM2$cluster_level12==3&TPM2$cluster_level13==8&TPM2$cluster_level14==3),39]<-17
			 		  TPM2[which(TPM2$cluster_level10==10&TPM2$cluster_level11==2&TPM2$cluster_level12==7&TPM2$cluster_level13==8&TPM2$cluster_level14==7),39]<-18



			
## creation of a supplemental cluster null with all the genes excluded from the clusters defined just previously

 ifelse(TPM2$final=="NA",19,TPM2$final)->TPM2$final

 length(subset(TPM2$final,TPM2$final==1))
 
 
subset(TPM2,TPM2$final==1)->C1
subset(TPM2,TPM2$final==2)->C2
subset(TPM2,TPM2$final==3)->C3
subset(TPM2,TPM2$final==4)->C4
subset(TPM2,TPM2$final==5)->C5
subset(TPM2,TPM2$final==6)->C6
subset(TPM2,TPM2$final==7)->C7
subset(TPM2,TPM2$final==8)->C8
subset(TPM2,TPM2$final==9)->C9
subset(TPM2,TPM2$final==10)->C10
subset(TPM2,TPM2$final==11)->C11
subset(TPM2,TPM2$final==12)->C12
subset(TPM2,TPM2$final==13)->C13
subset(TPM2,TPM2$final==14)->C14
subset(TPM2,TPM2$final==15)->C15
subset(TPM2,TPM2$final==16)->C16
subset(TPM2,TPM2$final==17)->C17
subset(TPM2,TPM2$final==18)->C18
subset(TPM2,TPM2$final==19)->C19


library(ggplot2)

C1[,22:33]->c1bis
C2[,22:33]->c2bis
C3[,22:33]->c3bis
C4[,22:33]->c4bis
C5[,22:33]->c5bis
C6[,22:33]->c6bis
C7[,22:33]->c7bis
C8[,22:33]->c8bis
C9[,22:33]->c9bis
C10[,22:33]->c10bis
C11[,22:33]->c11bis
C12[,22:33]->c12bis
C13[,22:33]->c13bis
C14[,22:33]->c14bis
C15[,22:33]->c15bis
C16[,22:33]->c16bis
C17[,22:33]->c17bis
C18[,22:33]->c18bis
C19[,22:33]->c19bis



cbind(c1bis[,1],"eye")->s1
cbind(c1bis[,2],"head_kidney")->s2
cbind(c1bis[,3],"nose")->s3
cbind(c1bis[,4],"pyloric_caeca")->s4
cbind(c1bis[,5],"skin")->s5
cbind(c1bis[,6],"liver")->s6
cbind(c1bis[,7],"gill")->s7
cbind(c1bis[,8],"gut")->s8
cbind(c1bis[,9],"heart")->s9
cbind(c1bis[,10],"kidney")->s10
cbind(c1bis[,11],"muscle")->s11
cbind(c1bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster1_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster1_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster1_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 1","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster1_k19wobrain.png",500,800)
heatmap(data.matrix(c1bis))
dev.off()
png("heatmap_cluster1_2_k19wobrain.png",500,500)
heatmap(data.matrix(c1bis),Colv=NA)
dev.off()
png("heatmap_cluster1_3_k19wobrain.png",500,800)
heatmap(data.matrix(C1[,8:19]))
dev.off()
png("heatmap_cluster1_4_k19wobrain.png",500,500)
heatmap(data.matrix(C1[,8:19]),Colv=NA)
dev.off()

## all ready for presentation

cbind(c2bis[,1],"eye")->s1
cbind(c2bis[,2],"head_kidney")->s2
cbind(c2bis[,3],"nose")->s3
cbind(c2bis[,4],"pyloric_caeca")->s4
cbind(c2bis[,5],"skin")->s5
cbind(c2bis[,6],"liver")->s6
cbind(c2bis[,7],"gill")->s7
cbind(c2bis[,8],"gut")->s8
cbind(c2bis[,9],"heart")->s9
cbind(c2bis[,10],"kidney")->s10
cbind(c2bis[,11],"muscle")->s11
cbind(c2bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster2_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster2_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster2_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 2","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster2_k19wobrain.png",500,800)
heatmap(data.matrix(c2bis))
dev.off()
png("heatmap_cluster2_2_k19wobrain.png",500,500)
heatmap(data.matrix(c2bis),Colv=NA)
dev.off()
png("heatmap_cluster2_3_k19wobrain.png",500,800)
heatmap(data.matrix(C2[,8:22]))
dev.off()
png("heatmap_cluster2_4_k19wobrain.png",500,500)
heatmap(data.matrix(C2[,8:22]),Colv=NA)
dev.off()


cbind(c3bis[,1],"eye")->s1
cbind(c3bis[,2],"head_kidney")->s2
cbind(c3bis[,3],"nose")->s3
cbind(c3bis[,4],"pyloric_caeca")->s4
cbind(c3bis[,5],"skin")->s5
cbind(c3bis[,6],"liver")->s6
cbind(c3bis[,7],"gill")->s7
cbind(c3bis[,8],"gut")->s8
cbind(c3bis[,9],"heart")->s9
cbind(c3bis[,10],"kidney")->s10
cbind(c3bis[,11],"muscle")->s11
cbind(c3bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster3_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster3_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster3_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 3","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster3_k19wobrain.png",500,800)
heatmap(data.matrix(c3bis))
dev.off()
png("heatmap_cluster3_2_k19wobrain.png",500,500)
heatmap(data.matrix(c3bis),Colv=NA)
dev.off()
png("heatmap_cluster3_3_k19wobrain.png",500,800)
heatmap(data.matrix(C3[,8:22]))
dev.off()
png("heatmap_cluster3_4_k19wobrain.png",500,500)
heatmap(data.matrix(C3[,8:22]),Colv=NA)
dev.off()


cbind(c41bis[,1],"eye")->s1
cbind(c4bis[,2],"head_kidney")->s2
cbind(c4bis[,3],"nose")->s3
cbind(c4bis[,4],"pyloric_caeca")->s4
cbind(c4bis[,5],"skin")->s5
cbind(c4bis[,6],"liver")->s6
cbind(c4bis[,7],"gill")->s7
cbind(c4bis[,8],"gut")->s8
cbind(c4bis[,9],"heart")->s9
cbind(c4bis[,10],"kidney")->s10
cbind(c4bis[,11],"muscle")->s11
cbind(c4bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster4_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster4_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster4_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 4","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster4_k19wobrain.png",500,800)
heatmap(data.matrix(c4bis))
dev.off()
png("heatmap_cluster4_2_k19wobrain.png",500,500)
heatmap(data.matrix(c4bis),Colv=NA)
dev.off()
png("heatmap_cluster4_3_k19wobrain.png",500,800)
heatmap(data.matrix(C4[,8:22]))
dev.off()
png("heatmap_cluster4_4_k19wobrain.png",500,500)
heatmap(data.matrix(C4[,8:22]),Colv=NA)
dev.off()


cbind(c51bis[,1],"eye")->s1
cbind(c5bis[,2],"head_kidney")->s2
cbind(c5bis[,3],"nose")->s3
cbind(c5bis[,4],"pyloric_caeca")->s4
cbind(c5bis[,5],"skin")->s5
cbind(c5bis[,6],"liver")->s6
cbind(c5bis[,7],"gill")->s7
cbind(c5bis[,8],"gut")->s8
cbind(c5bis[,9],"heart")->s9
cbind(c5bis[,10],"kidney")->s10
cbind(c5bis[,11],"muscle")->s11
cbind(c5bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster5_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster5_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster5_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 5","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster5_k19wobrain.png",500,800)
heatmap(data.matrix(c5bis))
dev.off()
png("heatmap_cluster5_2_k19wobrain.png",500,500)
heatmap(data.matrix(c5bis),Colv=NA)
dev.off()
png("heatmap_cluster5_3_k19wobrain.png",500,800)
heatmap(data.matrix(C5[,8:22]))
dev.off()
png("heatmap_cluster5_4_k19wobrain.png",500,500)
heatmap(data.matrix(C5[,8:22]),Colv=NA)
dev.off()


cbind(c6bis[,1],"eye")->s1
cbind(c6bis[,2],"head_kidney")->s2
cbind(c6bis[,3],"nose")->s3
cbind(c6bis[,4],"pyloric_caeca")->s4
cbind(c6bis[,5],"skin")->s5
cbind(c6bis[,6],"liver")->s6
cbind(c6bis[,7],"gill")->s7
cbind(c6bis[,8],"gut")->s8
cbind(c6bis[,9],"heart")->s9
cbind(c6bis[,10],"kidney")->s10
cbind(c6bis[,11],"muscle")->s11
cbind(c6bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster6_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster6_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster6_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 6","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster6_k19wobrain.png",500,800)
heatmap(data.matrix(c6bis))
dev.off()
png("heatmap_cluster6_2_k19wobrain.png",500,500)
heatmap(data.matrix(c6bis),Colv=NA)
dev.off()
png("heatmap_cluster6_3_k19wobrain.png",500,800)
heatmap(data.matrix(C6[,8:22]))
dev.off()
png("heatmap_cluster6_4_k19wobrain.png",500,500)
heatmap(data.matrix(C6[,8:22]),Colv=NA)
dev.off()



cbind(c7bis[,1],"eye")->s1
cbind(c7bis[,2],"head_kidney")->s2
cbind(c7bis[,3],"nose")->s3
cbind(c7bis[,4],"pyloric_caeca")->s4
cbind(c7bis[,5],"skin")->s5
cbind(c7bis[,6],"liver")->s6
cbind(c7bis[,7],"gill")->s7
cbind(c7bis[,8],"gut")->s8
cbind(c7bis[,9],"heart")->s9
cbind(c7bis[,10],"kidney")->s10
cbind(c7bis[,11],"muscle")->s11
cbind(c7bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster7_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster7_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster7_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 7","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster7_k19wobrain.png",500,800)
heatmap(data.matrix(c7bis))
dev.off()
png("heatmap_cluster7_2_k19wobrain.png",500,500)
heatmap(data.matrix(c7bis),Colv=NA)
dev.off()
png("heatmap_cluster7_3_k19wobrain.png",500,800)
heatmap(data.matrix(C7[,8:22]))
dev.off()
png("heatmap_cluster7_4_k19wobrain.png",500,500)
heatmap(data.matrix(C7[,8:22]),Colv=NA)
dev.off()


cbind(c8bis[,1],"eye")->s1
cbind(c8bis[,2],"head_kidney")->s2
cbind(c8bis[,3],"nose")->s3
cbind(c8bis[,4],"pyloric_caeca")->s4
cbind(c8bis[,5],"skin")->s5
cbind(c8bis[,6],"liver")->s6
cbind(c8bis[,7],"gill")->s7
cbind(c8bis[,8],"gut")->s8
cbind(c8bis[,9],"heart")->s9
cbind(c8bis[,10],"kidney")->s10
cbind(c8bis[,11],"muscle")->s11
cbind(c8bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster8_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster8_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster8_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 8","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster8_k19wobrain.png",500,800)
heatmap(data.matrix(c8bis))
dev.off()
png("heatmap_cluster8_2_k19wobrain.png",500,500)
heatmap(data.matrix(c8bis),Colv=NA)
dev.off()
png("heatmap_cluster8_3_k19wobrain.png",500,800)
heatmap(data.matrix(C8[,8:22]))
dev.off()
png("heatmap_cluster8_4_k19wobrain.png",500,500)
heatmap(data.matrix(C8[,8:22]),Colv=NA)
dev.off()

cbind(c9bis[,1],"eye")->s1
cbind(c9bis[,2],"head_kidney")->s2
cbind(c9bis[,3],"nose")->s3
cbind(c9bis[,4],"pyloric_caeca")->s4
cbind(c9bis[,5],"skin")->s5
cbind(c9bis[,6],"liver")->s6
cbind(c9bis[,7],"gill")->s7
cbind(c9bis[,8],"gut")->s8
cbind(c9bis[,9],"heart")->s9
cbind(c9bis[,10],"kidney")->s10
cbind(c9bis[,11],"muscle")->s11
cbind(c9bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster9_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster9_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster9_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 9","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster9_k19wobrain.png",500,800)
heatmap(data.matrix(c9bis))
dev.off()
png("heatmap_cluster9_2_k19wobrain.png",500,500)
heatmap(data.matrix(c9bis),Colv=NA)
dev.off()
png("heatmap_cluster9_3_k19wobrain.png",500,800)
heatmap(data.matrix(C9[,8:22]))
dev.off()
png("heatmap_cluster9_4_k19wobrain.png",500,500)
heatmap(data.matrix(C9[,8:22]),Colv=NA)
dev.off()


cbind(c10bis[,1],"eye")->s1
cbind(c10bis[,2],"head_kidney")->s2
cbind(c10bis[,3],"nose")->s3
cbind(c10bis[,4],"pyloric_caeca")->s4
cbind(c10bis[,5],"skin")->s5
cbind(c10bis[,6],"liver")->s6
cbind(c10bis[,7],"gill")->s7
cbind(c10bis[,8],"gut")->s8
cbind(c10bis[,9],"heart")->s9
cbind(c10bis[,10],"kidney")->s10
cbind(c10bis[,11],"muscle")->s11
cbind(c10bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster10_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster10_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster10_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 10","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster10_k19wobrain.png",500,800)
heatmap(data.matrix(c10bis))
dev.off()
png("heatmap_cluster10_2_k19wobrain.png",500,500)
heatmap(data.matrix(c10bis),Colv=NA)
dev.off()
png("heatmap_cluster10_3_k19wobrain.png",500,800)
heatmap(data.matrix(C10[,8:22]))
dev.off()
png("heatmap_cluster10_4_k19wobrain.png",500,500)
heatmap(data.matrix(C10[,8:22]),Colv=NA)
dev.off()



cbind(c11bis[,1],"eye")->s1
cbind(c11bis[,2],"head_kidney")->s2
cbind(c11bis[,3],"nose")->s3
cbind(c11bis[,4],"pyloric_caeca")->s4
cbind(c11bis[,5],"skin")->s5
cbind(c11bis[,6],"liver")->s6
cbind(c11bis[,7],"gill")->s7
cbind(c11bis[,8],"gut")->s8
cbind(c11bis[,9],"heart")->s9
cbind(c11bis[,10],"kidney")->s10
cbind(c11bis[,11],"muscle")->s11
cbind(c11bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster11_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster11_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster11_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 11","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster11_k19wobrain.png",500,800)
heatmap(data.matrix(c11bis))
dev.off()
png("heatmap_cluster11_2_k19wobrain.png",500,500)
heatmap(data.matrix(c11bis),Colv=NA)
dev.off()
png("heatmap_cluster11_3_k19wobrain.png",500,800)
heatmap(data.matrix(C11[,8:22]))
dev.off()
png("heatmap_cluster11_4_k19wobrain.png",500,500)
heatmap(data.matrix(C11[,8:22]),Colv=NA)
dev.off()


cbind(c12bis[,1],"eye")->s1
cbind(c12bis[,2],"head_kidney")->s2
cbind(c12bis[,3],"nose")->s3
cbind(c12bis[,4],"pyloric_caeca")->s4
cbind(c12bis[,5],"skin")->s5
cbind(c12bis[,6],"liver")->s6
cbind(c12bis[,7],"gill")->s7
cbind(c12bis[,8],"gut")->s8
cbind(c12bis[,9],"heart")->s9
cbind(c12bis[,10],"kidney")->s10
cbind(c12bis[,11],"muscle")->s11
cbind(c12bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster12_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster12_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster12_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 12","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster12_k19wobrain.png",500,800)
heatmap(data.matrix(c12bis))
dev.off()
png("heatmap_cluster12_2_k19wobrain.png",500,500)
heatmap(data.matrix(c12bis),Colv=NA)
dev.off()
png("heatmap_cluster12_3_k19wobrain.png",500,800)
heatmap(data.matrix(C12[,8:22]))
dev.off()
png("heatmap_cluster12_4_k19wobrain.png",500,500)
heatmap(data.matrix(C12[,8:22]),Colv=NA)
dev.off()


cbind(c13bis[,1],"eye")->s1
cbind(c13bis[,2],"head_kidney")->s2
cbind(c13bis[,3],"nose")->s3
cbind(c13bis[,4],"pyloric_caeca")->s4
cbind(c13bis[,5],"skin")->s5
cbind(c13bis[,6],"liver")->s6
cbind(c13bis[,7],"gill")->s7
cbind(c13bis[,8],"gut")->s8
cbind(c13bis[,9],"heart")->s9
cbind(c13bis[,10],"kidney")->s10
cbind(c13bis[,11],"muscle")->s11
cbind(c13bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster13_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster13_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster13_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 13","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster13_k19wobrain.png",500,800)
heatmap(data.matrix(c13bis))
dev.off()
png("heatmap_cluster13_2_k19wobrain.png",500,500)
heatmap(data.matrix(c13bis),Colv=NA)
dev.off()
png("heatmap_cluster13_3_k19wobrain.png",500,800)
heatmap(data.matrix(C13[,8:22]))
dev.off()
png("heatmap_cluster13_4_k19wobrain.png",500,500)
heatmap(data.matrix(C13[,8:22]),Colv=NA)
dev.off()


cbind(c14bis[,1],"eye")->s1
cbind(c14bis[,2],"head_kidney")->s2
cbind(c14bis[,3],"nose")->s3
cbind(c14bis[,4],"pyloric_caeca")->s4
cbind(c14bis[,5],"skin")->s5
cbind(c14bis[,6],"liver")->s6
cbind(c14bis[,7],"gill")->s7
cbind(c14bis[,8],"gut")->s8
cbind(c14bis[,9],"heart")->s9
cbind(c14bis[,10],"kidney")->s10
cbind(c14bis[,11],"muscle")->s11
cbind(c14bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster14_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster14_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster14_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 14","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster14_k19wobrain.png",500,800)
heatmap(data.matrix(c14bis))
dev.off()
png("heatmap_cluster14_2_k19wobrain.png",500,500)
heatmap(data.matrix(c14bis),Colv=NA)
dev.off()
png("heatmap_cluster14_3_k19wobrain.png",500,800)
heatmap(data.matrix(C14[,8:22]))
dev.off()
png("heatmap_cluster14_4_k19wobrain.png",500,500)
heatmap(data.matrix(C14[,8:22]),Colv=NA)
dev.off()


cbind(c15bis[,1],"eye")->s1
cbind(c15bis[,2],"head_kidney")->s2
cbind(c15bis[,3],"nose")->s3
cbind(c15bis[,4],"pyloric_caeca")->s4
cbind(c15bis[,5],"skin")->s5
cbind(c15bis[,6],"liver")->s6
cbind(c15bis[,7],"gill")->s7
cbind(c15bis[,8],"gut")->s8
cbind(c15bis[,9],"heart")->s9
cbind(c15bis[,10],"kidney")->s10
cbind(c15bis[,11],"muscle")->s11
cbind(c15bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster15_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster15_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster15_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 15","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster15_k19wobrain.png",500,800)
heatmap(data.matrix(c15bis))
dev.off()
png("heatmap_cluster15_2_k19wobrain.png",500,500)
heatmap(data.matrix(c15bis),Colv=NA)
dev.off()
png("heatmap_cluster15_3_k19wobrain.png",500,800)
heatmap(data.matrix(C15[,8:22]))
dev.off()
png("heatmap_cluster15_4_k19wobrain.png",500,500)
heatmap(data.matrix(C15[,8:22]),Colv=NA)
dev.off()

cbind(c16bis[,1],"eye")->s1
cbind(c16bis[,2],"head_kidney")->s2
cbind(c16bis[,3],"nose")->s3
cbind(c16bis[,4],"pyloric_caeca")->s4
cbind(c16bis[,5],"skin")->s5
cbind(c16bis[,6],"liver")->s6
cbind(c16bis[,7],"gill")->s7
cbind(c16bis[,8],"gut")->s8
cbind(c16bis[,9],"heart")->s9
cbind(c16bis[,10],"kidney")->s10
cbind(c16bis[,11],"muscle")->s11
cbind(c16bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster16_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster16_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster16_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 16","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster16_k19wobrain.png",500,800)
heatmap(data.matrix(c16bis))
dev.off()
png("heatmap_cluster16_2_k19wobrain.png",500,500)
heatmap(data.matrix(c16bis),Colv=NA)
dev.off()
png("heatmap_cluster16_3_k19wobrain.png",500,800)
heatmap(data.matrix(C16[,8:22]))
dev.off()
png("heatmap_cluster16_4_k19wobrain.png",500,500)
heatmap(data.matrix(C16[,8:22]),Colv=NA)
dev.off()

cbind(c17bis[,1],"eye")->s1
cbind(c17bis[,2],"head_kidney")->s2
cbind(c17bis[,3],"nose")->s3
cbind(c17bis[,4],"pyloric_caeca")->s4
cbind(c17bis[,5],"skin")->s5
cbind(c17bis[,6],"liver")->s6
cbind(c17bis[,7],"gill")->s7
cbind(c17bis[,8],"gut")->s8
cbind(c17bis[,9],"heart")->s9
cbind(c17bis[,10],"kidney")->s10
cbind(c17bis[,11],"muscle")->s11
cbind(c17bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster17_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster17_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster17_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 17","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster17_k19wobrain.png",500,800)
heatmap(data.matrix(c17bis))
dev.off()
png("heatmap_cluster17_2_k19wobrain.png",500,500)
heatmap(data.matrix(c17bis),Colv=NA)
dev.off()
png("heatmap_cluster17_3_k19wobrain.png",500,800)
heatmap(data.matrix(C17[,8:22]))
dev.off()
png("heatmap_cluster17_4_k19wobrain.png",500,500)
heatmap(data.matrix(C17[,8:22]),Colv=NA)
dev.off()

cbind(c18bis[,1],"eye")->s1
cbind(c18bis[,2],"head_kidney")->s2
cbind(c18bis[,3],"nose")->s3
cbind(c18bis[,4],"pyloric_caeca")->s4
cbind(c18bis[,5],"skin")->s5
cbind(c18bis[,6],"liver")->s6
cbind(c18bis[,7],"gill")->s7
cbind(c18bis[,8],"gut")->s8
cbind(c18bis[,9],"heart")->s9
cbind(c18bis[,10],"kidney")->s10
cbind(c18bis[,11],"muscle")->s11
cbind(c18bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster18_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster18_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster18_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 18","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster18_k19wobrain.png",500,800)
heatmap(data.matrix(c18bis))
dev.off()
png("heatmap_cluster18_2_k19wobrain.png",500,500)
heatmap(data.matrix(c18bis),Colv=NA)
dev.off()
png("heatmap_cluster18_3_k19wobrain.png",500,800)
heatmap(data.matrix(C18[,8:22]))
dev.off()
png("heatmap_cluster18_4_k19wobrain.png",500,500)
heatmap(data.matrix(C18[,8:22]),Colv=NA)
dev.off()

cbind(c19bis[,1],"eye")->s1
cbind(c19bis[,2],"head_kidney")->s2
cbind(c19bis[,3],"nose")->s3
cbind(c19bis[,4],"pyloric_caeca")->s4
cbind(c19bis[,5],"skin")->s5
cbind(c19bis[,6],"liver")->s6
cbind(c19bis[,7],"gill")->s7
cbind(c19bis[,8],"gut")->s8
cbind(c19bis[,9],"heart")->s9
cbind(c19bis[,10],"kidney")->s10
cbind(c19bis[,11],"muscle")->s11
cbind(c19bis[,12],"spleen")->s12
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)->s
write.table(s,"cluster19_k19_genes_boxplot_wobrain.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster19_k19_genes_boxplot_wobrain.txt")->s
lev1=c("eye","skin","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster19_k19wobrain.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 19","\n","(19 clusters)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster19_k19wobrain.png",500,800)
heatmap(data.matrix(c19bis))
dev.off()
png("heatmap_cluster19_2_k19wobrain.png",500,500)
heatmap(data.matrix(c19bis),Colv=NA)
dev.off()
png("heatmap_cluster19_3_k19wobrain.png",500,800)
heatmap(data.matrix(C19[,8:22]))
dev.off()
png("heatmap_cluster19_4_k19wobrain.png",500,500)
heatmap(data.matrix(C19[,8:22]),Colv=NA)
dev.off()


write.table(C1,"cluster1_method2.txt",sep="\t",row.names=FALSE)
write.table(C2,"cluster2_method2.txt",sep="\t",row.names=FALSE)
write.table(C3,"cluster3_method2.txt",sep="\t",row.names=FALSE)
write.table(C4,"cluster4_method2.txt",sep="\t",row.names=FALSE)
write.table(C5,"cluster5_method2.txt",sep="\t",row.names=FALSE)
write.table(C6,"cluster6_method2.txt",sep="\t",row.names=FALSE)
write.table(C7,"cluster7_method2.txt",sep="\t",row.names=FALSE)
write.table(C8,"cluster8_method2.txt",sep="\t",row.names=FALSE)
write.table(C9,"cluster9_method2.txt",sep="\t",row.names=FALSE)
write.table(C10,"cluster10_method2.txt",sep="\t",row.names=FALSE)
write.table(C11,"cluster11_method2.txt",sep="\t",row.names=FALSE)
write.table(C12,"cluster12_method2.txt",sep="\t",row.names=FALSE)
write.table(C13,"cluster13_method2.txt",sep="\t",row.names=FALSE)
write.table(C14,"cluster14_method2.txt",sep="\t",row.names=FALSE)
write.table(C15,"cluster15_method2.txt",sep="\t",row.names=FALSE)
write.table(C16,"cluster16_method2.txt",sep="\t",row.names=FALSE)
write.table(C17,"cluster17_method2.txt",sep="\t",row.names=FALSE)
write.table(C18,"cluster18_method2.txt",sep="\t",row.names=FALSE)
write.table(C19,"cluster19_method2.txt",sep="\t",row.names=FALSE)


 for (i in 1:16)
 {
 read.table(paste("cluster",i,"_method1.txt",sep=""),header=TRUE)->A
 write.table(A,"cluster_method1.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
 }
read.table("cluster_method1.txt")->B
colnames(A)->colnames(B)
write.table(B,"cluster_method1.txt",sep="\t",row.names=FALSE)


 for (i in 1:19)
 {
 read.table(paste("cluster",i,"_method2.txt",sep=""),header=TRUE)->A
 write.table(A,"cluster_method2.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
 }
read.table("cluster_method2.txt")->B
colnames(A)->colnames(B)
write.table(B,"cluster_method2.txt",sep="\t",row.names=FALSE)

read.table("cluster_method1.txt",header=TRUE)->method1
read.table("cluster_method2.txt",header=TRUE)->method2
colnames(method1)[41]<-"cluster_method1"
colnames(method2)[39]<-"cluster_method2"

 total<-merge(method1,method2,by=c("gene_id","chrom","start","end","strand","gff_gene_id","ncbi_gene_id","eyeTPM","HKTPM","noseTPM","PCTPM","SkinTPM","liverTPM","gillTPM","gutTPM","heartTPM","kidneyTPM","muscleTPM","spleenTPM"))
## x = values method 1
## y = values method 2
## method1 keep more genes to analyse than method2
# we have to extract them to be sure to dont lost them; these genes are only expressed in the brain for sure ==> we are going to create an artifact cluster 20 in method2 to conserve them

 setdiff(method1[,1],total[,1])->method1_only
 posi<-c()
 
for (g in 1:length(method1_only))
 {
 which(method1_only[g]==method1[,1])->posi[g]
 }

 method1[posi,]->method1_only
 cbind(method1_only,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20)->compl 
 colnames(total)->colnames(compl)
 rbind(total,compl)->all_data
 
 ## now we can look at the difference of the distribution between the two methods
 
 all_data$cluster_method1->nb1
ifelse(nb1==5,1,ifelse(nb1==6,2,ifelse(nb1==8,3,ifelse(nb1==3,4,ifelse(nb1==4,5,ifelse(nb1==12,6,ifelse(nb1==13,7,ifelse(nb1==9,8,ifelse(nb1==11,9,ifelse(nb1==1,10,ifelse(nb1==2,11,ifelse(nb1==7,12,ifelse(nb1==10,13, ifelse(nb1==14,14, ifelse(nb1==16,15,ifelse(nb1==15,16,"NA"))))))))))))))))->all_data[,62]

 
 all_data$cluster_method2->nb2
ifelse(nb2==11,1,ifelse(nb2==12,2,ifelse(nb2==10,3,ifelse(nb2==2,4,ifelse(nb2==1,5,ifelse(nb2==5,6,ifelse(nb2==3,7,ifelse(nb2==15,8,ifelse(nb2==4,9,ifelse(nb2==14,10,ifelse(nb2==6,11,ifelse(nb2==7,12,ifelse(nb2==17,13, ifelse(nb2==18,14, ifelse(nb2==8,15,ifelse(nb2==16,16,ifelse(nb2==9,17,ifelse(nb2==13,18,19))))))))))))))))))->all_data[,63]

 
all_data[,62:63]->allu
library(alluvial)

aggregate(allu,by=allu[,1:2],length)->C
C[,1:3]->C
colnames(C)<-c("method1","method2","freq")
as.numeric(C$method1)->C$method1
as.numeric(C$method2)->C$method2


library(RColorBrewer)
###

png("alluvial_plotmethods2.png",800,500)
alluvial(C[,1:2],freq=C[,3],
col = ifelse(C[,1] == 1, "#BF812D", ifelse(C[,1]==2,"#543005",ifelse(C[,1]==3,"#9E0142",
ifelse(C[,1]==4,"#D53E4F",ifelse(C[,1]==5,"#F46D43",ifelse(C[,1]==6,"#FDAE61",ifelse(C[,1]==7,"#FEE08B",ifelse(C[,1]==8,"#FFFFBF",ifelse(C[,1]==9,"#E6F598",
ifelse(C[,1]==10,"#ABDDA4",ifelse(C[,1]==11,"#66C2A5",ifelse(C[,1]==12,"#3288BD",ifelse(C[,1]==13,"#5E4FA2",ifelse(C[,1]==14, "#9970AB",ifelse(C[,1]==15,"#8E0152","#C51B7D"))))))))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()

 
 subset(C,C[,3]>200)->C
 
png("alluvial_plotmethods3.png",800,500)
alluvial(C[,1:2],freq=C[,3],
col = ifelse(C[,1] == 1, "#BF812D", ifelse(C[,1]==2,"#543005",ifelse(C[,1]==3,"#9E0142",
ifelse(C[,1]==4,"#D53E4F",ifelse(C[,1]==5,"#F46D43",ifelse(C[,1]==6,"#FDAE61",ifelse(C[,1]==7,"#FEE08B",ifelse(C[,1]==8,"#FFFFBF",ifelse(C[,1]==9,"#E6F598",
ifelse(C[,1]==10,"#ABDDA4",ifelse(C[,1]==11,"#66C2A5",ifelse(C[,1]==12,"#3288BD",ifelse(C[,1]==13,"#5E4FA2",ifelse(C[,1]==14, "#9970AB",ifelse(C[,1]==15,"#8E0152","#C51B7D"))))))))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()

 C->C_save
subset(C,C[,3]>600)->C 

png("alluvial_plotmethods4.png",800,500)
alluvial(C[,1:2],freq=C[,3],
col = ifelse(C[,1] == 1, "#BF812D", ifelse(C[,1]==2,"#543005",ifelse(C[,1]==3,"#9E0142",
ifelse(C[,1]==4,"#D53E4F",ifelse(C[,1]==5,"#F46D43",ifelse(C[,1]==6,"#FDAE61",ifelse(C[,1]==7,"#FEE08B",ifelse(C[,1]==8,"#FFFFBF",ifelse(C[,1]==9,"#E6F598",
ifelse(C[,1]==10,"#ABDDA4",ifelse(C[,1]==11,"#66C2A5",ifelse(C[,1]==12,"#3288BD",ifelse(C[,1]==13,"#5E4FA2",ifelse(C[,1]==14, "#9970AB",ifelse(C[,1]==15,"#8E0152","#C51B7D"))))))))))))))),
                 cex = 0.7,
				 border=NA
)
dev.off()

write.table(C,"cluster_over_600_genes.txt",sep="\t",row.names=FALSE)
write.table(all_data, "all_data_RNA_seq.txt",sep="\t",row.names=FALSE)

ifelse(all_data[,62]==1&all_data[,63]==1,1,ifelse(all_data[,62]==2&all_data[,63]==2,2,ifelse(all_data[,62]==3&all_data[,63]==3,3,ifelse(all_data[,62]==4&all_data[,63]==4,4,ifelse(all_data[,62]==5&all_data[,63]==5,5,ifelse(all_data[,62]==6&all_data[,63]==6,6,ifelse(all_data[,62]==14&all_data[,63]==7,7,ifelse(all_data[,62]==7&all_data[,63]==7,8,ifelse(all_data[,62]==7&all_data[,63]==8,9,ifelse(all_data[,62]==14&all_data[,63]==9,10,ifelse(all_data[,62]==8&all_data[,63]==9,11,ifelse(all_data[,62]==8&all_data[,63]==10,12,ifelse(all_data[,62]==9&all_data[,63]==11,13,ifelse(all_data[,62]==10&all_data[,63]==13,14,ifelse(all_data[,62]==11&all_data[,63]==13,15,ifelse(all_data[,62]==15&all_data[,63]==13,16,ifelse(all_data[,62]==12&all_data[,63]==15,17,ifelse(all_data[,62]==14&all_data[,63]==15,18,ifelse(all_data[,62]==13&all_data[,63]==16,19,ifelse(all_data[,62]==15&all_data[,63]==16,20,ifelse(all_data[,62]==14&all_data[,63]==17,21,ifelse(all_data[,62]==10&all_data[,63]==19,22,ifelse(all_data[,62]==11&all_data[,63]==19,23,ifelse(all_data[,62]==14&all_data[,63]==19,24,ifelse(all_data[,62]==15&all_data[,63]==19,25,0)))))))))))))))))))))))))->all_data$final

for(i in 0:25)
 {
 (subset(all_data,all_data$final==i)[,7])->gene
 as.data.frame(gene)->gene
 write.table(gene,paste("gene_cluster",i,".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
 }

for (i in 0:25)
{
read.table(paste("gene_cluster",i,".txt",sep=""))->gene_list
cbind(paste("cluster",i,sep=""),gene_list[,1])->gene_list2
rbind(D,gene_list2)->D
}
write.table(D,"multi_list_cluster_ontology.txt",sep="\t",row.names=FALSE,col.names=FALSE)


subset(all_data,all_data$final==1)->C1
subset(all_data,all_data$final==2)->C2
subset(all_data,all_data$final==3)->C3
subset(all_data,all_data$final==4)->C4
subset(all_data,all_data$final==5)->C5
subset(all_data,all_data$final==6)->C6
subset(all_data,all_data$final==7)->C7
subset(all_data,all_data$final==8)->C8
subset(all_data,all_data$final==9)->C9
subset(all_data,all_data$final==10)->C10
subset(all_data,all_data$final==11)->C11
subset(all_data,all_data$final==12)->C12
subset(all_data,all_data$final==13)->C13
subset(all_data,all_data$final==14)->C14
subset(all_data,all_data$final==15)->C15
subset(all_data,all_data$final==16)->C16
subset(all_data,all_data$final==17)->C17
subset(all_data,all_data$final==18)->C18
subset(all_data,all_data$final==19)->C19
subset(all_data,all_data$final==20)->C20
subset(all_data,all_data$final==21)->C21
subset(all_data,all_data$final==22)->C22
subset(all_data,all_data$final==23)->C23
subset(all_data,all_data$final==24)->C24
subset(all_data,all_data$final==25)->C25
subset(all_data,all_data$final==0)->C0


library(ggplot2)
C0[,23:35]->c0bis
C1[,23:35]->c1bis
C2[,23:35]->c2bis
C3[,23:35]->c3bis
C4[,23:35]->c4bis
C5[,23:35]->c5bis
C6[,23:35]->c6bis
C7[,23:35]->c7bis
C8[,23:35]->c8bis
C9[,23:35]->c9bis
C10[,23:35]->c10bis
C11[,23:35]->c11bis
C12[,23:35]->c12bis
C13[,23:35]->c13bis
C14[,23:35]->c14bis
C15[,23:35]->c15bis
C16[,23:35]->c16bis
C17[,23:35]->c17bis
C18[,23:35]->c18bis
C19[,23:35]->c19bis
C20[,23:35]->c20bis
C21[,23:35]->c21bis
C22[,23:35]->c22bis
C23[,23:35]->c23bis
C24[,23:35]->c24bis
C25[,23:35]->c25bis

C0[,44:55]->c0bis2
C1[,44:55]->c1bis2
C2[,44:55]->c2bis2
C3[,44:55]->c3bis2
C4[,44:55]->c4bis2
C5[,44:55]->c5bis2
C6[,44:55]->c6bis2
C7[,44:55]->c7bis2
C8[,44:55]->c8bis2
C9[,44:55]->c9bis2
C10[,44:55]->c10bis2
C11[,44:55]->c11bis2
C12[,44:55]->c12bis2
C13[,44:55]->c13bis2
C14[,44:55]->c14bis2
C15[,44:55]->c15bis2
C16[,44:55]->c16bis2
C17[,44:55]->c17bis2
C18[,44:55]->c18bis2
C19[,44:55]->c19bis2
C20[,44:55]->c20bis2
C21[,44:55]->c21bis2
C22[,44:55]->c22bis2
C23[,44:55]->c23bis2
C24[,44:55]->c24bis2
C25[,44:55]->c25bis2


cbind(c1bis[,1],"eye")->s1
cbind(c1bis[,2],"head_kidney")->s2
cbind(c1bis[,3],"nose")->s3
cbind(c1bis[,4],"pyloric_caeca")->s4
cbind(c1bis[,5],"skin")->s5
cbind(c1bis[,6],"liver")->s6
cbind(c1bis[,7],"brain")->s7
cbind(c1bis[,8],"gill")->s8
cbind(c1bis[,9],"gut")->s9
cbind(c1bis[,10],"heart")->s10
cbind(c1bis[,11],"kidney")->s11
cbind(c1bis[,12],"muscle")->s12
cbind(c1bis[,13],"spleen")->s13
rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13)->s
write.table(s,"cluster1_Zscores_method1_genes_boxplot.txt",sep="\t",row.names=FALSE,col.names=FALSE)
read.table("cluster1_Zscores_method1_genes_boxplot.txt")->s
lev1=c("eye","skin","brain","nose","liver","pyloric_caeca","gut","gill","kidney","head_kidney","muscle","heart","spleen")
factor(s[,2],levels=lev1)->s[,2]
png("boxplot_cluster1_method1.png",800,500)
ggplot(data=s,aes(x=s[,2],y=s[,1]))+geom_boxplot(aes(fill=s[,2]))+ 
labs(title=paste("Z-score distribution for the genes in the cluster 1","\n","(method 1)",sep=""),y = "Z-score", fill="Tissue",x="Tissue")
dev.off()
png("heatmap_cluster1_k16.png",500,800)
heatmap(data.matrix(c1bis))
dev.off()
png("heatmap_cluster1_2_k16.png",500,500)
heatmap(data.matrix(c1bis),Colv=NA)
dev.off()
png("heatmap_cluster1_3_k16.png",500,800)
heatmap(data.matrix(C1[,8:19]))
dev.off()
png("heatmap_cluster1_4_k16.png",500,500)
heatmap(data.matrix(C1[,8:19]),Colv=NA)
dev.off()


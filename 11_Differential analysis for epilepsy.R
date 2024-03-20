---
title: "epilepsy"
author: "Qi Xiao"
date: "2022/9/27"
output: html_document
---

rm(list = ls())
#install.packages("pacman")
pacman::p_load(readxl,magrittr,
               RColorBrewer,
               Rtsne, umap,
               pheatmap,vioplot,
               ggpubr, ggplot2,ggrepel,
               corrplot, stringr,reshape2,preprocessCore)

setwd("D:/Guomics/HBA/Data_analysis20220901")

#-Read------------------------------------------------------------

df <- read.csv('HBA_prot_matrix20220811.csv')
info <- read.csv("FIleInfo_HBA_20220816_withoutFALSE_and_pool.csv",header = T,sep = ",",stringsAsFactors = F)
# info <- info[-c(which(info$Use=="FALSE"),which(info$pathology=="pool")),]

#-----------------------------------------------------------------
row.names(df) <- df[,1]
df <- df[,-1]
# df[is.na(df)] <- min(df,na.rm = T)/2
df <- t(df) %>% as.data.frame()
info2 <- info[info$sample%in%names(df),]

# df1<-df[,info2$sample]
info3<-info2
info3$for_epilepsy<-info3$Lobe 
info3$for_epilepsy[which(info3$Gyrus=="Hippocampus_Hipp")]<-"Hippocampus_Hipp"
info3$for_epilepsy[which(info3$Gyrus=="Amygdala_Amyg")]<-"Amygdala_Amyg"
e <- unique(info3$for_epilepsy[grep("epilepsy",info3$pathology)])
print(e)

info4 <- info3[info3$W_G=="G",]
info4 <- info4[info4$for_epilepsy%in%e,]

df1 <- df[,info4$sample]


#--FC------------------------------------------------------------
df1 <- as.data.frame(t(df1))
df1$label <- info4$for_epilepsy
e_region<-unique(df1$label)
print(e_region)

df_epilepsy<-df1[info4$sample[which(info4$pathology=="epilepsy")],]
df_HBA<-df1[info4$sample[which(info4$pathology=="health")],]

df_FC <- as.data.frame(matrix(nrow = 9639,ncol = 6),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])
for (i in 1:6) {
   data1<-df_HBA[which(df_HBA$label==e_region[i]),-ncol(df_HBA)]
   data2<-df_epilepsy[which(df_epilepsy$label==e_region[i]),-ncol(df_epilepsy)]
   data_mean1 <- apply(data1,2,function(x){mean(as.numeric(x))})
   data_mean2 <- apply(data2,2,function(x){mean(as.numeric(x))})
   df_FC[i] <- data_mean2-data_mean1
}
names(df_FC)[1:6] <- paste("FC",e_region,sep = "_")

#---wilcox-test-------------------------------------------------------

df_p_w <- as.data.frame(matrix(nrow = 9639,ncol = 6),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])
df_p_t <- as.data.frame(matrix(nrow = 9639,ncol = 6),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])
for (i in 1:6) {
   data1<-df_HBA[which(df_HBA$label==e_region[i]),-ncol(df_HBA)]
   data2<-df_epilepsy[which(df_epilepsy$label==e_region[i]),-ncol(df_epilepsy)]
   data_c <- rbind(data1,data2)
   df_p_w[i] <- apply(data_c,2,function(x){wilcox.test(x[1:nrow(data1)],x[(1+nrow(data1)):(nrow(data1)+nrow(data2))],exact = F)$p.value})
   df_p_t[i]<-apply(data_c,2,function(x){p<-try(t.test(x[1:nrow(data1)],x[(1+nrow(data1)):(nrow(data1)+nrow(data2))]))
   if(is(p,"try-error")) {NA} else {p$p.value}})
}

names(df_p_w)[1:6] <- paste("p_value_w",e_region,sep = "_")
names(df_p_t)[1:6] <- paste("p_value_t",e_region,sep = "_")

df_p <- cbind(df_p_w,df_p_t)

df_ap<- as.data.frame(apply(df_p,2,function(x){p.adjust(x,method = "BH")}))
names(df_ap_w)[1:6] <- paste("ap_value",e_region,sep = "_")

df_data <- cbind(df_FC,df_p,df_ap)
write.csv(df_data,"HBA_epilepsy_fc_p_6regions20220926.csv")

#----Select-proteins-2FC--BH-corrected-P<0.05----------------

#select proteins BH-P<0.05 & log2(FC)>2

# dif_w <- as.data.frame(matrix(nrow = 9639,ncol = 6),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])
# 
# for (i in 1:nrow(df_p_w)){
#   for(j in 1:ncol(df_p_w))
#     if(df_ap[i,j]<0.05&df_FC[i,j]>2){dif[i,j]<-1}   
# }
# 
# names(dif) <- names(df_ap)
# 
# select_pro <- dif[which(apply(dif,1,function(x){sum(!is.na(x))>0})),]  #log2(fc)>3 3771 proteins 
# write.csv(select_pro,"normal_epilepsy_pairwise_fc_2_select_5558proteins.csv")
# 
# diff_pro_num <- as.data.frame(apply(select_pro,2,function(x){sum(x,na.rm = T)}))
# # FrontalLobe 1023
# # ParietalLobe 1062
# # Hippocampus_Hipp 1284
# # OccipitalLobe 898
# # TemporalLobe 1184
# # Amygdala_Amyg 1137
# 
# 
# dif_p <- as.data.frame(matrix(nrow = 9639,ncol = 12),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])
# 
# for (i in 1:nrow(df_p)){
#   for(j in 7:ncol(df_p))
#     if(df_ap[i,j]<0.05&df_FC[i,j-6]>2){dif_p[i,j]<-1}   
# }
# names(dif_p)[7:12] <- names(df_ap)[7:12]
# dif_p <- dif_p[,-c(1:6)]
# 
# select_pro_t <- dif_p[which(apply(dif_p,1,function(x){sum(!is.na(x))>0})),]  #log2(fc)>3 3771 proteins 
# # write.csv(select_pro,"normal_epilepsy_pairwise_fc_2_select_5558proteins.csv")
# 
# diff_pro_num_t <- as.data.frame(apply(select_pro_t,2,function(x){sum(x,na.rm = T)}))
# 
# # 1830
# # p_value_t_FrontalLobe 830/813
# # p_value_t_ParietalLobe 870/852
# # p_value_t_Hippocampus_Hipp 1361/1255
# # p_value_t_OccipitalLobe 927/682
# # p_value_t_TemporalLobe 1243/1184
# # p_value_t_Amygdala_Amyg 1216/1049
# 
# length(intersect(c(row.names(select_pro)[which(select_pro$ap_value_Amygdala_Amyg==1)]),c(row.names(select_pro_t)[which(select_pro_t$p_value_t_Amygdala_Amyg==1)])))

# decide use wilcox test
#--Select-proteins-3FC--BH-corrected-P<0.05-------------------

dif_w_5fc <- as.data.frame(matrix(nrow = 9639,ncol = 6),row.names = names(df_epilepsy)[-ncol(df_epilepsy)])

for (i in 1:nrow(df_p_w)){
  for(j in 1:ncol(df_p_w))
   if(df_ap[i,j]<0.05& df_FC[i,j] > 2.32) {dif_w_5fc[i,j]<-1}else if(df_ap[i,j]<0.05& df_FC[i,j] < -2.32){dif_w_5fc[i,j]<-2}
     # if(df_ap[i,j]<0.05& abs(df_FC[i,j]) > 2.32) {dif_w_3fc[i,j]<-1}
}

names(dif_w_5fc) <- e_region

select_pro_w_5fc <- dif_w_5fc[which(apply(dif_w_5fc,1,function(x){sum(!is.na(x))>0})),]  #log2(fc)>3 3771 proteins 
write.csv(select_pro_w_5fc,"epilepsy_pairwise_fc_5_select_3210proteins20220926.csv")

select_pro_w_5fc$prot <- row.names(select_pro_w_5fc)
sum_w_5fc <- melt(select_pro_w_5fc,id=c("prot"))
sum_w_5fc <- sum_w_5fc[,-1]
sum_fc <- as.data.frame(table(sum_w_5fc))
# diff_pro_num_w_3fc <- as.data.frame(apply(select_pro_w_3fc,2,function(x){sum(!is.na(x),na.rm = T)}))

#--annotation----------------------------------------------------

anno <- read.csv("D:/share 114/HBA_xiaoqi/DataAnalysis20220712/Drug_taget/Database/all_human_drug_target.csv",header = T,sep = ",",stringsAsFactors = F)
anno <- tidyr::separate_rows(anno,'UniProt.Accessions',sep = '\\|')

data_anno <- select_pro_w_5fc
data_anno$prot <- sapply(row.names(data_anno),function(x){strsplit(x,"\\.")[[1]][1]})
data_anno$type <- anno$Type.1[match(data_anno$prot,anno$UniProt.Accessions)]
data_anno <- data_anno[!is.na(data_anno$type),-c(which(names(data_anno)=="prot"))]
data_anno <- melt(data_anno,id="type")
data_anno <- data_anno[!is.na(data_anno$value),]
data_anno <- as.data.frame(table(data_anno),stringsAsFactors = F)
data_anno$Freq[which(data_anno$value==2)] <- (data_anno$Freq[which(data_anno$value==2)])*-1
data_anno$text <- abs(data_anno$Freq)

# data_anno[is.na(data_anno)] <- 0

setwd("D:/Guomics/HBA/Data_analysis20220901/epilepsy")

for (i in unique(data_anno$variable)) {
   
   pdf(paste0(i,"_epilepsy_annotation_barplot.pdf"))
   data <- data_anno[which(data_anno$variable==i),-which(names(data_anno)=="variable")]
   
   p <- ggplot(data,aes(x=type,y=Freq)) +
   geom_bar(stat = 'identity',aes(fill=value),width = 0.8,color='white')+
   geom_text(aes(label=text),size=3 ,position = position_dodge(width = 0.8),
             vjust=-0.3)+
   theme(text = element_text(size=15),axis.text.x = element_text( 
            angle = 45,vjust = .5))+
   ggtitle(paste0(i,"_drug_target"))+
   scale_fill_manual(values=c("#E69F00","#999999"))
   print(p)
   dev.off()
}

#--------------------------------------------------------------
select_pro_w_5fc <- select_pro_w_5fc[,-ncol(select_pro_w_5fc)]
select_pro_w_5fc[select_pro_w_5fc==2] <- 1

for(i in e_region){
  region<-row.names(select_pro_w_5fc)[which(select_pro_w_5fc[grep(i,e_region)]==1)]
  assign(i,region)
}
region_pro<-list(FrontalLobe,ParietalLobe,Hippocampus_Hipp,OccipitalLobe,TemporalLobe,Amygdala_Amyg)
# inter<-intersect(TL_G,intersect(Hipp_G,intersect(Amyg_G,intersect(FL_G,intersect(PL_G,OL_G)))))
inter<-as.data.frame(Reduce(intersect,region_pro))  #617proteins
colnames(inter)<-"proteins"
# TL_unique<-setdiff(TL_G,Reduce(union,list(Hipp_G,Amyg_G,FL_G,PL_G,OL_G)))
write.csv(inter,"all_e_region_intersection_proteins20220926.csv")

#install.packages("UpSetR")
library(UpSetR) 
#data: col_region;row_proteins
# for_upset<-data.frame(t(select_pro))
for_upset <- select_pro_w_5fc
for_upset[is.na(for_upset)] <- 0

upset(for_upset,mb.ratio = c(0.55,0.45),nsets = 6,nintersects=100,number.angles = 30,point.size = 4,line.size = 1,mainbar.y.label = "intersect_number",sets.x.label = "pro_number",text.scale = c(2,2,2,2,2,2),order.by = c('freq', 'degree'), decreasing = c(TRUE, TRUE),queries = list(list(query = intersects, params = c("FrontalLobe","ParietalLobe","Hippocampus_Hipp","OccipitalLobe","TemporalLobe","Amygdala_Amyg"), color = 'darksalmon')))

#--summary--------------------------------------------------------

df_na <- read.csv('HBA_prot_matrix_11746pros_2824samples_202200909.csv')
row.names(df_na) <- df_na[,1]
df_na <- df_na[,-1]
info_e <- info4[which(info4$pathology=="epilepsy"),]
df_na <- df_na[,info_e$sample]

summary_e <- as.data.frame(matrix(nrow = length(unique(info_e$for_epilepsy)),ncol = 5),row.names = unique(info_e$for_epilepsy))
names(summary_e)=c("iden","N","n","UpRegulated","DownRegulated")

n_e <- as.data.frame(table(info_e$for_epilepsy))

for (i in unique(info_e$for_epilepsy)) {
   tmp <- df_na[,which(info_e$for_epilepsy==i)]
   summary_e[i,1] <- length(which(apply(tmp,1,function(x){sum(!is.na(x))>0})))
   summary_e[i,2] <- length(unique(info_e$Donor_ID[which(info_e$for_epilepsy==i)]))
   summary_e[i,3] <- length(which(info_e$for_epilepsy==i))
}

summary_e <- summary_e[e_region,]
summary_e$UpRegulated <- sum_fc$Freq[1:6]
summary_e$DownRegulated <- sum_fc$Freq[7:12]
summary_e$total <- apply(summary_e,1,function(x){sum(x[4],x[5])})

write.csv(summary_e,"epilepsy_summary20220926.csv")
sum_plot <- summary_e[,c(1,4,5)]
sum_plot$region <- row.names(sum_plot)
sum_plot <- melt(sum_plot,id="region",value.name = "num")

p <- ggplot(sum_plot,aes(x=region,y=num)) +
 
 geom_bar(stat = 'identity',position="dodge",aes(fill=variable),width = 0.8,color='white')+
 geom_text(aes(label=num),size=3 ,position = position_dodge(width = 0.8), 
             vjust=-0.3)+
 scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

print(p)

#--volcano-plot---------------------------------------------------

df_v_fc <- df_FC
df_v_ap <- df_ap[1:6]
df_vol <- data.frame(df_v_fc[6],df_v_ap[6])
df_vol$threshold = factor(ifelse(df_vol[2] < 0.05 & abs(df_vol[1]) >= 2.32, ifelse(df_vol[1]>= 2.32 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df_vol$Gene <- row.names(df_vol)

ggplot(df_vol,aes(x=FC_Amygdala_Amyg,y=-log10(p_value_w_Amygdala_Amyg),color=threshold))+
   geom_point()+
   scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
   geom_text_repel(
      data = df_vol[df_vol[2]<0.00001&abs(df_vol[1])>2.23,],
      aes(label = Gene),
      size = 3,
      segment.color = "black", show.legend = FALSE )+
   theme_bw()+
   theme(
      legend.title = element_blank()
   )+
   ylab('-log10 (p-adj)')+
   xlab('log2 (FoldChange)')+
   geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
   geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
rm(list = ls())
options(stringsAsFactors = F)
#install.packages("pacman")
pacman::p_load(readxl,magrittr,
               RColorBrewer,
               Rtsne, umap,
               pheatmap,vioplot,
               ggpubr, ggplot2,
               corrplot, stringr)

setwd("D:/Guomics/HBA/Data_analysis20220901/ts_score")
#-Read------------------------------------------------------------

df <- read.csv("D:/Guomics/HBA/Data_analysis20220901/HBA_prot_matrix20220811.csv",header = T,sep = ",",stringsAsFactors = F)
info <- read.csv("D:/Guomics/HBA/Data_analysis20220901/FIleInfo_HBA_20220816_withoutFALSE_and_pool.csv",header = T,sep = ",",stringsAsFactors = F)
# df_hk_all <- read.csv("C:/Users/A/Desktop/Ullen_class_xiaoqi20210804/New_matrix/Observed_atleast_80%_27regions_4596p.csv")

#Each region average ---------------------------------------------

row.names(df) <- df$X
df <- df[,-1]
df <- as.data.frame(t(df))

info2 <- info[which(info$pathology=="health"),]
df <- df[,info2$sample]
df_G <- df[,info2$sample[which(info2$W_G=="G")]]
df_W <- df[,info2$sample[which(info2$W_G=="W")]]

#TS_score_function------------------------------------------------

source('./AdaTiSS_fn.R')
AdaTiSS = function(X, tiss.abd=NULL) {
  
  pop.fit.mx = matrix(NA, nrow(X), 6)
  rownames(pop.fit.mx) = rownames(X)
  colnames(pop.fit.mx) = c("n.observed", "gam.sel", "mu0.hat", "sd0.hat" , "pi0.hat" ,"crt")
  id.ls.1 = rownames(X)[rowSums(!is.na(X)) >= 20]
  length(id.ls.1)
  for (i in 1:length(id.ls.1)) { 
    if(i %% 500 == 0) print(i)
    id = id.ls.1[i]
    x.0 = as.numeric(X[id, ]) ## as.numeric by yujing
    gam.limit = as.numeric(ifelse(sum(!is.na(x.0)) <= 100, 1, 3)) ## as.numeric by yujing
    result.x = adapt.gam.rob.fit.fn(x.0, gam.seq=seq(0,gam.limit,by=0.1), bin.num=round(length(x.0)/10))
    pop.fit.mx[id,] = result.x[["est.hat"]]
  }
  
  id.ls.2 = setdiff(rownames(X), id.ls.1)
  length(id.ls.2)
  pop.info.2 = apply(X[id.ls.2, ], 1, function(x) 
    c(sum(!is.na(x)), median(x, na.rm=TRUE), mad(x, na.rm=TRUE), 
      sum(abs(x-median(x, na.rm=TRUE)) <= 2*mad(x, na.rm=TRUE), na.rm=TRUE )/sum(!is.na(x)) ))
  pop.info.2 = t(pop.info.2)
  pop.fit.mx[id.ls.2, c("n.observed", "mu0.hat", "sd0.hat" , "pi0.hat")] = pop.info.2
  pop.fit.mx[, 'sd0.hat'] = pmax(pop.fit.mx[, 'sd0.hat'], 0.01)
  
  ada.s = (X - outer(pop.fit.mx[rownames(X), "mu0.hat"], rep(1, ncol(X))))/
    outer(pop.fit.mx[rownames(X), "sd0.hat"], rep(1, ncol(X)))
  ##null:
  if (!is.null(tiss.abd)) {
    ada.z = (tiss.abd - outer(pop.fit.mx[rownames(tiss.abd), "mu0.hat"], rep(1, ncol(tiss.abd))))/
      outer(pop.fit.mx[rownames(tiss.abd), "sd0.hat"], rep(1, ncol(tiss.abd)))
  } else {
    ada.z = NULL
  }
  
  return(list(ada.s = ada.s, ada.z=ada.z, pop.fit.mx=pop.fit.mx))
}

#TS_score_calculate-make-a-fresh-start------------------------------

# Tss_G <- AdaTiSS(df_G)
# df_tss_G <- Tss_G$ada.s
# write.csv(df_tss_G,"protein_matrix_G_afterTss20220915.csv")
# 
# Tss_W <- AdaTiSS(df_W)
# df_tss_W <- Tss_W$ada.s
# write.csv(df_tss_W,"protein_matrix_W_afterTss20220915.csv")
# 
# df_tss <- df_tss_G[which(apply(df_tss_G,1,function(x){sum(!is.na(x))>0})),]
# df_tss <- df_tss_W[which(apply(df_tss_W,1,function(x){sum(!is.na(x))>0})),]

#TS_score_calculate-middle-out------------------------------------

df_tss_W <- read.csv("D:/Guomics/HBA/Data_analysis20220901/ts_score/protein_matrix_G_afterTss20220915.csv",stringsAsFactors = F)
row.names(df_tss_W) <- df_tss_W[,1]
df_tss_W <- df_tss_W[,-1]
df_tss <- df_tss_W[which(apply(df_tss_W,1,function(x){sum(!is.na(x))>0})),]

info <- read.csv("D:/Guomics/HBA/Data_analysis20220901/FIleInfo_HBA_20220816_withoutFALSE_and_pool.csv",header = T,sep = ",",stringsAsFactors = F)
info2 <- info[which(info$pathology=="health"),]

#--Overlapping-with-mRNA-------------------------------------------

gene2prot <- read.csv("D:/Guomics/HBA/Data_analysis20220901/mRNA/gene2protein_reviewed_uniprot.csv")
gene2prot <- tidyr::separate_rows(gene2prot, 'Gene.Names', sep = ' ')
mRNA_df <- read.csv("D:/Guomics/HBA/Data_analysis20220901/mRNA/Allen_mRNA_sub_matrix_9162prot_2113sample20220920.csv")
prot <- unique(gene2prot$Entry[which(gene2prot$Gene.Names%in%mRNA_df$X)])

df_tss <- df_tss[which(sapply(row.names(df_tss),function(x){strsplit(x,"\\.")[[1]][1]})%in%prot),]
df_tss$label <- sapply(row.names(df_tss),function(x){strsplit(x,"\\.")[[1]][1]})

dup <- unique(df_tss$label[duplicated(df_tss$label)])
for (i in dup){
  # region <- info2$File_name[which(info2$Region2 == i & info2$W_G =="G")]
  data <- df_tss[which(df_tss$label==i),]
  data <- data[,-ncol(data)]
  data_mean <- apply(data,2,function(x){mean(as.numeric(x))})
  df_tss <- df_tss[-which(df_tss$label==i),]
  df_tss <- rbind(df_tss,c(data_mean,i))
}
name_r <- df_tss$label
df_tss <- df_tss[,-ncol(df_tss)]
df_tss <- apply(df_tss,2,as.numeric)%>%as.data.frame()
row.names(df_tss) <- name_r 
#-----------------------------------------------------------------

plot(density(as.matrix(df_tss)))

df <- t(df_tss) %>% as.data.frame()
info2 <- info[info$sample%in%row.names(df),]
df <- df[info2$sample,]
# ann_col <- data.frame(region=info2$region_test,row.names = row.names(df))

#--Lobe-analysis---------------------------------------------------

# df$label <- paste(info2$Lobe,info2$W_G,sep = "_")
# 
# for (i in unique(df$label)){
#   # region <- info2$File_name[which(info2$Region2 == i & info2$W_G =="G")]
#   data <- df[which(df$label==i),]
#   data <- data[,-ncol(data)]
#   data_mean <- apply(data,2,function(x){median(as.numeric(x))})
#   df <- df[-which(df$label==i),]
#   df <- rbind(df,c(data_mean,i))
# }
# 
# row.names(df) <- df$label
# tissue_df <- df[,-ncol(df)]
# tissue_df2 <- apply(tissue_df,2,as.numeric)%>%as.data.frame()
# row.names(tissue_df2) <- row.names(tissue_df)

#----------------------------------------------------------------


df$label <- paste(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions,info2$W_G,sep = "_")

for (i in unique(df$label)){
  # region <- info2$File_name[which(info2$Region2 == i & info2$W_G =="G")]
  data <- df[which(df$label==i),]
  data <- data[,-ncol(data)]
  data_mean <- apply(data,2,function(x){median(as.numeric(x))})
  df <- df[-which(df$label==i),]
  df <- rbind(df,c(data_mean,i))
}

row.names(df) <- df$label
tissue_df <- df[,-ncol(df)]
tissue_df2 <- apply(tissue_df,2,as.numeric)%>%as.data.frame()
row.names(tissue_df2) <- row.names(tissue_df)
# tissue_df <- log2(tissue_df2)
# write.csv(tissue_df2,"protein_gray_ts_matrix_121region_median20230328.csv")

# density_tissue_df <- data.frame(x = density(as.matrix(tissue_df))$x, y = density(as.matrix(tissue_df))$y)
# for(i in 2:(nrow(density_tissue_df) - 1)){
#   density_tissue_df$isExtreLarge[i] <- (density_tissue_df$y[i+1] > density_tissue_df$y[i]) & (density_tissue_df$y[i-1] > density_tissue_df$y[i])
# }
# density_tissue_df_extre <- density_tissue_df[which(density_tissue_df$isExtreLarge ), ]
# plot(density_tissue_df_extre$x, density_tissue_df_extre$y)
# 0.60
tissue_df_lt2 <- tissue_df2 %>% .[, apply(., 2, function(x){ any(x > 2,na.rm = T) })]
tissue_df_lt2 %>% apply(., 2, max) %>% min

tissue_df_lt2.5 <- tissue_df2 %>% .[, apply(., 2, function(x){ any(x > 2.5,na.rm = T) })]
tissue_df_lt2.5 %>% apply(., 2, max) %>% min

tissue_df_lt3 <- tissue_df2 %>% .[, apply(., 2, function(x){ any(x > 3,na.rm = T) })]
tissue_df_lt3 %>% apply(., 2, max) %>% min

tissue_df_lt4 <- tissue_df2 %>% .[, apply(., 2, function(x){ any(x > 4,na.rm = T) })]
tissue_df_lt4 %>% apply(., 2, max) %>% min

# x1 <- c(1,2,3,4,5)
# x2 <- c(5,6,7,8,9)
# setdiff(x1, x2)
# setdiff(x2, x1)
# x1 %>% setdiff(x2)
# x1 %>% setdiff(., x2)
# x1 %>% setdiff(x2, .)

#---G-or-W-or-All-------------------------------------------------

plot(density(as.matrix(tissue_df2)))
tissue_df <- tissue_df2

#----------------------------------------------------------------------

for (i in c("ts","tn","hk")){
  df_tp <-as.data.frame(matrix(0,ncol=ncol(tissue_df),nrow = nrow(tissue_df)))
  colnames(df_tp) = colnames(tissue_df)
  row.names(df_tp) <- paste(row.names(tissue_df),i,sep = "_")
  assign(i,df_tp)
}

#-get_tissue_specific-"ts"----------------------------------------------

# for (i in 1:ncol(tissue_df)){
#   tissue_df <- tissue_df[order(tissue_df[,i],decreasing = T),]
#   ts <- ts[row.names(tissue_df),]
#   if(tissue_df[1,i]> 3 & tissue_df[2,i]< 2){ts[1,i]<-5}}

for (i in 1:ncol(tissue_df)){
  tmp <- tissue_df[i]
  m <- max(tmp)
  m1 <- max(tmp[tmp!=max(tmp)])
  if(m>3 & m1< 2){ts[which.max(unlist(tmp)),i]<-5}}

#-get_tissue_entiched_but_not_specific-"tn"-----------------------------

# for (i in 1:ncol(tissue_df)){
#   tissue_df <- tissue_df[order(tissue_df[,i],decreasing = F),]
#   tn <- tn[row.names(tissue_df),]
#   for (v in 1:nrow(tissue_df)) { if(tissue_df[v,i]>2.5){tn[v:nrow(tissue_df),i]<-4}
#   }
# }

for (i in 1:ncol(tissue_df)){
  tmp <- tissue_df[i]
  t1 <- which(tmp>2.5)
  tn[t1,i] <- 4
}

#-get_house_keeping-"hk"------------------------------------------------

# df <- read.csv("//172.16.13.114/share/project/HBA_xiaoqi/20211209_20211213_yuque_Reanalysis/HBA_spnlib_matrix_after_preprocess_7536pro_2180sam_20211211.csv",header = T,sep = ",",stringsAsFactors = F)
# df <- df[,info2$SampleName]
# df <- as.data.frame(t(df))
# df$label <- 
# for (i in 1:ncol(tissue_df)){
# if(max(tissue_df[i])<0.6){hk[i]<-1}}
# hk[hk==0] <- NA
# hk <-hk[,apply(hk, 2, function(x) sum(!is.na(x)))>0]
# write.csv(hk,"house_keeping_max_0.6_4065pro_20210821.csv")
# hk_p <- intersect(names(hk),df_hk_all$X)
# write.csv(hk_p,"house_keeping_max_0.6_1062pro_20210821.csv")
# 
# hk <- read.csv("house_keeping_max_0.6_4065pro_20210821.csv")

#-sum-"ts"and"tn"-------------------------------------------------------

df_all <- rbind(ts,tn)
glabel <- sapply(row.names(df_all),function(x){gsub("_tn","",x)})
df_all$label <- sapply(glabel,function(x){gsub("_ts","",x)})
df_all[is.na(df_all)] <- 0

for (i in unique(df_all$label)) {
  data <- df_all[which(df_all$label==i),]
  data <- data[,-ncol(data)]
  data_max <- apply(data, 2,function(x){max(as.numeric(x),na.rm=T)} )
  df_all <- df_all[-which(df_all$label==i),]
  df_all <- rbind(df_all,c(data_max,i))
}

row.names(df_all) <- df_all$label
df_all <- df_all[,-ncol(df_all)]  # 0 means no groups
#df_all <- data.frame(as.matrix(df_all)%>%as.numeric())

df_all %<>% apply(c(1, 2), as.numeric) %>% as.data.frame()

df_all[df_all==0] <- NA
df_all_v <-df_all[,apply(df_all, 2, function(x) sum(!is.na(x)))>0] 
write.csv(df_all_v,"HBA_ts_9639pro_1521sam_Gray_Lobe_3_2_2.5_20230306.csv")

# df_all_t <- t(df_all) %>% as.data.frame()
# df_all_t[is.na(df_all_t)] <- 0

summary <- as.data.frame(matrix(ncol=2,nrow = nrow(df_all)),row.names = row.names(df_all))
names(summary) <- c("tissue_entiched_but_not_specific","tissue_specific")

for ( i in 4:5) {
  a<-function(x) sum(x==i,na.rm = T)
  num <- apply(df_all,1,a)
  summary[,i-3] <- num
}
sample_n <- as.data.frame(table(info2$Lobe))
summary$sample_n <- sample_n$Freq [match(row.names(summary),paste(sample_n$Var1,"G",sep = "_"))]
write.csv(summary,"HBA_ts_9639pro_1521sam_Gray_Lobe_3_2_2.5_20230306_summary.csv")
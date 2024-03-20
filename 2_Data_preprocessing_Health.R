rm(list = ls())
pacman::p_load(readxl,magrittr,
               RColorBrewer,
               Rtsne, umap,
               pheatmap,vioplot,
               ggpubr, ggplot2,
               corrplot, stringr,reshape2,preprocessCore)
require(magrittr)
require(plyr)
# Other used packages: readr, data.table, tibble, preprocessCore

# Step 1
t1 <- proc.time()
cat("1) File reading: \n"); print(proc.time() - t1)

pep_data <- data.table::fread('//172.16.13.136/share/members/xieyuting/HBA_combine_P1P2/combine_peptides_matrix_20220906.txt', sep = '\t', quote = '', check.names = F, stringsAsFactors = F, data.table = F, fill = T)

# pep_data[is.infinite(pep_data)] <- NA # assign infinity as NA
pep_data[is.na(pep_data)] <- NA
pep_data <- pep_data[complete.cases(pep_data[, 1]), ] # remove NA peptides
pep_data <- pep_data[!grepl("^1/CON", pep_data[, 2], fixed = F), ] # remove contaminants
pep_data[pep_data == 0] <- NA # assign zeros as NA


# Step 2
cat("2) log2 transformation: \n"); print(proc.time() - t1)
pep_data_log2 <- log2(pep_data[, 3:ncol(pep_data), drop = F]) %>% tibble::add_column(., prot = pep_data$prot, .before = 1)
rownames(pep_data_log2) <- pep_data[, 1]


# Step 3
cat("3) Quantile normalization: \n"); print(proc.time() - t1)
pep_data_log2_qn <- preprocessCore::normalize.quantiles(as.matrix(pep_data_log2[, 2:ncol(pep_data_log2)]))
colnames(pep_data_log2_qn) <- colnames(pep_data_log2)[-1]
rownames(pep_data_log2_qn) <- rownames(pep_data_log2)

data_tech_rep <- cbind(pep_data[, 1:2], pep_data_log2_qn)
#is.null(tech_rep_f)
#is.null(batchf)

rm(list = ls()[grep('^pep_data', ls())])


# Step 4
cat("4) Arrangement: \n"); print(proc.time() - t1)
data <- data_tech_rep
colnames(data)[1:2] <- c("tg", "prot")

n <- ncol(data)
pep2 <- apply(data[, -c(1, 2), drop = F], 1, function(x) { # log2 then mean of all files
  NAs <- length(which(is.na(x)))
  meanexpr1 <- sum(as.numeric(x), na.rm = TRUE) / (n - NAs)
  meanexpr2 <- sum(as.numeric(x), na.rm = TRUE) / n
  d <- c(NAs, meanexpr1, meanexpr2)
  return(d)
})
pep2 <- t(pep2)
colnames(pep2) <- c("NAs", "meanexpr1", "meanexpr2")
pep_expr <- cbind(data[, 1], pep2, data[, c(-1)])

#order by pg ,#NA,intesity
pep_order <- pep_expr[order(pep_expr[, 5], pep_expr[, 2], -pep_expr[, 3]), ]
colnames(pep_order)[1] <- "tg"
pep_order2 <- pep_order[, c(-2, -3, -4)]

rm(list = ls()[grep('^data', ls())])
rm(list = c('pep2', 'pep_expr', 'pep_order'))
gc() # Garbage Collection

# Step 5
cat("5) Top 3 peptides selection: \n"); print(proc.time() - t1)
pep_order2_dup1 <- pep_order2[duplicated(pep_order2$prot), ] #1764 protein only one peptide
pep_order2_dup2 <- pep_order2_dup1[duplicated(pep_order2_dup1$prot), ]
pep_order2_dup3 <- pep_order2_dup2[duplicated(pep_order2_dup2$prot), ]
pep_order2_top3 <- pep_order2[!(rownames(pep_order2) %in% rownames(pep_order2_dup3)), ]

rm(list = ls()[grep('^pep_order2_dup', ls())])


# Step 6
cat("6) Protein matrix generation: \n"); print(proc.time() - t1)
pep_order2_top3 <- pep_order2_top3[c("prot", "tg", colnames(pep_order2_top3)[3:ncol(pep_order2_top3)])]
pep_order2_top3[pep_order2_top3 == 0] <- NA

lr_top3 <- "top3"
if(lr_top3 == "top3"){ # mean of top 3
  top3_mean <- ddply(pep_order2_top3,
                     .variables = "prot",
                     .fun = function(df_sub){
                       mean_ls <- colMeans(df_sub[, -c(1, 2), drop = F], na.rm = T)
                       return(round(mean_ls, 2))
                     })
  
  readr::write_csv(top3_mean, 'prot_matrix_top3.csv', na = '')
  
}else{ # LR
  prot_matrix <- pep2prot(pep_order2_top3)
  prot_matrix <- prot_matrix[, -2]
  
  readr::write_csv(prot_matrix, 'prot_matrix_lr.csv', na = '')
}

#---------------------------------------------------
df <- top3_mean %>% tibble::column_to_rownames('prot')
# df <- data.table::fread('//172.16.13.114/share/project/HBA_xiaoqi/DataAnalysis20220712/DIANN_combine_part1_part2_20220804/prot_matrix_top3.csv', sep = ',', quote = '', check.names = F, stringsAsFactors = F, data.table = F, fill = T)
# row.names(df) <- df$prot
# df <- df[,-1]

info <- read.csv("D:/Guomics/HBA/Data_analysis20220901/FIleInfo_HBA_20220816_withoutFALSE_and_pool.csv",header = T,sep = ",",stringsAsFactors = F)
source('D:/Guomics/HBA/Data_analysis_record/20210517_DIANN_no_inference_test/result/Fragpipe/NA_deal.R')
#stringsAsFactors = F 

info1 <- info[which(info$Use=="TRUE"),]
# info2 <- info1[-which(info1$Lobe =="pool"|info1$pathology=="epilepsy"),]#2583 samples
# info2 <- info1[-which(info1$sample=="b28_406"|info1$sample=="b44_636"),]
info2 <- info1
# this two sample only identified 1700 proteins

df1 <- df[!(grepl(";", row.names(df))),] ## no protein group
names(df1)  <-  sapply(names(df1),function(x){str_split(x,"_30min_2")[[1]][1]})
names(df1)  <-  sapply(names(df1),function(x){str_split(x,"HBA_")[[1]][2]})


df2 <- df1[,info2$sample] # 2824 samples
# tmp <- df1[,grep("pool",names(df))]
# df2 <- cbind(df2,tmp)

tmp1 <- df2[apply(df2, 1, function(x) sum(!is.na(x)))>0,] #delete 0 NA in all samples
df3 <- tmp1[,apply(tmp1, 2, function(x) sum(!is.na(x)))>0] # 11746 protein *2824 samples
NA_threshold_table(df3)
write.csv(df3,"HBA_prot_matrix_11746pros_2824samples_202200909.csv")

#-----------------------------------------------
non_miss <- c()
for (i in unique(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions)){
  region <- info2$sample[which(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions ==i & info2$W_G =="G")]
  df4 <- df3[,region]
  non_miss_tmp <- row.names(df4)[which(apply(df4,1,function(x){sum(!is.na(x))>0.2*length(region)}))]
  if(i== unique(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions)[1]){non_miss=non_miss_tmp}else{non_miss <- union(non_miss,non_miss_tmp)}
}

for (i in unique(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions[which(info2$W_G=="W")])){
  region <- info2$sample[which(info2$BBA_Anatomical.and.modified.Cyto.architectonic.descriptions[which(info2$W_G=="W")] == i)]
  df4 <- df3[,region]
  non_miss_tmp <- row.names(df4)[which(apply(df4,1,function(x){sum(!is.na(x))>0.2*length(region)}))]
  non_miss <- union(non_miss,non_miss_tmp)
}

non_miss_epilepsy <- row.names(df3)[apply(df3[,which(info2$pathology =="epilepsy")],1,function(x){sum(!is.na(x)) > 0.2*89})]
non_miss_t <- union(non_miss,non_miss_epilepsy)

df5 <- df3[non_miss_t,]
NA_threshold_table(df5) #"100%" "9639" "55.02%"
write.csv(df5,"df5_HBA_peptide2protein_20220811.csv",quote = F)

#-pool-sample-----------------------------------------

df_pool <- df5[,grepl("pool",names(df5))]#154 pool sample
tmp2 <-df_pool[apply(df_pool, 1, function(x) sum(!is.na(x)))>0,]
df_pool2 <-tmp2[,apply(tmp2, 2, function(x) sum(!is.na(x)))>0]  #8213 proteins
write.csv(df_pool2,"df_pool_test20220809.csv",quote = F)
NA_threshold_table(df_pool2) #"100%" "8213" "40.59%"

#-----------------------------------------------------

# intensity <- data.frame(intensity=apply(df5,2,function(x){sum(x,na.rm = T)}),row.names = names(df5)) # min 41214.07 max 117387.1

df7 <- df5
df7[is.na(df7)]  <-  min(df7,na.rm = T)-1
df7_qn <- preprocessCore::normalize.quantiles(as.matrix(df7))
df8 <- as.data.frame(df7_qn)
names(df8) <- names(df7)
row.names(df8) <- row.names(df7)
info3 <- info2[-which(info2$pathology == "epilepsy"|grepl("be",info2$sample)),]
info4 <- info2[which(info2$pathology == "epilepsy"|grepl("be",info2$sample)),]
df_else <- df8[,info4$sample]
df8 <- df8[,info3$sample]
intensity <- data.frame(intensity=apply(df8,2,function(x){sum(x,na.rm = T)}),row.names = names(df8)) # min 119187.3 max 135700.8

#-Limma-------------------------------------------------

library(limma)

info3$batchID<-info3$sample
info3$batchID<-sapply(info3$batchID,function(x){str_split(x,"_")[[1]][1]})
rownames(info3)<-info3$sample
info3<-info3[colnames(df8),]

limma_df<-limma::removeBatchEffect(df8,batch=info3$batchID)

df9 <- cbind(limma_df, df_else)

write.csv(df9,"HBA_prot_matrix20220811.csv",quote = F)
write.csv(info2,"FIleInfo_HBA20220811.csv",quote = F)
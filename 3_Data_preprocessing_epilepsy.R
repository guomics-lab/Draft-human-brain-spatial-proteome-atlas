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

setwd("Z:/Project/HBA/epilepsy_re-libsearch_and_data_analysis20230703/combine_part1_part2")
pep_data <- data.table::fread('combine_peptide_matrix_epilepsy.txt', sep = '\t', quote = '', check.names = F, stringsAsFactors = F, data.table = F, fill = T)

pep_data[is.infinite(pep_data)] <- NA # assign infinity as NA
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

n <- ncol(data)  #ncol(data)-2 
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

  readr::write_csv(top3_mean, 'P2P_prot_matrix_top3.csv', na = '')

}else{ # LR
  prot_matrix <- pep2prot(pep_order2_top3)
  prot_matrix <- prot_matrix[, -2]
  
  readr::write_csv(prot_matrix, 'prot_matrix_lr.csv', na = '')
}

#----all epilepsy samples 20230705-------
# rm(list = ls())
setwd("Z:/Project/HBA/epilepsy_re-libsearch_and_data_analysis20230703/combine_part1_part2")
de <- read.csv('P2P_prot_matrix_top3.csv', row.names = 1)
de1 <- de[!(grepl(";", row.names(de))),] ## no protein group
names(de1)  <-  sapply(names(de1),function(x){str_split(x,"_30min_2")[[1]][1]})
names(de1)  <-  sapply(names(de1),function(x){str_split(x,"HBA_")[[1]][2]})
# write.csv(de1,"P2P_epilepsy_protein_matrix_with_NA.csv")
de1 <- de1[,-grep("pool", colnames(de1))]
# de3 <- de2[!apply(de2, 1, function(x){sum(is.na(x))>0.8*ncol(de2)}), ]

setwd("Z:/Project/HBA/epilepsy_re-libsearch_and_data_analysis20230703")
# df <- read.csv("df5_HBA_peptide2protein_20220811.csv", row=1) 
df <- data.table::fread('HBA_prot_matrix_11746pros_2824samples_202200909.csv', stringsAsFactors = F, data.table = F, fill = T)
info <- read.csv("FIleInfo_HBA_20220816_withoutFALSE_and_pool_update_20230706.csv")
rownames(df) <- df[,1]
df <- df[,-1]

sample <- data.frame(samples=colnames(df))
sample$pathology <- info$pathology[match(sample$samples, info$sample)]
sample1 <- sample[-grep("epilepsy", sample$pathology),]
# sample1 <- sample1[-grep("be", sample1$samples),]
sample1 <- na.omit(sample1) #delete pool samples

df1 <- df[,sample1$samples]
df1$pro <- rownames(df1)

# ProHE <- intersect(rownames(de1), rownames(df1))
de2 <- de1  #[ProHE, ]
# de2 <- de2[, -grep("pool", colnames(de2))] #don't delete pool samples
de2$pro <- rownames(de2)

df_combine <- merge(df1, de2, by="pro", all=T)
rownames(df_combine) <- df_combine[,1] 
df_combine <- df_combine[,-1]
df_combine1 <- df_combine[apply(df_combine, 1, function(x){sum(!is.na(x))>0}), ]

setwd("Z:/Project/HBA/data_and_code_20230726/1_library_quantify")
write.csv(df_combine1, "HBA_protein_matrix_after_p2p_11761pro_2681_samples_without_pool_20230727.csv")

#normalize.quantiles
df_combine1[is.na(df_combine1)]  <-  min(df_combine1,na.rm = T)-1
df_combine_nor <- preprocessCore::normalize.quantiles(as.matrix(df_combine1))
df_combine_nor <- as.data.frame(df_combine_nor)
rownames(df_combine_nor) <- rownames(df_combine1)
colnames(df_combine_nor) <- colnames(df_combine1)

datENor <- df_combine_nor[,c(grep("be", colnames(df_combine_nor)), grep("b112_e55", colnames(df_combine_nor)))]
datENor$pro <- rownames(datENor)

setwd("D:/1_guomics_xyt/1Project/hba_research_workflow/9_DataAnalysis_20220901/WGCNA_analysis")
HNor<-data.table::fread('HBA_prot_matrix_normolized_20220811.csv',data.table = F)
rownames(HNor)<-HNor[,1]
HNor<-HNor[,-1]

HNor$pathology <- info$pathology[match(rownames(HNor), info$sample)]
HNor <- HNor[grep("health", HNor$pathology), ] 
HNor1 <- HNor[,-ncol(HNor)]
HNor2 <- as.data.frame(t(HNor1))
HNor2$pro <- rownames(HNor2)

NewE_OldH <- merge(datENor, HNor2, by="pro") #9638*2681 
setwd("Z:/Project/HBA/epilepsy_re-libsearch_and_data_analysis20230703/combine_part1_part2")
write.csv(NewE_OldH, "HBA_pro_matrix_new_epilepsy_nor_data_old_healthy_nor_data_without_pool.csv")

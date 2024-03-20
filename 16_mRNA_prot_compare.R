rm(list = ls())
options(stringsAsFactors = F)
#install.packages("pacman")
pacman::p_load(readxl,magrittr,
               RColorBrewer,
               Rtsne, umap,
               pheatmap,vioplot,
               ggpubr, ggplot2,
               corrplot, stringr,gridExtra)

#-Read------------------------------------------------------------

setwd("D:/Guomics/HBA/Data_analysis20220901/ts_score/")

mRNA <- read.csv("D:/Guomics/HBA/Data_analysis20220901/mRNA/HBA_Allen_ts_18415prot_2113sample_Gray_BBA_3_2_2.5_20221011.csv",stringsAsFactors = F)
row.names(mRNA) <- mRNA$X
mRNA <- mRNA[,-1]

protein <- read.csv("HBA_ts_9639pro_1521sam_Grey_BBA_3_2_2.5_20220921.csv",stringsAsFactors = F)
row.names(protein) <- sapply(protein$X,function(x){str_sub(x,1,length(x)-4)})
protein <- protein[,-1]

protein_sub <- read.csv("HBA_Allen_ts_8929pro_1521sam_Gray_BBA_3_2_2.5_20220922.csv",stringsAsFactors = F)
row.names(protein_sub) <- sapply(protein_sub$X,function(x){str_sub(x,1,length(x)-4)})
protein_sub <- protein_sub[,-1]

#------------------------------------------------------------
# protein <- protein_sub
  
protein <- protein[which(row.names(protein)%in%row.names(mRNA)),]
mRNA <- mRNA[which(row.names(mRNA)%in%row.names(protein)),]
protein <- protein[row.names(mRNA),]
names(protein) <- sapply(names(protein),function(x){strsplit(x,"\\.")[[1]][1]})
protein_sub <- protein_sub[row.names(mRNA),]

gene2prot <- read.csv("D:/Guomics/HBA/Data_analysis20220901/mRNA/gene2protein_reviewed_uniprot.csv")
gene2prot <- tidyr::separate_rows(gene2prot, 'Gene.Names', sep = ' ')
# protein <- protein[,which(names(protein)%in%Pro2geneSymbol$From)]
names(mRNA) <- gene2prot$Entry[match(names(mRNA),gene2prot$Gene.Names)]


int <- as.data.frame(matrix(nrow = nrow(protein),ncol = 3))
row.names(int) <- row.names(protein)
names(int)[1:3] <-c("prot","mRNA","inter") 

# length(unique(Pro2geneSymbol$To[which(Pro2geneSymbol$From %in%names(protein))]))

for (i in 1:nrow(int)) {
  d_m <- t(mRNA[i,])
  m_n <- row.names(d_m)[which(!is.na(d_m))]
  d_p <- t(protein[i,])
  p_n <- row.names(d_p)[which(!is.na(d_p))]
  int[i,1] <- length(setdiff(p_n,m_n))
  int[i,2] <- length(setdiff(m_n,p_n))
  int[i,3] <- length(intersect(m_n,p_n))
}

int_name <- list()
for (i in row.names(int)) {
  d_m <- t(mRNA[i,])
  m_n <- row.names(d_m)[which(!is.na(d_m))]
  d_p <- t(protein[i,])
  p_n <- row.names(d_p)[which(!is.na(d_p))]
  int_name[[c(paste0(i,"-prot"))]] <- setdiff(p_n,m_n)
  int_name[[c(paste0(i,"-mRNA"))]] <- setdiff(m_n,p_n)
  int_name[[(paste0(i,"-inter"))]]  <- intersect(m_n,p_n)
}
int_name <- as.data.frame(unlist(int_name))
int_name$type <- sapply(row.names(int_name),function(x){strsplit(x,"-")[[1]][2]})
int_name$type <- gsub("\\d","",int_name$type)
int_name$region <- sapply(row.names(int_name),function(x){strsplit(x,"-")[[1]][1]})

setwd("D:/Guomics/HBA/Data_analysis20220901/ts_score/pie_charts/")
write.csv(int,"TS_score_mRNA_protein_intersect20220923.csv")
write.csv(int_name,"TS_score_mRNA_protein_intersect_list20220923.csv")

#--pie-plot-----------------------------------------------

library(ggplot2)
library(dplyr)

int_pie <- int
# row.names(int_pie) <- sapply(row.names(int_pie),function(x){strsplit(x,"\\,")[[1]][2]})
# Sys.setlocale('LC_ALL','C')
row.names(int_pie) <- gsub("\\/","",row.names(int_pie))
size <- apply(int_pie,1,sum)

setwd("D:/Guomics/HBA/Data_analysis20220901/ts_score/pie_charts/")
for (i in row.names(int_pie)) {
  data_tmp <- as.data.frame(t(int_pie[i,]))
  names(data_tmp) <- "value"
  row.names(data_tmp) <- paste(row.names(data_tmp),data_tmp$value,sep = "_")
  data_tmp <- data_tmp %>% 
    arrange(desc(row.names(data_tmp))) %>%
    mutate(prop = value / sum(data_tmp$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  pdf(paste0(i,"_PieGraph.pdf"),width = size[i]/100+2,height = size[i]/100+2)
  p=ggplot(data_tmp, aes(x="", y=prop, fill=row.names(data_tmp))) +
    geom_bar(stat="identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=c("#e63a46", "#cae3d0","#4c84a8"))+ 
    # inter/mRNA/prot
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = row.names(data_tmp)), color = "white", size=6) 
    
  print(p)
  dev.off()
}

#--all-in-one------------------------------------------------------------
pl <- list()
# int_pie$sum <- apply(int_pie,1,sum)
int_pie$size = -(int_pie$sum/300+5)+14
for (i in row.names(int_pie)) {
  data_tmp <- as.data.frame(t(int_pie[i,c(1,2,3)]))
  size_use <- int_pie[i,5]
  names(data_tmp) <- "value"
  row.names(data_tmp) <- paste(row.names(data_tmp),data_tmp$value,sep = "_")
  data_tmp <- data_tmp %>% 
    arrange(desc(row.names(data_tmp))) %>%
    mutate(prop = value / sum(data_tmp$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  p=ggplot(data_tmp, aes(x="", y=prop, fill=row.names(data_tmp))) +
    geom_bar(stat="identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=c("#e63a46", "#cae3d0","#4c84a8"))+ 
    # inter/mRNA/prot
    ggtitle(i)+
    # geom_text(aes(y = ypos, label = row.names(data_tmp)), color = "white", size=6) +
    theme_void() + 
    theme(legend.position="none",
    plot.margin=unit(c(size_use,size_use,size_use,size_use), "mm"),
    plot.title = element_text(size = 2)) 
    pl[[i]] <- ggplotGrob(p)
}

pdf('test.pdf',width = 4,height = 50)
grid.arrange(grobs=pl,ncol=2)
dev.off()


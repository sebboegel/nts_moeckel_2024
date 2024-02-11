library("tximport")
library("DESeq2")
library("ggplot2")
library(dplyr)
library("tidyr")
library(tidyverse)
library(reshape2)
library(stringr)
#BiocManager::install('EnhancedVolcano')
#library(EnhancedVolcano)
#library(readr)

##################################
#load biomart database 


#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "mmusculus_gene_ensembl",
#                         host = "ensembl.org")
#ttg <- biomaRt::getBM(
#  attributes = c("ensembl_transcript_id", "transcript_version",
#                 "ensembl_gene_id", "external_gene_name", "description",
#                 "transcript_biotype"),
#  mart = mart)

###################################

dir="."

#1.) differential expression analysis Baff-R vs control
samples <- read.table(file.path(dir,"metadata_baff-r.txt"), header=TRUE) #Baff-R vs control


rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")
#load biomart gene ID vonversion table
ttg=readRDS(file.path(dir,"biomart.rds"))
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))
#read DESEQ2 results
txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)


abundance_table=round(txi$abundance,0)
colnames(abundance_table)=samples$run_accession
write.csv(abundance_table,file=file.path(dir,"BAFF-R_vs_WT.csv"))
rownames(samples) <- samples$alias
#perform differential expression analysis
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~strain )
dds=DESeq(ddsTxi)
res <- results(dds,contrast = c("strain","BAFF-R-ko","WT"))
write.csv(res,file=file.path(dir,"BAFF-R_WT-diff_exp.csv"))


#2.) differential expression analysis Baff vs control

samples <- read.table(file.path(dir,"metadata_baff.txt"), header=TRUE) #Baff vs control

rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")
#load biomart gene ID vonversion table
ttg=readRDS(file.path(dir,"biomart.rds"))
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))
#read DESEQ2 results
txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)


abundance_table=round(txi$abundance,0)
colnames(abundance_table)=samples$run_accession
write.csv(abundance_table,file=file.path(dir,"BAFF_vs_WT.csv"))

#perform differential expression analysis
ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design = ~strain )
dds=DESeq(ddsTxi)
res <- results(dds,contrast = c("strain","BAFF-ko","WT"))
write.csv(res,file=file.path(dir,"BAFF_WT-diff_exp.csv"))

###3.) Figures
samples <- read.table(file.path(dir,"metadata_baff_baff-r.txt"), header=TRUE) #all samples


rownames(samples) <- samples$run_accession
files=file.path(dir, samples$run_accession, "abundance.h5")
#load biomart gene ID vonversion table
ttg=readRDS(file.path(dir,"biomart.rds"))
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id','ext_gene'))
#read DESEQ2 results
txi <- tximport(files, type="kallisto", tx2gene=ttg, ignoreTxVersion=TRUE)

abundance_table=round(txi$abundance,0)
colnames(abundance_table)=samples$run_accession

all_mouse <- cbind(rownames(abundance_table), data.frame(abundance_table, row.names=NULL))
colnames(all_mouse)[1]="gene"

m <- melt(all_mouse)
m <- m %>% 
  mutate(variable = str_replace_all(variable, "\\.", "-"))
m <- na.omit(m)
colnames(samples)=c("variable","strain")
test=merge(m,samples,by="variable")


#Figure 2
gene_set = c("Col3a1")
test$facet = factor(test$gene, levels = gene_set)

dataMedian_Baff <- summarise(group_by(test[test$strain=="BAFF-ko",], gene), MD = median(value))
dataMedian_BaffR <- summarise(group_by(test[test$strain=="BAFF-R-ko",], gene), MD = median(value))
dataMedian_WT <- summarise(group_by(test[test$strain=="WT",], gene), MD = median(value))

dataMedian_Baff[dataMedian_Baff$gene=="Col3a1",]
dataMedian_BaffR[dataMedian_BaffR$gene=="Col3a1",]
dataMedian_WT[dataMedian_WT$gene=="Col3a1",]

ylabels=c("BAFF ko","BAFF-R ko","wt")
test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=strain,y=value,fill=factor(strain,levels = c("BAFF-ko", "WT" , "BAFF-R-ko"),labels=c("BAFF ko", "wt","BAFF-R ko")))) +
  geom_boxplot(outlier.shape = NA) + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=2,strip.position="right",scales = "free") +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=1,strip.position="right") +
  scale_x_discrete(labels= ylabels)

#Figure 3
gene_set = c("Gpx3","Igfbp7","Ccn2","Kap","Umod","Ren1","Txnip")
test$facet = factor(test$gene, levels = gene_set)
ylabels=c("BAFF ko","BAFF-R ko","wt")
test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=strain,y=value,fill=factor(strain,levels = c("BAFF-ko", "WT" , "BAFF-R-ko"),labels=c("BAFF ko", "wt","BAFF-R ko")))) +
  geom_boxplot(outlier.shape = NA) + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=2,strip.position="right",scales = "free") +
  scale_x_discrete(labels= ylabels)

#Figure 4
gene_set = c("Tnfrsf12a","Tnfsf12")
test$facet = factor(test$gene, levels = gene_set)
ylabels=c("BAFF ko","BAFF-R ko","wt")
test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=strain,y=value,fill=factor(strain,levels = c("BAFF-ko", "WT" , "BAFF-R-ko"),labels=c("BAFF ko", "wt","BAFF-R ko")))) +
  geom_boxplot(outlier.shape = NA) + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=2,strip.position="right",scales = "free") +
  scale_x_discrete(labels= ylabels)

#Figure S1
gene_set = c("Tnfsf13b","Tnfrsf13c","Tnfrsf17","Tnfrsf13b")
test$facet = factor(test$gene, levels = gene_set)
ylabels=c("BAFF ko","BAFF-R ko","wt")
test %>% 
  filter(gene %in% gene_set) %>%
  ggplot(aes(x=strain,y=value,fill=factor(strain,levels = c("BAFF-ko", "WT" , "BAFF-R-ko"),labels=c("BAFF ko", "wt","BAFF-R ko")))) +
  geom_boxplot(outlier.shape = NA) + 
  labs(fill = "strain",y="TPM") +  
  geom_point(position=position_jitterdodge()) +
  theme_bw(base_size = 16) + facet_wrap(~facet,ncol=2,strip.position="right",scales = "free") +
  scale_x_discrete(labels= ylabels)

library(DESeq2)
library(dplyr)

count_table <- read.csv2("counts.txt",header=TRUE,sep='\t')
count_table_oneSample_samnii <- count_table %>% select(lst of samples) %>% filter(str_detect(Geneid, "CONPKAHK")) #changes for each taxon
rownames(count_table_oneSample_samnii) <- count_table_oneSample_samnii[,1]
count_table_oneSample_samnii <- count_table_oneSample_samnii[,-1]
count_table_oneSample_samnii <-  count_table_oneSample_samnii[rowSums(count_table_oneSample_samnii) > 1000 , ]
count_table_oneSample_samnii <-  count_table_oneSample_samnii[, colSums(count_table_oneSample_samnii) > 1000]
col_data_oneSample <- read.csv2("term_status.txt",header=FALSE,sep='\t')
colnames(col_data_oneSample) <- c("sampleID","condition")
col_data_oneSample_samnii <- col_data_oneSample %>% filter(sampleID %in% names(count_table_oneSample_samnii))
rownames(col_data_oneSample_samnii) <- col_data_oneSample_samnii[,1]
col_data_oneSample_samnii <- col_data_oneSample_samnii %>%  select(condition)
dataset_oneSample_samnii <- DESeqDataSetFromMatrix(countData = count_table_oneSample_samnii,colData = col_data_oneSample_samnii, design = ~condition)
dataset_oneSample_samnii$condition <- factor(dataset_oneSample_samnii$condition, levels = c("fullterm","preterm"))
dds_oneSample_samnii <- DESeq(dataset_oneSample_samnii)
res_oneSample_samnii <- results(dds_oneSample_samnii)
write.table(res_oneSample_samnii,"deseq_out.txt",sep="\t")

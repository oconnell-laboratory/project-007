library(data.table)
library(DESeq2)

###read in count data###

raw_counts<-fread(file="whole_blood_raw_counts.csv")
ensembl_id<-raw_counts$Ensembl_ID
gene_name<-raw_counts$Gene
gene_biotype<-raw_counts$Gene_biotype
raw_counts$Ensembl_ID<-NULL
raw_counts$Gene<-NULL
raw_counts$Gene_biotype<-NULL
raw_counts<-as.matrix(raw_counts)
rownames(raw_counts)<-ensembl_id

###generate scaling factors###

total_read_count<-apply(raw_counts, 2, sum)
scaling_factors<-total_read_count/min(total_read_count)

###filter for protein coding genes###

protein_coding_index<-which(gene_biotype=="protein_coding")
filtered_raw_counts<-raw_counts[protein_coding_index,]
ensembl_id<-ensembl_id[protein_coding_index]
gene_name<-gene_name[protein_coding_index]
gene_biotype<-gene_biotype[protein_coding_index]

###read in sample meta data###

meta_data<-read.csv(file="whole_blood_sample_meta_data.csv", stringsAsFactors=F)
rownames(meta_data)<-meta_data$Sample_ID
meta_data$Race<-as.factor(meta_data$Race)
meta_data$Ethnicity<-as.factor(meta_data$Ethnicity)
meta_data$Diabetes<-as.factor(meta_data$Diabetes)

###filter count data and meta data for Hispanic and White non-Hispanic samples###

group_1_index<-which(meta_data$Ethnicity=="Hispanic")
group_2_index<-which(meta_data$Ethnicity=="Not_hispanic")
group_2_index<-group_2_index[which(meta_data$Race[group_2_index]=="White")]
filtered_sample_index<-c(group_1_index, group_2_index)
meta_data<-meta_data[filtered_sample_index,]
filtered_raw_counts<-filtered_raw_counts[,filtered_sample_index]
scaling_factors<-scaling_factors[filtered_sample_index]

###scale and center continuous meta data###

meta_data$Age<-scale(meta_data$Age, center=T, scale=T)
meta_data$Actual_neutrophil_count<-scale(meta_data$Actual_neutrophil_count, center=T, scale=T)
meta_data$Actual_lymphocyte_count<-scale(meta_data$Actual_lymphocyte_count, center=T, scale=T)
meta_data$Actual_monocyte_count<-scale(meta_data$Actual_monocyte_count, center=T, scale=T)
meta_data$Actual_eosinophil_count<-scale(meta_data$Actual_eosinophil_count, center=T, scale=T)
meta_data$Inferred_neutrophil_count<-scale(meta_data$Inferred_neutrophil_count, center=T, scale=T)
meta_data$Inferred_lymphocyte_count<-scale(meta_data$Inferred_lymphocyte_count, center=T, scale=T)
meta_data$Inferred_monocyte_count<-scale(meta_data$Inferred_monocyte_count, center=T, scale=T)
meta_data$Inferred_eosinophil_count<-scale(meta_data$Inferred_eosinophil_count, center=T, scale=T)

###differential expression uncorrected for cell counts###

DEseqDS<-DESeqDataSetFromMatrix(countData=filtered_raw_counts, colData=meta_data, design = ~ Age + Diabetes + Ethnicity)

sizeFactors(DEseqDS)<-scaling_factors

DEseqDS<-DESeq(DEseqDS)

results<-results(DEseqDS, independentFiltering=F)
results<-cbind(ensembl_id, gene_name, gene_biotype, results)
results<-results[order(results$pvalue),]

write.csv(results, file="uncorrected_de_results.csv")

###differential expression corrected by actual cell counts###

DEseqDS<-DESeqDataSetFromMatrix(countData=filtered_raw_counts, colData=meta_data, design = ~ Age + Diabetes + Actual_neutrophil_count + Actual_lymphocyte_count + Actual_monocyte_count + Actual_eosinophil_count + Ethnicity)

sizeFactors(DEseqDS)<-scaling_factors

DEseqDS<-DESeq(DEseqDS)

results<-results(DEseqDS, independentFiltering=F)
results<-cbind(ensembl_id, gene_name, gene_biotype, results)
results<-results[order(results$pvalue),]

write.csv(results, file="actual_count_corrected_de_results.csv")

###differential expression corrected by inferred cell counts###

DEseqDS<-DESeqDataSetFromMatrix(countData=filtered_raw_counts, colData=meta_data, design = ~ Age + Diabetes + Inferred_neutrophil_count + Inferred_lymphocyte_count + Inferred_monocyte_count + Inferred_eosinophil_count + Ethnicity)

sizeFactors(DEseqDS)<-scaling_factors

DEseqDS<-DESeq(DEseqDS)

results<-results(DEseqDS, independentFiltering=F)
results<-cbind(ensembl_id, gene_name, gene_biotype, results)
results<-results[order(results$pvalue),]

write.csv(results, file="inferred_count_corrected_de_results.csv")

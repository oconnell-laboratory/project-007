###IMPORTANT: modified CIBERSORT source code must be loaded first###

library(data.table)

###define funtions###

unity_scaling<-function(x){
x_min<-min(x)
x_range<-max(x)-min(x)
(x-x_min)/x_range
}

rpm_norm<-function(x){
(x/sum(x))*1000000
}

###read in whole blood count data###

raw_counts<-fread(file="whole_blood_raw_counts.csv")
mix_ensembl_id<-raw_counts$Ensembl_ID
raw_counts$Ensembl_ID<-NULL
raw_counts$Gene<-NULL
raw_counts$Gene_biotype<-NULL
raw_counts<-as.matrix(raw_counts)
rownames(raw_counts)<-mix_ensembl_id

###normalize raw counts###

rpm_values<-apply(raw_counts, 2, rpm_norm)

###read in CIBERSORT reference signature matrix###

cibersort_reference_sig<-read.csv(file="cibersort_reference_signature_matrix.csv", stringsAsFactors=F)
ref_ensembl_id<-cibersort_reference_sig$Ensembl_ID
cibersort_reference_sig$Ensembl_ID<-NULL
cibersort_reference_sig<-as.matrix(cibersort_reference_sig)
rownames(cibersort_reference_sig)<-ref_ensembl_id

####read in WBC4.126 marker gene list###

WBC4_126_marker_genes<-read.csv(file="WBC4_126_marker_genes.csv", header=T, stringsAsFactors=F)
marker_ensembl_id<-WBC4_126_marker_genes$Ensembl_ID
cell_type<-WBC4_126_marker_genes$Cell_type

###read in actual leukocyte counts###

actual_cell_counts<-read.csv(file="actual_cell_counts.csv", header=T, stringsAsFactors=F)

###correlation between inferred and actual counts via each pipline###

WBC4_126_marker_rpm_values<-rpm_values[match(marker_ensembl_id, mix_ensembl_id),]
collapsed_data<-WGCNA::collapseRows(WBC4_126_marker_rpm_values, rowGroup=cell_type, rowID=marker_ensembl_id, method="ME", connectivityBasedCollapsing=FALSE)
WBC4_126_inferred_cell_counts<-collapsed_data$datETcollapsed
WBC4_126_inferred_cell_counts<-t(WBC4_126_inferred_cell_counts)
WBC4_126_inferred_cell_counts<-apply(WBC4_126_inferred_cell_counts, 2, unity_scaling)
WBC4_126_inferred_cell_counts<-as.data.frame(WBC4_126_inferred_cell_counts)

pca_accuracy<-c(cor(actual_cell_counts$Neutrophil_count, WBC4_126_inferred_cell_counts$Neutrophil, method="spearman"),
cor(actual_cell_counts$Lymphocyte_count, WBC4_126_inferred_cell_counts$Lymphocyte, method="spearman"),
cor(actual_cell_counts$Monocyte_count, WBC4_126_inferred_cell_counts$Monocyte, method="spearman"),
cor(actual_cell_counts$Eosinophil_count, WBC4_126_inferred_cell_counts$Eosinophil, method="spearman")
)

cibersort_inferred_cell_counts<-CIBERSORT(cibersort_reference_sig, rpm_values, perm=0, QN=F)
cibersort_inferred_cell_counts<-as.data.frame(cibersort_inferred_cell_counts)

cibersort_accuracy<-c(cor(actual_cell_counts$Neutrophil_count, cibersort_inferred_cell_counts$Neutrophil, method="spearman"),
cor(actual_cell_counts$Lymphocyte_count, cibersort_inferred_cell_counts$Lymphocyte, method="spearman"),
cor(actual_cell_counts$Monocyte_count, cibersort_inferred_cell_counts$Monocyte, method="spearman"),
cor(actual_cell_counts$Eosinophil_count, cibersort_inferred_cell_counts$Eosinophil, method="spearman")
)

###generate bootstrap samples###

n_boot<-1000

neutrophil<-rep(NA, n_boot)
lymphocyte<-rep(NA, n_boot)
monocyte<-rep(NA, n_boot)
eosinophil<-rep(NA, n_boot)

pca_boot_accuracy<-data.frame(neutrophil, lymphocyte, monocyte, eosinophil) 
cibersort_boot_accuracy<-data.frame(neutrophil, lymphocyte, monocyte, eosinophil) 

for(i in 1:n_boot){

sample_index<-sample(1:dim(rpm_values)[2], dim(rpm_values)[2], replace=T)
sample_rpm_values<-rpm_values[,sample_index]
sample_actual_cell_counts<-actual_cell_counts[sample_index,]

sample_WBC4_126_marker_rpm_values<-sample_rpm_values[match(marker_ensembl_id, mix_ensembl_id),]
colnames(sample_WBC4_126_marker_rpm_values)<-NULL
collapsed_data<-WGCNA::collapseRows(sample_WBC4_126_marker_rpm_values, rowGroup=cell_type, rowID=marker_ensembl_id, method="ME", connectivityBasedCollapsing=FALSE)
sample_WBC4_126_inferred_cell_counts<-collapsed_data$datETcollapsed
sample_WBC4_126_inferred_cell_counts<-t(sample_WBC4_126_inferred_cell_counts)
sample_WBC4_126_inferred_cell_counts<-apply(sample_WBC4_126_inferred_cell_counts, 2, unity_scaling)
sample_WBC4_126_inferred_cell_counts<-as.data.frame(sample_WBC4_126_inferred_cell_counts)

pca_boot_accuracy[i,]<-c(cor(sample_actual_cell_counts$Neutrophil_count, sample_WBC4_126_inferred_cell_counts$Neutrophil, method="spearman"),
cor(sample_actual_cell_counts$Lymphocyte_count, sample_WBC4_126_inferred_cell_counts$Lymphocyte, method="spearman"),
cor(sample_actual_cell_counts$Monocyte_count, sample_WBC4_126_inferred_cell_counts$Monocyte, method="spearman"),
cor(sample_actual_cell_counts$Eosinophil_count, sample_WBC4_126_inferred_cell_counts$Eosinophil, method="spearman")
)

sample_cibersort_inferred_cell_counts<-CIBERSORT(cibersort_reference_sig, sample_rpm_values, perm=0, QN=F)
sample_cibersort_inferred_cell_counts<-as.data.frame(sample_cibersort_inferred_cell_counts)

cibersort_boot_accuracy[i,]<-c(cor(sample_actual_cell_counts$Neutrophil_count, sample_cibersort_inferred_cell_counts$Neutrophil, method="spearman"),
cor(sample_actual_cell_counts$Lymphocyte_count, sample_cibersort_inferred_cell_counts$Lymphocyte, method="spearman"),
cor(sample_actual_cell_counts$Monocyte_count, sample_cibersort_inferred_cell_counts$Monocyte, method="spearman"),
cor(sample_actual_cell_counts$Eosinophil_count, sample_cibersort_inferred_cell_counts$Eosinophil, method="spearman")
)

}

###calculate bootstrap statistics###

pca_95_low<-c(quantile(pca_boot_accuracy$neutrophil, probs=0.025),
quantile(pca_boot_accuracy$lymphocyte, probs=0.025),
quantile(pca_boot_accuracy$monocyte, probs=0.025),
quantile(pca_boot_accuracy$eosinophil, probs=0.025)
)
 
pca_95_hi<-c(quantile(pca_boot_accuracy$neutrophil, probs=0.975),
quantile(pca_boot_accuracy$lymphocyte, probs=0.975),
quantile(pca_boot_accuracy$monocyte, probs=0.975),
quantile(pca_boot_accuracy$eosinophil, probs=0.975)
)

cibersort_95_low<-c(quantile(cibersort_boot_accuracy$neutrophil, probs=0.025),
quantile(cibersort_boot_accuracy$lymphocyte, probs=0.025),
quantile(cibersort_boot_accuracy$monocyte, probs=0.025),
quantile(cibersort_boot_accuracy$eosinophil, probs=0.025)
)

cibersort_95_hi<-c(quantile(cibersort_boot_accuracy$neutrophil, probs=0.975),
quantile(cibersort_boot_accuracy$lymphocyte, probs=0.975),
quantile(cibersort_boot_accuracy$monocyte, probs=0.975),
quantile(cibersort_boot_accuracy$eosinophil, probs=0.975)
)

boot_accuracy_difference<-pca_boot_accuracy-cibersort_boot_accuracy

p_value<-c(min(length(which(boot_accuracy_difference$neutrophil>0))/length(boot_accuracy_difference$neutrophil), 
length(which(boot_accuracy_difference$neutrophil<0))/length(boot_accuracy_difference$neutrophil))*2,
min(length(which(boot_accuracy_difference$lymphocyte>0))/length(boot_accuracy_difference$lymphocyte), 
length(which(boot_accuracy_difference$lymphocyte<0))/length(boot_accuracy_difference$lymphocyte))*2,
min(length(which(boot_accuracy_difference$monocyte>0))/length(boot_accuracy_difference$monocyte), 
length(which(boot_accuracy_difference$monocyte<0))/length(boot_accuracy_difference$monocyte))*2,
min(length(which(boot_accuracy_difference$eosinophil>0))/length(boot_accuracy_difference$eosinophil), 
length(which(boot_accuracy_difference$eosinophil<0))/length(boot_accuracy_difference$eosinophil))*2
)

results<-data.frame(pca_accuracy, pca_95_low, pca_95_hi, cibersort_accuracy, cibersort_95_low, cibersort_95_hi, p_value)
row.names(results)<-c("neutrophils", "lymphocytes", "monocytes", "eosinophils")
write.csv(results, file="pipeline_comparision_results.csv")

###IMPORTANT: modified CIBERSORT source code must be loaded first###

library(data.table)

###read in whole blood count data###

raw_counts<-fread(file="whole_blood_raw_counts.csv")
mix_ensembl_id<-raw_counts$Ensembl_ID
raw_counts$Ensembl_ID<-NULL
raw_counts$Gene<-NULL
raw_counts$Gene_biotype<-NULL
raw_counts<-as.matrix(raw_counts)
rownames(raw_counts)<-mix_ensembl_id

###normalize raw counts###

rpm_norm<-function(x){
(x/sum(x))*1000000
}

rpm_values<-apply(raw_counts, 2, rpm_norm)

###read in reference signature matrix###

reference_sig<-read.csv(file="cibersort_reference_signature_matrix.csv", stringsAsFactors=F)
ref_ensembl_id<-reference_sig$Ensembl_ID
reference_sig$Ensembl_ID<-NULL
reference_sig<-as.matrix(reference_sig)
rownames(reference_sig)<-ref_ensembl_id

###perform decovolution###

inferred_cell_counts<-CIBERSORT(reference_sig, rpm_values, perm=100, QN=F)
inferred_cell_counts<-as.data.frame(inferred_cell_counts)
write.csv(inferred_cell_counts, file="cibersort_inferred_cell_counts.csv")

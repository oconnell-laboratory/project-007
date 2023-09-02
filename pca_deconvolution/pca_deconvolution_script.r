library(data.table)
library(WGCNA)

###read in count data###

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

###filter normalized counts for marker genes###

marker_genes<-read.csv(file="WBC4_126_marker_genes.csv", header=T, stringsAsFactors=F)

marker_ensembl_id<-marker_genes$Ensembl_ID
cell_type<-marker_genes$Cell_type
marker_rpm_values<-rpm_values[match(marker_ensembl_id, mix_ensembl_id),]

###generate inferred cell counts###

collapsed_data<-collapseRows(marker_rpm_values, rowGroup=cell_type, rowID=marker_ensembl_id, method="ME", connectivityBasedCollapsing=FALSE)
inferred_cell_counts<-collapsed_data$datETcollapsed
inferred_cell_counts<-t(inferred_cell_counts)

unity_scaling<-function(x){
x_min<-min(x)
x_range<-max(x)-min(x)
(x-x_min)/x_range
}

inferred_cell_counts<-apply(inferred_cell_counts, 2, unity_scaling)
inferred_cell_counts<-as.data.frame(inferred_cell_counts)
write.csv(inferred_cell_counts, file="pca_inferred_cell_counts.csv")

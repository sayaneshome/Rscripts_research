library(Biobase)
library(preprocessCore)
library("pheatmap")
library("RColorBrewer")
library('xlsx')
d <- read.xlsx("dataset_cult_BM.xlsx",1)
#this row is empty

d <- d[-1538:-1660,]

colnames(cult_BM_miRNA_norm) <- names(d)[-1]
row.names(cult_BM_miRNA_norm) <- d[,1]



#names can't be converted to numeric so pull them out
whole_names1 <- d[,1]
whole_names2 <- colnames(d)[-1]
#transpose will now be numeric
cult_BM_miRNA_matrix <- as.matrix.data.frame(t(d[,-1]))

#row ids become column names
colnames(cult_BM_miRNA_matrix) <- trnames1
row.names(cult_BM_miRNA_matrix) <- whole_names2

#Compare this with read.csv("mirnaclean.csv")
str(cult_BM_miRNA_matrix)

#now you can normalize
cult_BM_miRNA_matrix_norm <- normalize.quantiles(cult_BM_miRNA_matrix[2:5,2:1537])
colnames(cult_BM_miRNA_matrix_norm) <- trnames1
row.names(cult_BM_miRNA_matrix_norm) <- whole_names2

write.csv(miRNA_norm, "mirna_norm.csv")
write.csv(miRNA_transp_norm, "mirna_transp_norm.csv")

test <- read.delim('mc10_goodhits_o6.out',header = FALSE,sep = )
test2 <- read.delim('release94.accession2geneid',header = FALSE,sep = '\t')
test1 <- read.delim('All_Data.gene_info',header = TRUE,sep = '\t')
colnames(test2)[3] <- "refseq_accession"
colnames(test)[2] <- "refseq_accession"
test2$V4 <- NULL
results_file <- merge(test,test2,by = "refseq_accession")
colnames(results_file)[14] <- "GeneID"
colnames(results_file)[13] <- "tax_id"
results_file1 <- merge(results_file,test1,by = "GeneID")
write.csv(results_file1,'result_output.csv')
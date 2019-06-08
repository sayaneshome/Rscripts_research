test <- read.delim('mc_test.out',header = FALSE,sep = '\t')
test2 <- read.delim('release94.accession2geneid',header = FALSE,sep = '\t')
test1 <- read.delim('All_Data.gene_info',header = TRUE,sep = '\t')
test2$V3 <- NULL
colnames(test2)[3] <- "refseq_accession"
colnames(test)[2] <- "refseq_accession"

test <- test[order(test$refseq_accession),]
test2 <- test2[order(test2$refseq_accession),]

results_file <- merge(test,test2,by = "refseq_accession")
write.csv(results_file2,'results_chk.csv')
colnames(results_file)[14] <- "GeneID"
colnames(results_file)[13] <- "tax_id"
results_file1 <- merge(results_file,test1,by = "GeneID")
write.csv(results_file1,'result_output.csv')
#test_trinity <- read.delim('trinity_contigs_inreference.sam',header = FALSE,sep = '\t')
#test_tr <- test_trinity[-c(1),]
#test_tr1 <- data.frame(test_tr$V1,test_tr$V2,test_tr$V3,test_tr$V4)
#test_tr2 <- data.frame(test_tr1[grep("TRINITY", test_tr1$test_tr.V1),])
#variant <- read.delim('trinity_var.txt',header = TRUE,sep = '\t')
#colnames(variant)[2] <- "position_in_contig"
#colnames(variant)[1] <- "Trinity_ID"
#colnames(test_tr2)[1] <- "Trinity_ID"
#colnames(test_tr2)[4] <- "Starting_position_inreferencegenome"
#results_file2 <- merge(variant,test_tr2,by = "Trinity_ID")
#results_file2$variantposition_reference <- results_file2$position_in_contig+results_file2$Starting_position_inreferencegenome
#write.csv(results_file2,'variants_with_trinityreference.csv')
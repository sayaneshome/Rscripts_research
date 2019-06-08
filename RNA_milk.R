rawData_otherRNA<-read.xlsx("otherRNA_foranalysis.xlsx",1);
View(rawData_otherRNA)
count(rawData_otherRNA, "antisense_RNA")
sum(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
nrow(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
sum(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "antisense_RNA",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "lincRNA",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "protein_coding",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "processed_transcript",])
View(compared_gene_targets)
View(rawData)

#
#split multiple-mapped miRNAs listed in X1 column into seperate rows
raw_mirna_norm_updated <- strsplit(raw_mirna_norm$X1, split = "/")
data.frame(X1 = unlist(raw_mirna_norm_updated),Raw = rep(raw_mirna_norm$Raw, sapply(raw_mirna_norm_updated,length)))
raw_mirna_norm_seperated <- data.frame(X1 = unlist(raw_mirna_norm_updated),Raw = rep(raw_mirna_norm$Raw, sapply(raw_mirna_norm_updated,length)))
View(raw_mirna_norm_seperated)

#rename the column name from X1 to miRNA
colnames(raw_mirna_norm_seperated)[1] <- "miRNA"

#import the compared gene target list as excelsheet.
library("xlsx")
library(stringr)
dat = readLines("stable.txt")
dat = as.data.frame(do.call(rbind, strsplit(dat, split=" {2,10}")))
] 
#FINAL CODE USED;kept the above lines just in case if required 
#will remove all miRNAs after '/' character.
raw_mirna_norm$X1 <- sub('(/).*$', '', raw_mirna_norm$X1, perl=TRUE)
#filter all pir into new dataframe
pir_Dataset <- subset(raw_mirna_norm[grep("piR", raw_mirna_norm$X1), ],)
#import as excelsheet
write.xlsx(pir_Dataset, "pirDataset.xlsx")
#filter all miRNAs into new dataframe
miRNA_Dataset <- raw_mirna_norm[-grep("piR", raw_mirna_norm$X1), ]
#import as excelsheet
write.xlsx(miRNA_Dataset, "miRNADataset.xlsx")

#rename column name X1 to miRNA
colnames(miRNA_Dataset)[1] <- "miRNA"
#compare miRNAs from our data with miRTarbase and list details by using natural join

#import the details obtained of miRNAs from the database as excelsheet
write.xlsx(compared_gene_targets, "gene_targets_miRNA.xlsx")
#list all miRNAs that are novel from our data using 


#Link for reference : http://www.dummies.com/programming/r/how-to-use-the-merge-function-with-data-sets-in-r/
require('dplyr')
novel_miRNA <- anti_join(miRNA_Dataset,miRTarBase_MTI)
write.xlsx(novel_miRNA,"novel_miRNA.xlsx")

#yet to be reframed.
rawData_otherRNA<-read.xlsx("otherRNA_foranalysis.xlsx",1);
View(rawData_otherRNA)
count(rawData_otherRNA, "antisense_RNA")
sum(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
nrow(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
sum(which(rawData_otherRNA$Gene.type == "antisense_RNA"))
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "antisense_RNA",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "lincRNA",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "protein_coding",])
nrow(rawData_otherRNA[rawData_otherRNA$Gene.type == "processed_transcript",])
View(compared_gene_targets)
View(rawData)

pir_reference1 <- read.table("sperm_sncRNA_expression_table_reference.txt")
pir_reference1$V5 <- sub('\\|.*','',pir_reference1$V5)
View(pir_reference1)

pir_reference2 <- read.table("PiR_base.txt")




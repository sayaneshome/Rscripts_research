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



#t-test analysis









Control <- read_xlsx('Raw_wholemilk_normalised.xlsx')
Past_homo_Skim <- read_xlsx('Past_homo_skim_normalised.xlsx')
Homo_HT_Skim <- read_xlsx('Homo_HT_Skim_normalised.xlsx')
Post_Skim <- read_xlsx('Post_cream_normalised.xlsx')
Raw_HT_Skim <- read_xlsx('Raw_HT_skim_normalised.xlsx')
SC_BM <- read_xlsx('SC_BM_normalised.xlsx')
Control <- Control[, -c(3:6)]

Past_homo_Skim_vsControl <- merge(Control,Past_homo_Skim,by="miRNA")
Homo_HT_Skim_vsControl <- merge(Control,Homo_HT_Skim,by="miRNA")
Post_Skim_vsControl <- merge(Control,Post_Skim,by="miRNA")
SC_BMvsControl <- merge(Control,SC_BM,by="miRNA")
Raw_HTSKim_vsControl <- merge(Control,Raw_HT_Skim,by="miRNA")

SC_BMvsControl1 <- SC_BMvsControl[-grep("piR", SC_BMvsControl$miRNA), ]
Raw_HTSKim_vsControl1 <- Raw_HTSKim_vsControl[-grep("piR", Raw_HTSKim_vsControl$miRNA), ]
Post_SKim_vsControl1 <- Post_Skim_vsControl[-grep("piR", Post_Skim_vsControl$miRNA), ]
Past_homo_Skim_vsControl1 <- Past_homo_Skim_vsControl[-grep("piR", Past_homo_Skim_vsControl$miRNA), ]
Homo_HT_Skim_vsControl1 <- Homo_HT_Skim_vsControl[-grep("piR", Homo_HT_Skim_vsControl$miRNA), ]


t.test(Homo_HT_Skim_vsControl1$Control_R1, Homo_HT_Skim_vsControl1$Homo_HTSK_R1, paired = TRUE, alternative = "two.sided")
t.test(Homo_HT_Skim_vsControl1$Control_R2, Homo_HT_Skim_vsControl1$Homo_HTSK_R2, paired = TRUE, alternative = "two.sided")
t.test(Homo_HT_Skim_vsControl1$Control_R3, Homo_HT_Skim_vsControl1$Homo_HTSK_R3, paired = TRUE, alternative = "two.sided")
t.test(Homo_HT_Skim_vsControl1$Control_R4, Homo_HT_Skim_vsControl1$Homo_HTSK_R4, paired = TRUE, alternative = "two.sided")
t.test(Past_homo_Skim_vsControl1$Control_R1, Past_homo_Skim_vsControl1$Past_hoSK_R1, paired = TRUE, alternative = "two.sided")
t.test(Past_homo_Skim_vsControl1$Control_R2, Past_homo_Skim_vsControl1$Past_hoSK_R2, paired = TRUE, alternative = "two.sided")
t.test(Past_homo_Skim_vsControl1$Control_R3, Past_homo_Skim_vsControl1$Past_hoSK_R3, paired = TRUE, alternative = "two.sided")
t.test(Past_homo_Skim_vsControl1$Control_R4, Past_homo_Skim_vsControl1$Past_hoSK_R4, paired = TRUE, alternative = "two.sided")
t.test(Post_SKim_vsControl1$Control_R1, Post_SKim_vsControl1$Post_SK_R1, paired = TRUE, alternative = "two.sided")
t.test(Post_SKim_vsControl1$Control_R2, Post_SKim_vsControl1$Post_SK_R2, paired = TRUE, alternative = "two.sided")
t.test(Post_SKim_vsControl1$Control_R3, Post_SKim_vsControl1$Post_SK_R3, paired = TRUE, alternative = "two.sided")
t.test(Post_SKim_vsControl1$Control_R4, Post_SKim_vsControl1$Post_SK_R4, paired = TRUE, alternative = "two.sided")
View(Raw_HTSKim_vsControl1)
t.test(Raw_HTSKim_vsControl1$Control_R1, Raw_HTSKim_vsControl1$Raw_HTSK_R1, paired = TRUE, alternative = "two.sided")
t.test(Raw_HTSKim_vsControl1$Control_R2, Raw_HTSKim_vsControl1$Raw_HTSK_R2, paired = TRUE, alternative = "two.sided")
t.test(Raw_HTSKim_vsControl1$Control_R3, Raw_HTSKim_vsControl1$Raw_HTSK_R3, paired = TRUE, alternative = "two.sided")
t.test(Raw_HTSKim_vsControl1$Control_R4, Raw_HTSKim_vsControl1$Raw_HTSK_R4, paired = TRUE, alternative = "two.sided")
View(SC_BMvsControl1)
t.test(SC_BMvsControl1$Control_R1, SC_BMvsControl1$SCBM_R1, paired = TRUE, alternative = "two.sided")
t.test(SC_BMvsControl1$Control_R2, SC_BMvsControl1$SCBM_R2, paired = TRUE, alternative = "two.sided")
t.test(SC_BMvsControl1$Control_R3, SC_BMvsControl1$SCBM_R3, paired = TRUE, alternative = "two.sided")
t.test(SC_BMvsControl1$Control_R4, SC_BMvsControl1$SCBM_R4, paired = TRUE, alternative = "two.sided")
write.csv(Past_homo_Skim_vsControl1,"Control_past_homo_skim.csv")
write.csv(Post_SKim_vsControl1,"Control_post_skim.csv")
write.csv(Raw_HTSKim_vsControl1,"Control_raw_HTSKIM.csv")
write.csv(SC_BMvsControl1,"Control_SC_BM.csv")





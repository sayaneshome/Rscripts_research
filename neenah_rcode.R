neenah_chipseq <- read_xlsx('ChIPseq raw data.xlsx',2)
neenah_rnaseq <- read_xlsx('RNAseq data_edger3_rn6.xlsx',2)
neenah_chipseq$...12 <- NULL
neenah_chipseq$...11 <- NULL
neenah_chipseq$pro <- NULL
neenah_chipseq$npro <- NULL
neenah_chipseq1 <- neenah_chipseq[!is.na(neenah_chipseq$genes),]
neenah_rnaseq$...13 <- NULL
neenah_rnaseq$FC <- NULL
neenah_rnaseq$FDR...14 <- NULL
colnames(neenah_rnaseq)[1] <- "genes"
neenah_chipseq1$genes <- sub(".*0\\|(.*)\\|[+|-]", "\\1",neenah_chipseq1$genes,perl = TRUE)
neenah_rnaseq$genes <- sub(".*0\\|(.*)\\|[+|-]", "\\1",neenah_rnaseq$genes,perl = TRUE)
neenah_chipseq1$chr <- NULL
neenah_chipseq1$start <- NULL
neenah_chipseq1$end <- NULL
chipseq_rnaseq_merge <- merge(neenah_chipseq1,neenah_rnaseq,by = 'genes')
chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'decr' & chipseq_rnaseq_merge$expr_CG == 'down'), ]
down_both_chipseq_rnaseq <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'decr' & chipseq_rnaseq_merge$expr_CG == 'down'), ]
down_both_chipseq_rnaseq <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'incr' & chipseq_rnaseq_merge$expr_CG == 'up'), ]
down_both_chipseq_rnaseq <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'decr' & chipseq_rnaseq_merge$expr_CG == 'down'), ]
up_both_chipseq_rnaseq <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'incr' & chipseq_rnaseq_merge$expr_CG == 'up'), ]
chipseq_incr_rnaseq_down <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'incr' & chipseq_rnaseq_merge$expr_CG == 'down'), ]
chipseq_decr_rnaseq_up <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'decr' & chipseq_rnaseq_merge$expr_CG == 'up'), ]
write.xlsx(down_both_chipseq_rnaseq,'Chipseq_rnaseq_both_Downregulated_tfam.xlsx')
write.csv(down_both_chipseq_rnaseq,'Chipseq_rnaseq_both_Downregulated_tfam.csv')
write.csv(up_both_chipseq_rnaseq,'Chipseq_rnaseq_both_UPregulated_tfam.csv')
write.csv(chipseq_incr_rnaseq_down,'Chipseq_up_rnaseq_down_regulated_tfam.csv')
write.csv(chipseq_decr_rnaseq_up,'Chipseq_down_rnaseq_up_regulated_tfam.csv')
chipseq_rnaseq_nchg <- chipseq_rnaseq_merge[(chipseq_rnaseq_merge$change == 'nchg' & chipseq_rnaseq_merge$expr_CG == 'nchg'), ]
write.csv(chipseq_rnaseq_nchg,'Chipseq_rnaseq_nochange_tfam.csv')
write.csv(neenah_rnaseq,'rnaseq_Results_simplified_.csv')
write.csv(neenah_chipseq1,'chipseq_Results_simplified_tfam.csv')
write.csv(list_pyri_rnaseq,'pyri_rnaseq_list_for_go.csv')
list_tfam_rnaseq <- neenah_rnaseq$genes
write.csv(list_tfam_rnaseq,'tfam_rnaseq_list_for_go.csv')
treatment_1 <- read.table('pantherChart (5).txt',sep = "")
treatment_1 <- read.table('pantherChart (5).txt')
treatment_1 <- read_xlsx('pyri_pathway_go.xlsx')
treatment_12 <- read_xlsx('tfam_pathway_go.xlsx')
treatment_2 <- read_xlsx('tfam_pathway_go.xlsx')
View(treatment_1)
View(treatment_2)
colnames(treatment_2)[1] <- 'Signalling pathways'
treatment_12 <- merge(treatment_1,treatment_2,by = 'Signalling pathways')
View(treatment_12)
write.csv(treatment_12,'signalling_pathways_tfam_pyri.csv')

plot(chipseq_rnaseq_merge$logFC,-log10(as.numeric(chipseq_rnaseq_merge$PValue)),pch = 20,xlab = 'FoldChange',ylab = 'P-value',xlim = c(-6,6))
with(subset(chipseq_rnaseq_merge, logFC > 0), points(logFC,-log10(as.numeric(PValue)),pch = 20, col="red"))
with(subset(chipseq_rnaseq_merge,logFC < 0), points(logFC,-log10(as.numeric(PValue)),pch = 20, col="blue"))
with(subset(chipseq_rnaseq_merge,genes=='Tnxa-ps1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)


with(subset(chipseq_rnaseq_merge,genes=='Adipor1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Arhgap24'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Cd24'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Dync2li1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Edem3'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Fam129a'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Fgf5'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Fxyd6'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Glrx2'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Lnc134'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Rnf43'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Swt1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Samd5'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Tmem47'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Tnxa-ps1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="yellow"))
with(subset(chipseq_rnaseq_merge,genes=='Arhgap8'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Asns'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Ctgf'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Eya4'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Fbn1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Fjx1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Fst'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Gja1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Lyn'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Mnda'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Myo7a'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Pccb'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Prr5'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Prrx1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Prtfdc1'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Rassf9'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))
with(subset(chipseq_rnaseq_merge,genes=='Sema3d'),points(logFC,-log10(as.numeric(PValue)),pch = 20, col="orange"))




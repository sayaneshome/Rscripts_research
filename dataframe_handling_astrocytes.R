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
write.csv(neenah_rnaseq,'rnaseq_Results_simplified_tfam.csv')
write.csv(neenah_chipseq1,'chipseq_Results_simplified_tfam.csv')


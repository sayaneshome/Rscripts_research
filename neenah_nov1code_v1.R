#newcode for VC plots
setwd("~/Desktop/neenah_analysis_ensemble")
chip_u_rna_u_pyri <- read.csv('Chipseq_rnaseq_both_UPregulated.csv')
chip_u_rna_u_tfam <- read.csv('Chipseq_rnaseq_both_UPregulated_tfam.csv')
chip_d_rna_d_pyri <- read.csv('Chipseq_rnaseq_both_Downregulated.csv')
chip_d_rna_d_tfam <- read.csv('Chipseq_rnaseq_both_Downregulated_tfam.csv')
chip_d_rna_u_tfam <- read.csv('Chipseq_down_rnaseq_up_regulated_tfam.csv')
chip_d_rna_u_pyri <- read.csv('Chipseq_down_rnaseq_up_regulated.csv')
chip_u_rna_d_pyri <- read.csv('Chipseq_up_rnaseq_down_regulated.csv')
chip_u_rna_d_tfam <- read.csv('Chipseq_up_rnaseq_down_regulated_tfam.csv')
chip_rna_nc_pyri <- read.csv('Chipseq_rnaseq_nochange.csv')
chip_rna_nc_tfam <- read.csv('Chipseq_rnaseq_nochange_tfam.csv')
plot(chip_rna_nc_pyri$logFC,-log10(as.numeric(chip_rna_nc_pyri$PValue)),pch = 20,xlab = 'Log2 FoldChange',ylab = '-Log10 p-value',col = 'grey',xlim = c(-8,6),ylim = c(0,100))
points(chip_u_rna_u_pyri$logFC,-log10(as.numeric(chip_u_rna_u_pyri$PValue)),col = "red")
points(chip_d_rna_u_pyri$logFC,-log10(as.numeric(chip_d_rna_u_pyri$PValue)),col = "red")
points(chip_d_rna_d_pyri$logFC,-log10(as.numeric(chip_d_rna_d_pyri$PValue)),col = "blue")
points(chip_u_rna_d_pyri$logFC,-log10(as.numeric(chip_u_rna_d_pyri$PValue)),col = "blue")

plot(chip_rna_nc_tfam$logFC,-log10(as.numeric(chip_rna_nc_tfam$PValue)),pch = 20,xlab = 'Log2 FoldChange',ylab = '-Log10 p-value',col = 'grey',xlim = c(-8,6),ylim = c(0,100))
points(chip_u_rna_u_tfam$logFC,-log10(as.numeric(chip_u_rna_u_tfam$PValue)),col = "red")
points(chip_d_rna_u_tfam$logFC,-log10(as.numeric(chip_d_rna_u_tfam$PValue)),col = "red")
points(chip_d_rna_d_tfam$logFC,-log10(as.numeric(chip_d_rna_d_tfam$PValue)),col = "blue")
points(chip_u_rna_d_tfam$logFC,-log10(as.numeric(chip_u_rna_d_tfam$PValue)),col = "blue")

pyri_list <- do.call("rbind", list(chip_d_rna_d_pyri,chip_u_rna_u_pyri,chip_u_rna_d_pyri,chip_d_rna_u_pyri,chip_rna_nc_pyri))
tfam_list <- do.call("rbind", list(chip_d_rna_d_tfam,chip_u_rna_u_tfam,chip_u_rna_d_tfam,chip_d_rna_u_tfam,chip_rna_nc_tfam))

counts <- table(go_graph$PANTHER.Pathways,go_graph$Genecount_T1, go_graph$Genecount_T2)
barplot(counts, main="Test",xlab="Number of Gears", col=c("darkblue","red"),legend = rownames(counts))




pdb_native <- read.pdb('p53_native_prod_stride500.pdb')
pdb_c141_m133 <- read.pdb('p53_c141y_m133l_stride500.pdb')
pdb_native <- read.pdb('p53_native_prod_stride500.pdb',multi = TRUE)
pdb_c141_m133 <- read.pdb('p53_c141y_m133l_stride500.pdb',multi = TRUE)
pdb_c141 <- read.pdb('p53_c141y_stride500_full.pdb',multi = TRUE)
pdb_g266 <- read.pdb('p53_g266e_stride500.pdb',multi = TRUE)
pdb_g266_r267 <- read.pdb('p53_g266e_r267h_stride500.pdb')
pdb_g266_r267 <- read.pdb('p53_g266e_r267h_stride500.pdb',multi = TRUE)
pdb_r158 <- read.pdb('p53_r158h_stride500_full.pdb',multi = TRUE)
ca.inds_p53 <- atom.select(pdb_native,elety = 'CA')
ca.inds_p53_c141 <- atom.select(pdb_c141,elety = 'CA')
ca.inds_p53_c141_m133 <- atom.select(pdb_c141_m133,elety = 'CA')
ca.inds_p53_g266 <- atom.select(pdb_g266,elety = 'CA')
ca.inds_p53_g266_r267 <- atom.select(pdb_g266_r267,elety = 'CA')
ca.inds_p53_r158 <- atom.select(pdb_r158,elety = 'CA')
xyz_p53 <- pdbfit(pdb_native,inds = ca.inds_p53)
xyz_p53_c141 <- pdbfit(pdb_c141,inds = ca.inds_p53_c141)
xyz_p53_c141_m133 <- pdbfit(pdb_c141_m133,inds = ca.inds_p53_c141_m133)
xyz_p53_g266 <- pdbfit(pdb_g266,inds = ca.inds_p53_g266)
xyz_p53_g266_r267 <- pdbfit(pdb_g266_r267,inds = ca.inds_p53_g266_r267)
xyz_p53_r158 <- pdbfit(pdb_r158,inds = ca.inds_p53_r158)
par(mfrow=c(2,3))
rf <- data.frame(rmsf(xyz_p53[,ca.inds_p53$xyz]))
rf$residue_no <- rownames(rf)
rf$residue_no <- (as.numeric(rf$residue_no))+100
colnames(rf)[1] <- "RMSF_values"
plot(rf$residue_no,rf$RMSF_values, ylab="RMSF", xlab="Residue Position_p53", typ="l")
rf3 <- data.frame(rmsf(xyz_p53_c141_m133[,ca.inds_p53_c141_m133$xyz]))
rf3$residue_no <- rownames(rf3)
rf3$residue_no <- (as.numeric(rf3$residue_no))+100
colnames(rf3)[1] <- "RMSF_values"
plot(rf3$residue_no,rf3$RMSF_values, ylab="RMSF", xlab="Residue Position_p53_c141_m133", typ="l")
rf5 <- data.frame(rmsf(xyz_p53_g266_r267[,ca.inds_p53_g266_r267$xyz]))
rf5$residue_no <- rownames(rf5)
rf5$residue_no <- (as.numeric(rf5$residue_no))+100
colnames(rf5)[1] <- "RMSF_values"
plot(rf5$residue_no,rf5$RMSF_values, ylab="RMSF", xlab="Residue Position_g266_r267", typ="l")
rf6 <- data.frame(rmsf(xyz_p53_r158[,ca.inds_p53_r158$xyz]))
rf6$residue_no <- rownames(rf6)
rf6$residue_no <- (as.numeric(rf6$residue_no))+100
colnames(rf6)[1] <- "RMSF_values"
plot(rf6$residue_no,rf6$RMSF_values, ylab="RMSF", xlab="Residue Position_r158", typ="l")
rf2 <- data.frame(rmsf(xyz_p53_c141[,ca.inds_p53_c141$xyz]))
rf2$residue_no <- rownames(rf2)
rf2$residue_no <- (as.numeric(rf2$residue_no))+100
colnames(rf2)[1] <- "RMSF_values"
plot(rf2$residue_no,rf2$RMSF_values, ylab="RMSF", xlab="Residue Position_c141", typ="l")
rf4 <- data.frame(rmsf(xyz_p53_g266[,ca.inds_p53_g266$xyz]))
rf4$residue_no <- rownames(rf4)
rf4$residue_no <- (as.numeric(rf4$residue_no))+100
colnames(rf4)[1] <- "RMSF_values"
plot(rf4$residue_no,rf4$RMSF_values, ylab="RMSF", xlab="Residue Position_g266", typ="l")

par(mfrow=c(2,3))
cij<-dccm(xyz_p53[,ca.inds_p53$xyz],cutoff.cij=0.3)
test<-cna(cij,cm=cm)
cm<- cmap(xyz_p53[,ca.inds_p53$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test,pdb_native)


cij1<-dccm(xyz_p53_c141[,ca.inds_p53_c141$xyz],cutoff.cij=0.3)
test1<-cna(cij1,cm=cm)
cm1<- cmap(xyz_p53_c141[,ca.inds_p53_c141$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test1, pdb_c141)

cij2<-dccm(xyz_p53_c141_m133[,ca.inds_p53_c141_m133$xyz],cutoff.cij=0.3)
test2<-cna(cij2,cm=cm)
cm2<- cmap(xyz_p53_c141_m133[,ca.inds_p53_c141_m133$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test2, pdb_c141_m133)

cij3<-dccm(xyz_p53_g266[,ca.inds_p53_g266$xyz],cutoff.cij=0.3)
test3<-cna(cij3,cm=cm)
cm3<- cmap(xyz_p53_g266[,ca.inds_p53_g266$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test3, pdb_g266)

cij4<-dccm(xyz_p53_g266_r267[,ca.inds_p53_g266_r267$xyz],cutoff.cij=0.3)
test4<-cna(cij4,cm=cm)
cm4<- cmap(xyz_p53_g266_r267[,ca.inds_p53_g266_r267$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test4, pdb_g266_r267)

cij5<-dccm(xyz_p53_r158[,ca.inds_p53_r158$xyz],cutoff.cij=0.3)
test5<-cna(cij5,cm=cm)
cm5<- cmap(xyz_p53_r158[,ca.inds_p53_r158$xyz],dcut=10, pcut=0.75, mask.lower=FALSE)
plot(test5, pdb_r158)




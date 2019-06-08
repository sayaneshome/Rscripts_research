                                                                                                   


 for (i in 1:52){
  a <- torsion.xyz(pdb_tauro$xyz[1,],1)
  b <- torsion.xyz(pdb_tauro$xyz[2,],1)
  d <- wrap.tor(a-b)
  difference_in_torsion <- rbind(d)
  i = i+1
}
Torsion_change <- t(difference_in_torsion)

jpeg(filename = 'tauro_native_transport.jpeg',units = "in",width = 12 ,height = 8,pointsize = 25,res = 200)
plot(Torsion_change,xlab = 'Atom_residue_index',typ='h',col='purple',lwd=5)
axis(side = 1,lwd = 2)
axis(side = 2, lwd = 2)
dev.off()
p1 <- read.csv('go_plot.csv')
COLS = c("Biological_Processes","Genecount","fold.Enrichment","P.value")
df1 = data.frame(p1[,1:4])
colnames(df1) = COLS
df1$Treatment = "Tfam"
df2 = data.frame(p1[,c(1,5:7)])
colnames(df2) = COLS
df2$Treatment = "Pyri"
df = rbind(df1,df2)
df <- arrange(df, desc(Genecount), Treatment)
df$Biological_Processes <- factor(df$Biological_Processes, levels=as.character(unique(df$Biological_Processes)))

ggplot(df,aes(x=Treatment,y=Genecount,fill=-log10(P.value))) + 
  geom_col(position="dodge",width=0.4) +
  coord_flip() + scale_fill_viridis(trans='log10',option="B")+
  facet_grid(Biological_Processes~.)+
  theme(strip.text.y = element_text(angle = 0))
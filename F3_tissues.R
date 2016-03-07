#Tissue specific work
setwd("~/Desktop/HorseAnnotation")
library(ggplot2)
#read in dataset with expressino values for each term and a dataset with annotation for each gene
annotated <- read.csv("Annotated_GeneTable.csv", header=T)
tissues <- read.csv("allTissues_geneTPM.csv", header=T)
#merge and match to get expression values from one and gene names from the other
m <- merge(annotated, tissues, by.x="gene.ID",by.y="geneName" )
m_subset <- subset(m, select=c("gene.ID","BrainStem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord","chr","strand","ref_gene_id.y"))
write.csv(m_subset, file = "demoMerge.csv")

#assessing fraction of total transcriptional output for each tissue
t_out <- aggregate(m_subset[,sapply(m_subset,is.numeric)],m_subset["chr"],sum)
write.csv(t_out, file="chr_sum.csv")
#calculated totals in excel: total sum and total nuclear
#now separate chr1-32 from chrM and divide their sums by the total
calc <- read.csv("chr_sum.csv", header=T)
totals <- calc[37:38,1:9]
totals <- t(totals)
totals <- totals[2:9,1:2]
colnames(totals)[colnames(totals)==c("37","38")] <- c("nuclear","mt")

#making the graph
library(reshape2)
library(dplyr)
#reshaping data for geom_bar
df.melt <- melt(totals, id="chr")
bar <- group_by(df.melt, Var2, Var1)
#plotting as stacked geom_bar
ggplot(bar, aes(x=Var1, y=value, fill=factor(Var2)))+
  geom_bar(position="stack", stat="identity")

#adjusting legends and labels
ggplot(bar, aes(x=Var1, y=value, fill=factor(Var2)))+
  geom_bar(position="stack", stat="identity") + 
  ylab("proportion of transcriptional output") + xlab("tissue") +
  guides(fill=guide_legend(title="gene origin",
                           labels=c("nuclear", "mt")))


tissue <- c("BrainStem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord")





library(ggplot2)
setwd("~/Desktop/HorseAnnotation")
#making genes vs isoforms scatter
data<-read.csv("HT_isoformsvsGenes.csv", header=T)
ggplot(data,aes(Number_genes,Number_isoforms)) + geom_point(aes(color=Tissue)) + xlab("Number of genes") + ylab("Number of Isoforms")

#Making genes vs isoforms histo
data2<-read.csv("HT_isoformsvsGenes_histo.csv", header=T)
#for unique isoform ratio, normalized
ggplot(data2, aes(x=Tissue,y=Normalized_total_isoforms_noU)) + geom_bar(stat="identity", aes(fill=Normalized_total_Uisoforms), position = "stack") + ylab("Normalized total isoforms")
#for non normalized unique isoform ratio
ggplot(data2, aes(x=Tissue,y=X.isoforms)) + geom_bar(stat="identity",aes(fill=X.Uisoforms),position = "stack")
#for uniquely present and absent isoforms
data3<-read.csv("HT_isoformsvsGenes_histo2.csv", header=T)
ggplot(data3, aes(x=Tissue,y=Normalized_Uisoforms)) + geom_bar(stat="identity",aes(fill=Normalized_absent_Uisoforms),position = "stack") + xlab("Tissue") + ylab("Normalized Unique Isoforms")

#more isoform plots
data4<-read.csv("HT_isoforms.csv", header=T)
keeps <- c("Nisoforms", "NUisoforms", "Absent_Uisoforms","Tissue")
isoforms <- c("Nisoforms", "NUisoforms", "Absent_Uisoforms")
newdata4 <- data4[keeps]
isoform_type <- newdata4[isoforms]
Nisoforms <- "Nisoforms"
Nisoforms_data <- subset(newdata4,select=c(Nisoforms,Tissue))
NUisoforms <- "NUisoforms"
NUisoforms_data <- subset(newdata4,select=c(NUisoforms, Tissue))
Absent_Uisoforms <- "Absent_Uisoforms"
Absent_Uisoforms_data <- subset(newdata4,select=c(Absent_Uisoforms,Tissue))
#need to reshape data to plot it in geom_bar
library(reshape2)
library(dplyr)
df.melt <- melt(newdata4, id="Tissue")
bar <- group_by(df.melt, variable, Tissue)
ggplot(bar, aes(x=Tissue, y=value, fill=factor(variable)))+
  geom_bar(position="stack", stat="identity")
#can't see NUisoforms and Absent_Uisoforms with above code
dat1 <- subset(bar,value >= 0)
dat2 <- subset(bar,value < 0)
ggplot() + 
  geom_bar(data = dat1, aes(x=Tissue, y=value, fill=factor(variable)),stat = "identity") +
  geom_bar(data = dat2, aes(x=Tissue, y=value, fill=factor(variable)),stat = "identity") +
  scale_fill_brewer(type = "seq", palette = 16) + ylab("Number of isoforms")
#Can see absent_Uisoforms from above code, but NUisoforms is too small
#Another way to make same plot:
ggplot() +
  geom_bar(data=Nisoforms_data, aes(x=Tissue,y=Nisoforms),stat = "identity") +
  geom_bar(data=NUisoforms_data, aes(x=Tissue,y=NUisoforms, fill=factor(NUisoforms),stat = "identity") +
  geom_bar(data=Absent_Uisoforms_data, aes(x=Tissue,y=Absent_Uisoforms, color=Absent_Uisoforms),stat = "identity") +
  scale_fill_brewer(type = "seq", palette = 2) + ylab("Number of isoforms")

# making new unique isoforms into heatmap format for better visualization and with absent isoforms
ggplot(data2, aes(x=Tissue,y=X.isoforms)) + geom_bar(stat="identity",aes(fill=X.Uisoforms),position = "stack") + 
  geom_bar(data=Absent_Uisoforms_data, aes(x=Tissue,y=Absent_Uisoforms),stat = "identity") + ylab("Number of isoforms")

#To change legend
ggplot() +
  geom_bar(data=data2, aes(x=Tissue,y=X.isoforms, fill=X.Uisoforms), stat="identity") +
  geom_bar(data=Absent_Uisoforms_data, aes(x=Tissue,y=Absent_Uisoforms),stat = "identity") + 
  ylab("Number of isoforms") + guides(fill=guide_legend(title="Number of Unique Isoforms"))

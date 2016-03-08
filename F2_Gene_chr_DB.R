setwd("~/Desktop/HorseAnnotation")
library(ggplot2)
###transcripts per chrs
data_chr<- read.csv("Annotated_GeneTable.csv", header=T)

#To try and subset and stack u transcipts in with total transcripts
library(plyr)
library(ggplot2)
ggplot(subset(data_chr,class_code.x.1 %in% c("j","=", "u"))) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + guides(fill=guide_legend(title="class code"))
#genes for NCBI
keeps_NCBI <- c("major_iso_id.x", "chr", "class_code.x.1")
newdata_chr_NCBI <- data_chr[keeps_NCBI]
unique_newdata_chr_NCBI <- newdata_chr_NCBI[!duplicated(newdata_chr_NCBI), ]
j_equal_c_NCBI <- subset(unique_newdata_chr_NCBI,class_code.x.1 %in% c("j","=","c"))
u_NCBI <- subset(unique_newdata_chr_NCBI,class_code.x.1 %in% "u")
j_NCBI<- subset(unique_newdata_chr_NCBI,class_code.x.1 %in% "j")
j_u_equal_c_NCBI <- subset(unique_newdata_chr_NCBI,class_code.x.1 %in% c("j","=","u","c"))
#genes for ENSEMBL
keeps_EN <- c("major_iso_id.x", "chr", "class_code.y.1")
newdata_chr_EN <- data_chr[keeps_EN]
unique_newdata_chr_EN <- newdata_chr_EN[!duplicated(newdata_chr_EN), ]
j_equal_c_EN <- subset(unique_newdata_chr_EN,class_code.y.1 %in% c("j","=","c"))
u_EN <- subset(unique_newdata_chr_EN,class_code.y.1 %in% "u")
j_EN <- subset(unique_newdata_chr_EN,class_code.y.1 %in% "j")
j_u_equal_c_EN <- subset(unique_newdata_chr_EN,class_code.y.1 %in% c("j","=","u","c"))
#combine datasets from NCBI and EN so you can remove duplicates
j_equal_both <- data.frame(rbind(j_equal_NCBI,j_equal_EN))
j_equal_c_both <- data.frame(rbind(j_equal_c_NCBI,j_equal_c_EN))
u_both <- data.frame(rbind(u_NCBI,u_EN))
#removing duplicates
j_equal_both[duplicated(j_equal_both),]
j_equal_c_both[duplicated(j_equal_both),]
j_equal_both <- j_equal_both[!duplicated(j_equal_both), ]
u_both <- u_both[!duplicated(u_both), ]
u_both[duplicated(u_both),]

#now for other horse transcriptomes
#genes for NCBI
keeps_Hestand <- c("major_iso_id.x", "chr", "class_code.x")
newdata_chr_Hestand <- data_chr[keeps_Hestand]
unique_newdata_chr_Hestand <- newdata_chr_Hestand[!duplicated(newdata_chr_Hestand), ]
j_equal_c_Hestand <- subset(unique_newdata_chr_Hestand,class_code.x %in% c("j","=","c"))
u_Hestand <- subset(unique_newdata_chr_Hestand,class_code.x %in% "u")
j_Hestand<- subset(unique_newdata_chr_Hestand,class_code.x %in% "j")
j_u_equal_c_Hestand <- subset(unique_newdata_chr_Hestand,class_code.x %in% c("j","=","u","c"))
#genes for ISME
keeps_ISME <- c("major_iso_id.x", "chr", "class_code.y")
newdata_chr_ISME <- data_chr[keeps_ISME]
unique_newdata_chr_ISME <- newdata_chr_ISME[!duplicated(newdata_chr_ISME), ]
j_equal_c_ISME <- subset(unique_newdata_chr_ISME,class_code.y %in% c("j","=","c"))
u_ISME <- subset(unique_newdata_chr_ISME,class_code.y %in% "u")
j_ISME <- subset(unique_newdata_chr_ISME,class_code.y %in% "j")
j_u_equal_c_ISME <- subset(unique_newdata_chr_ISME,class_code.y %in% c("j","=","u","c"))
#combine datasets from NCBI and EN so you can remove duplicates
j_equal_both_horse <- data.frame(rbind(j_equal_Hestand,j_equal_ISME))
j_equal_c_both_horse <- data.frame(rbind(j_equal_c_Hestand,j_equal_c_ISME))
u_both <- data.frame(rbind(u_Hestand,u_ISME))
#removing duplicates
j_equal_both_horse[duplicated(j_equal_both_horse),]
j_equal_c_both_horse[duplicated(j_equal_both_horse),]
j_equal_both_horse <- j_equal_both_horse[!duplicated(j_equal_both_horse), ]
u_both_horse <- u_both_horse[!duplicated(u_both_horse), ]
u_both_horse[duplicated(u_both_horse),]
#merging the classcodes by major_iso_id.x
r <-Reduce(function(x, y) merge(x, y, all=TRUE), list(j_u_equal_c_NCBI,j_u_equal_c_EN,j_u_equal_c_ISME,j_u_equal_c_Hestand))
r_noNA <- r[with(r,class_code.x.1 %in% c("j","=", "u") & class_code.y.1 %in% c("j","=", "u") & class_code.y %in% c("j","=", "u") & class_code.x %in% c("j","=", "u")), ]
library(reshape2)
r_noNA_class <- r_noNA[ ,3:6]
colnames(r_noNA_class)[colnames(r_noNA_class)==c("class_code.x.1","class_code.y.1","class_code.y","class_code.x")] <- c("NCBI","ENSEMBL","ISME","Hestand")
flip<-t(r_noNA_class)
df.melt <- melt(flip, id=rownames(flip))
head(df.melt)
#The bar graph (Venn substitution)
ggplot(df.melt,aes(Var1, fill=value)) + 
  geom_bar(position="dodge") + xlab("database") + ylab("transcript count") +
  guides(fill=guide_legend(title="class code",
                           labels=c("match", "similar","novel")))


#making the genes/chr plot with chr in order
#Making larger dataset for input
keeps_all <- c("major_iso_id.x", "chr", "class_code.x.1", "class_code.y.1")
newdata_chr <- data_chr[keeps_all]


chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX", "chrM")
b<- ggplot(subset(unique_newdata_chr,class_code.x.1 %in% c("j","=", "u"))) + scale_x_discrete(limits = chrs) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + theme(legend.position="top") + ylab("gene count") + scale_fill_discrete(name="class code",labels=c("in NCBI","similar to NCBI entry", "novel"))
b
d<- ggplot(subset(unique_newdata_chr,class_code.x.1 %in% c("j","=", "u"))) + scale_x_discrete(limits = chrs) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + theme(legend.position="top") + ylab("gene count") + 
  scale_fill_discrete(name="class code",labels=c("in NCBI","similar to NCBI entry", "novel"))
d


# Making step line for chr size
datachrN<-read.csv("equCab2.chrom.sizes.csv", header=T)
reordered_chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX", "chrM")
datachrN_new<- data.frame(x=reordered_chr, y=datachrN$size)
#NEED the group=1 here for some reason?!
c<- ggplot(datachrN_new, aes(x,y, group=1)) + scale_x_discrete(limits = chrs) + geom_step(aes(x,y)) + xlab("chr") + ylab("bp size")
c

# With chr trend line
library(Rmisc)
multiplot(b,c,cols=1)
library(gridExtra)
grid.arrange(b, c, ncol=1)
ggplot() + geom_bar(subset(unique_newdata_chr_NCBI,class_code.x.1 %in% c("j","=", "u")),aes(chr, group=class_code.x.1, colour=class_code.x.1)) +
  geom_bar(subset(unique_newdata_chr_EN,class_code.x.1 %in% c("j","=", "u")),aes(chr, group=class_code.x.1, colour=class_code.x.1)) + 
  scale_x_discrete(limits = chrs) + guides(fill=guide_legend(title="Class code"))

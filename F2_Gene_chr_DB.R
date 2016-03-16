setwd("~/Desktop/HorseAnnotation")
library(ggplot2)
###transcripts per chrs
data_chr<- read.csv("Annotated_GeneTable.csv", header=T)

###subsetting, parsing out different class codes for each database
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

#other horse transcriptomes
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

#To try and subset and stack u transcripts in with total transcripts
library(plyr)
library(ggplot2)
ggplot(subset(data_chr,class_code.x.1 %in% c("j","=", "u"))) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + guides(fill=guide_legend(title="class code"))

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
b<- ggplot(subset(unique_newdata_chr_NCBI,class_code.x.1 %in% c("j","=", "u"))) + scale_x_discrete(limits = chrs) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + theme(legend.position="top") + ylab("gene count") + 
  guides(fill=guide_legend(title="class code",
                           labels=c("match", "similar","novel")))
b

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

#getting novel genes
keeps_novel <- c("major_iso_id.x", "chr", "class_code.x.1","class_code.y.1","class_code.x","class_code.y","strand","len.x", "gene.ID")
newdata_novel <- data_chr[keeps_novel]
head(newdata_novel)
#subsetting data for any entries that got class code "u" in any dataset
u_newdata <- subset(newdata_novel,class_code.x.1 %in% "u" | class_code.y.1 %in% "u" | class_code.x %in% "u" | class_code.y %in% "u")
write.csv(u_newdata, file="class_code_u_transcripts.csv")
#subsetting data with class code "u" in every dataset
u <- split(u_newdata, with(u_newdata, interaction(class_code.x.1,class_code.y.1,class_code.x,class_code.y)), drop = TRUE)
u_newdata_limited <- u$u.u.u.u
write.csv(u_newdata_limited, file="novel_transcripts.csv")

#merging the novel genes with their expression in each tissue
tissues <- read.csv("allTissues_isoformTPM.csv", header=T)
m <- merge(u_newdata_limited, tissues, by.x="major_iso_id.x",by.y="isoformName" )
write.csv(m, file="novelisoform_expression.csv")
#merging start/stop information
#visualizing expression of novel isoforms across tissues




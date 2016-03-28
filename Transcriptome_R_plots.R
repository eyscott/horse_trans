library(ggplot2)
setwd("~/Desktop/HorseAnnotation")
####Figure 2: general transcriptome stats
###comparing databases figure
#upload meta-data cuffcompare file between all databases
data_chr<- read.csv("Annotated_GeneTable.csv", header=T)

##subsetting/parsing out different class codes (cuffcompare) for each database
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

#merging the classcodes by "major_iso_id.x" id, with each column representing a db...then "metling" data so it can be plotted
r <-Reduce(function(x, y) merge(x, y, all=TRUE), list(j_u_equal_c_NCBI,j_u_equal_c_EN,j_u_equal_c_ISME,j_u_equal_c_Hestand))
r_noNA <- r[with(r,class_code.x.1 %in% c("j","=", "u") & class_code.y.1 %in% c("j","=", "u") & class_code.y %in% c("j","=", "u") & class_code.x %in% c("j","=", "u")), ]
library(reshape2)
r_noNA_class <- r_noNA[ ,3:6]
colnames(r_noNA_class)[colnames(r_noNA_class)==c("class_code.x.1","class_code.y.1","class_code.y","class_code.x")] <- c("NCBI","ENSEMBL","ISME","Hestand")
flip<-t(r_noNA_class)
df.melt <- melt(flip, id=rownames(flip))
head(df.melt)
#The bar graph comparing dbs (Venn substitution)
ggplot(df.melt,aes(Var1, fill=value)) + 
  geom_bar(position="dodge") + xlab("database") + ylab("transcript count") +
  scale_fill_manual(values=c("tomato","green3","dodgerblue1"),
                    name="class code",
                    breaks=c("=", "j", "u"),
                    labels=c("match", "similar","novel"))



##making the genes/chr/class code plot with chr in order
#making input dataset for class code's distribution over chrs plot
keeps_all <- c("major_iso_id.x", "chr", "class_code.x.1", "class_code.y.1")
newdata_chr <- data_chr[keeps_all]
#Making the plot just against NCBI class codes
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX", "chrM")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X", "M")
c<- ggplot(subset(unique_newdata_chr_NCBI,class_code.x.1 %in% c("j","=", "u"))) + scale_x_discrete(limits = chrs, labels = chrs_N) + 
  geom_bar(aes(chr, group=class_code.x.1, colour=class_code.x.1)) + theme(legend.position="top") + ylab("gene count") + 
  scale_colour_discrete(name  ="class code",
                        labels=c("match","similar","novel"),
                        expand=2) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))
c

# Making step line for chr size
datachrN<-read.csv("equCab2.chrom.sizes.csv", header=T)
reordered_chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX", "chrM")
datachrN_new<- data.frame(x=reordered_chr, y=datachrN$size)
#NEED the group=1 here for some reason?!
d<- ggplot(datachrN_new, aes(x,y, group=1)) + scale_x_discrete(limits = chrs, labels = chrs_N) + geom_step(aes(x,y)) + xlab("chr") + ylab("bp size") +
  theme(axis.text = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))
d
# Making a figure with class_code/chr plot and chr size trend line
library(Rmisc)
multiplot(c,d,cols=1)

### Extracting information about the novel genes
keeps_novel <- c("major_iso_id.x", "chr", "class_code.x.1","class_code.y.1","class_code.x","class_code.y","strand","len.x", "gene.ID")
newdata_novel <- data_chr[keeps_novel]
head(newdata_novel)
#subsetting novel genes data for any entries that have a class code "u" in any dataset
u_newdata <- subset(newdata_novel,class_code.x.1 %in% "u" | class_code.y.1 %in% "u" | class_code.x %in% "u" | class_code.y %in% "u")
write.csv(u_newdata, file="class_code_u_transcripts.csv")
#subsetting data with class code "u" in every dataset
u <- split(u_newdata, with(u_newdata, interaction(class_code.x.1,class_code.y.1,class_code.x,class_code.y)), drop = TRUE)
u_newdata_limited <- u$u.u.u.u
write.csv(u_newdata_limited, file="novel_transcripts.csv")

##merging the novel genes with their expression (in TPM) in each tissue
#start by uploading TPM values for all isoforms
tissues <- read.csv("allTissues_isoformTPM.csv", header=T)
#merge the TPM values with the novel gene by "major_iso_id.x"
m <- merge(u_newdata_limited, tissues, by.x="major_iso_id.x",by.y="isoformName" )
write.csv(m, file="novelisoform_expression.csv")
#Add exon # information to this table too by merge function
exon <- read.csv("all_tissues_frac0.05.nonGuided_Cufflinks.bed.titles.csv", header=T)
m2 <- merge(m,exon, by.x="major_iso_id.x",by.y="name")
write.csv(m2, file="novelisoform_exonNum.csv")

##visualizing expression of novel isoforms across tissues with exon number information
data3 <- read.csv("novelisoform_exonNum.csv", header=T)
desired <- c("major_iso_id.x","BrainStem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord","exon_number")
data3 <- data3[desired]
#reshaping data for geom_bar
df.melt <- melt(data3, id=c("major_iso_id.x","exon_number"))
## Plot sum of TPM of novel genes with exon number subset per tissue
#initial plot with no cap on exon # shown
ggplot(subset(df.melt, exon_number %in% c(1:10)), aes(variable,group=exon_number,fill=exon_number)) + 
  geom_bar(aes(weight=value),position="stack") + xlab("tissue") + 
  ylab("TPM") + scale_fill_gradientn(colours=rainbow(10)) +  
  guides(fill=guide_legend(title="exon number"))
#with up to 10 exons (but proportion of anything above 5 exons is not visible)
ggplot(subset(df.melt, exon_number %in% c(1:10)), aes(variable,group=exon_number,fill=exon_number)) + 
  geom_bar(aes(weight=value),position="stack") + xlab("tissue") + 
  ylab("sum(TPM)") + scale_fill_gradientn(colours=rainbow(10),
                                          breaks=c(1,2,3,4,5,6,7,8,9,10)) +
  guides(fill=guide_legend(title="number of exons"))
#with up to 5 exons
ggplot(subset(df.melt, exon_number %in% c(1:5)), aes(variable,group=exon_number,fill=exon_number)) + 
  geom_bar(aes(weight=value),position="stack") + xlab("tissue") + 
  ylab("sum(TPM)") + scale_fill_gradientn(colours=rainbow(10),
                                          breaks=c(1,2,3,4,5)) +
  guides(fill=guide_legend(title="number of exons")) +
  theme(axis.text = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))



####Figure 4 Tissue specificty figure
#read in dataset with expression values for each term and a dataset with annotation for each gene
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
#plotting nuclear vs mitochondrial origin transcripts per tissue
ggplot(bar, aes(x=Var1, y=value, fill=factor(Var2)))+
  geom_bar(position="stack", stat="identity") + 
  ylab("proportion of transcriptional output") + xlab("tissue") +
  guides(fill=guide_legend(title="gene origin",
                           labels=c("nuclear", "mt")))

##making the PCA/MDS plot for tissue expression
data<-read.csv("allTissues_geneTPM.csv", header=T)
##PCA for transcriptome data
#getting basic PCA statistics
data_PCA <- prcomp(data[-1,-1], scale = T)
plot(prcomp(data[-1,-1], scale = T))
summary(prcomp(data[-1,-1], scale = T))
biplot(prcomp(data[-1,-1], scale = T))
#PCA- new flavour with factoextra program...visuals are one step above basic
library("factoextra")
eig.val <- get_eigenvalue(data_PCA)
head(eig.val)
# Eigenvalues
eig <- (data_PCA$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.tissue <- data.frame(eig = eig, variance = variance,
                         cumvariance = cumvar)
head(eig.tissue)
#plot showing variance captured
barplot(eig.tissue[, 2], names.arg=1:nrow(eig.tissue), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
var <- get_pca_var(data_PCA)
var
var$coord[, 1:4]
#2D MDS plot based on coord values
coord <- data.frame(var$coord[ ,1:4])
pc1.1 <- qplot(x=Dim.1, y=Dim.2, data=coord, colour=factor(rownames(rotation))) +
  theme(legend.position="right")
pc1.1
##3D MDS plot of pc1:pc3
#need to make separate column with these tisue type so they can be variables with levels
tissue <-c("Brainstem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord")
coord_1 <- cbind(coord,tissue)
#setting up colour palatte
colors<- c("cadetblue","turquoise1", "yellow4","yellowgreen","chocolate2","pink2","peachpuff2","blueviolet")
colors <- colors[as.numeric(coord_1$tissue)]
#setting shape palatte
shapes <- c(1,2,3,4,5,6,7,8)
shapes <-shapes[as.numeric(coord_1$tissue)]
#3d plot with colour labels and legend
scatterplot3d(coord[,1:3],
              xlab = "Dim.1(32%)",
              ylab = "",
              zlab = "Dim.3(13%)",
              color=colors,
              angle = 18,box=FALSE,pch=16,scale.y=0.8, mar=c(2.8,2.8,0.2,0.05))
legend("topright",legend = levels(coord_1$tissue),
       col =  c("cadetblue","turquoise1", "yellow4","yellowgreen","chocolate2","pink2","peachpuff2","blueviolet"),
       pch=16,xpd = T,horiz = F)
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])
y <- dims[3]+ 0.08*diff(dims[3:4])
text(x,y,"Dim.2(23%)",srt=20)
##above is the money plot with adjusted margins, rotated y-axis label and positioned legend!
#3d plot with shapes and legend
scatterplot3d(coord[,1:3],
              xlab = "Dim.1(32%)",
              ylab = "Dim.2(23%)",
              zlab = "Dim.3(13%)",
              angle = 55,box=FALSE,pch=shapes)
legend("right", legend = levels(coord_1$tissue), pch=shapes)

####Figure 5: Isoform plots
###Making genes vs isoforms scatter
#matching colours to 3D MDS plot plot above
data<-read.csv("HT_isoformsvsGenes.csv", header=T)
colors<- c("cadetblue","turquoise1", "yellow4","yellowgreen","chocolate2","pink2","peachpuff2","blueviolet")
data$label<-c("Brainstem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord")
names(colors)<-levels(data$label)
colScale <- scale_colour_manual(name = "tissue",values = colors)
##The genes vs isoforms plot
g<- ggplot(data,aes(Number_genes,Number_isoforms)) + geom_point(aes(color=label)) +
  xlab("Number of genes") + ylab("Number of Isoforms") + 
  scale_color_discrete(name  =NULL,
                       labels=c("Brainstem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord"))
g + colScale

###Making genes vs isoforms histo with unique absent and present genes per tissue
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
df.melt <- melt(newdata4, id="Tissue")
bar <- group_by(df.melt, variable, Tissue)
#Can see absent_Uisoforms from above code, but NUisoforms is too small
# Making new unique isoforms into heatmap format for better visualization and with absent isoform proportion underneath
ggplot() +
geom_bar(data=data2, aes(x=Tissue,y=X.isoforms, fill=X.Uisoforms), stat="identity") +
geom_bar(data=Absent_Uisoforms_data, aes(x=Tissue,y=Absent_Uisoforms),stat = "identity") + 
ylab("Number of isoforms") + guides(fill=guide_legend(title="Number of Unique Isoforms"))
            
#To have just unique present and absent isoforms (if you don't really care about total # isoforms/tissue)
ggplot() +
geom_bar(data=NUisoforms_data, aes(x=Tissue,y=NUisoforms,color="aliceblue",label=NUisoforms), stat="identity") +
geom_bar(data=Absent_Uisoforms_data, aes(x=Tissue,y=Absent_Uisoforms,color="red",label=Absent_Uisoforms),stat = "identity") + 
ylab("Number of isoforms") + scale_color_discrete(name="Unique isoforms",
                                                               labels=c("present","absent"))
           
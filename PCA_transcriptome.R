#Setting up working environment
setwd("~/Desktop/HorseAnnotation")
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
##do 3d plot of pc1:pc3
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
              angle = 22,box=FALSE,pch=16,scale.y=0.8, mar=c(2.8,2.8,0.2,0.05))
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


#PCA plots based on tissue type
fviz_pca_var(data_PCA)
fviz_pca_var(data_PCA, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=55) + theme_minimal()
ind <- get_pca_ind(data_PCA)
ind
fviz_pca_ind(data_PCA)
#PCA plot  with variables and individuals
fviz_pca_biplot(data_PCA,  geom = "text") +
  theme_minimal()


fviz_pca_biplot(data_PCA, label ="var")

# Control automatically the color of individuals using the cos2
fviz_pca_biplot(data_PCA, label ="var", col.ind="cos2") +
  theme_minimal()
# Color individuals by groups
b <-as.factor(c("BrainStem","Cerebellum","SpinalCord","Retina","Skin","Muscle", "Embryo.ICM","Embryo.TE"))
fviz_pca_ind(data_PCA, label="none", habillage=b)


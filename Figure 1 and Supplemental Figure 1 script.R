###########################################################

# PCA of function

###########################################################

pc<-prcomp(dat.rat.na.log[,c(5:17)], scale=TRUE)

Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

dat.rat.na.log$Color<-dat.rat.na.log$Trophic
dat.rat.na.log$Color[dat.rat.na.log$Color == "biter"]<-"black"
dat.rat.na.log$Color[dat.rat.na.log$Color == "suction"]<-"red"


###########################
#   Plot Figure 1
###########################

pdf(file="PCA_function_labels.pdf", width=20, height=20)
plot(pc$x[,1], pc$x[,2], pch=16, col=dat.rat.na$Color,
     xlab="PC1", ylab="PC2", main="PC1 vs. PC2")
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Trophic == "biter"], ycoord= pc$x[,2][dat.rat.na.log$Trophic == "biter"], lcolor="black" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Trophic == "suction"], ycoord= pc$x[,2][dat.rat.na.log$Trophic == "suction"], lcolor="red" )
dev.off()

###########################
#   Plot Supplemental Figure 1
###########################

pdf(file="PCA_function_no labels.pdf", width=20, height=20)
plot(pc$x[,1], pc$x[,2], pch=16, col=dat.rat.na$Color,
     xlab="PC1", ylab="PC2", main="PC1 vs. PC2")
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Trophic == "biter"], ycoord= pc$x[,2][dat.rat.na.log$Trophic == "biter"], lcolor="black" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Trophic == "suction"], ycoord= pc$x[,2][dat.rat.na.log$Trophic == "suction"], lcolor="red" )
text(pc$x[,1], pc$x[,2], labels=dat.rat.na.log$Species)
dev.off()
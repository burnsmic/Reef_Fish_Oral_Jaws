###########################################
#     PCA of Diet
###########################################


pc<-prcomp(dat.rat.na.log[,c(5:17)], scale=TRUE)


dat.rat.na.log$Color<-dat.rat.na.log$Diet
dat.rat.na.log$Color[dat.rat.na.log$Color == "GC"]<-"#F9A21E"
dat.rat.na.log$Color[dat.rat.na.log$Color == "OM"]<-"#F06A95"
dat.rat.na.log$Color[dat.rat.na.log$Color == "HD"]<-"#11619B"
dat.rat.na.log$Color[dat.rat.na.log$Color == "MI"]<-"#27B79B"
dat.rat.na.log$Color[dat.rat.na.log$Color == "SI"]<-"#47ADE1"
dat.rat.na.log$Color[dat.rat.na.log$Color == "PK"]<-"#9B519E" 
dat.rat.na.log$Color[dat.rat.na.log$Color == "D"]<-"#E21C4D" 

gc.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#F9A21E"])
gc.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#F9A21E"])
om.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#F06A95"])
om.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#F06A95"])
hd.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#11619B"])
hd.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#11619B"])
mis.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#27B79B"])
mis.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#27B79B"])
pk.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#9B519E"])
pk.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#9B519E"])
si.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#47ADE1"])
si.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#47ADE1"])
mib.x<-mean(pc$x[,1][dat.rat.na.log$Color == "#E21C4D"])
mib.y<-mean(pc$x[,2][dat.rat.na.log$Color == "#E21C4D"])


########################################
#          Plot Figure 2
########################################

#pdf(file="PCA_PC1 vs PC2_diet_species labels.pdf", width=20, height=20)        
plot(pc$x[,1], pc$x[,2], col=dat.rat.na.log$Color, pch=16)
points(gc.x,gc.y, pch = 1, col="#F9A21E")
points(om.x,om.y, pch = 1, col="#F06A95")
points(hd.x,hd.y, pch = 1, col="#11619B")
points(mis.x,mis.y, pch = 1, col="#27B79B")
points(si.x,si.y, pch = 1, col="#47ADE1")
points(pk.x,pk.y, pch = 1, col="#9B519E")
points(mib.x,mib.y, pch = 1, col="#E21C4D")
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#F9A21E"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#F9A21E"], lcolor="#F9A21E" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#F06A95"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#F06A95"], lcolor="#F06A95" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#11619B"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#11619B"], lcolor="#11619B" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#27B79B"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#27B79B"], lcolor="#27B79B")
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#47ADE1"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#47ADE1"], lcolor="#47ADE1" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#9B519E"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#9B519E"], lcolor="#9B519E" )
Plot_ConvexHull(xcoord=pc$x[,1][dat.rat.na.log$Color == "#E21C4D"], ycoord= pc$x[,2][dat.rat.na.log$Color == "#E21C4D"], lcolor="#E21C4D" )
#text(pc$x[,1], pc$x[,2], labels=dat.rat.na$Species)
#dev.off()

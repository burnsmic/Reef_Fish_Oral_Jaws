library(evobiR)
dd<-dat.rat.na.log

######################################################
# Reorder data to match the order of the phylogeny
######################################################

rownames(dd)<-dat.rat.na.log$Species

dat.rat.na.order<-ReorderData(tree.final, dd, taxa.names="row names")

#######################################################
# Make a vector of each trait
#######################################################

gape<- dat.rat.na.order$log.Gape.head.length
names(gape)<-dat.rat.na.order$Species

Jaw.Protrusion<-dat.rat.na.order$log.prod.head.length
names(Jaw.Protrusion)<-dat.rat.na.order$Species

maxilla.length<-dat.rat.na.order$log.max.head.length
names(maxilla.length)<-dat.rat.na.order$Species

premax.length<-dat.rat.na.order$log.premax.head.length
names(premax.length)<-dat.rat.na.order$Species

lower.jaw<-dat.rat.na.order$log.low.jaw.lever
names(lower.jaw)<-dat.rat.na.order$Species

close.in<-dat.rat.na.order$log.close.out
names(close.in)<-dat.rat.na.order$Species

open.in<-dat.rat.na.order$log.open.out
names(open.in)<-dat.rat.na.order$Species

adductor.mass<-dat.rat.na.order$log.add.mass.body.mass
names(adductor.mass)<-dat.rat.na.order$Species

head.height<-dat.rat.na.order$log.head.height.head.length
names(head.height)<-dat.rat.na.order$Species

x.pos.jaw<-dat.rat.na.order$log.x.low.jaw.head.length
names(x.pos.jaw)<-dat.rat.na.order$Species

y.pos.jaw<-dat.rat.na.order$log.y.low.jaw.head.length
names(y.pos.jaw)<-dat.rat.na.order$Species

x.pal<-dat.rat.na.order$log.x.palatine.head.length
names(x.pal)<-dat.rat.na.order$Species

y.pal<-dat.rat.na.order$log.y.palatine.head.length
names(y.pal)<-dat.rat.na.order$Species


############################################################################
#  Create a dataframe where empirical correlation results will be stored
############################################################################

corr_results <- setNames(data.frame(matrix(ncol = 13, nrow = 13)), c("gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", 
                                                                     "Premax Length", "Maxilla Length", "Lower Jaw Length", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                                                                     "x.palatine", "y.palatine"))
rownames(corr_results)<-colnames(corr_results)

# Making diagonal values equal to 1

corr_results[1,1]<-1
corr_results[2,2]<-1
corr_results[3,3]<-1
corr_results[4,4]<-1
corr_results[5,5]<-1
corr_results[6,6]<-1
corr_results[7,7]<-1
corr_results[8,8]<-1
corr_results[9,9]<-1
corr_results[10,10]<-1
corr_results[11,11]<-1
corr_results[12,12]<-1
corr_results[13,13]<-1

#############################################################
#  Calculating correlations and storing the values in a data frame
#############################################################


#################################
#      Gape Vs Jaw Protrusion
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(Jaw.Protrusion, tree.final)

picModel.gape.pro<-cor.test(hPic,aPic,method="pearson")


corr_results[1,2]<-picModel.gape.pro$estimate
corr_results[2,1]<-picModel.gape.pro$estimate

#################################
#      Gape Vs adductor 
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.gape.add <- cor.test(hPic,aPic,method="pearson")

corr_results[1,3]<-picModel.gape.add$estimate
corr_results[3,1]<-picModel.gape.add$estimate


#################################
#      Gape Vs head
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(head.height, tree.final)
picModel.gape.head <- cor.test(hPic,aPic,method="pearson")

corr_results[1,4]<-picModel.gape.head$estimate
corr_results[4,1]<-picModel.gape.head$estimate

#################################
#      Gape Vs PreMaxilla
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(premax.length, tree.final)
picModel.gape.premax <- cor.test(hPic,aPic,method="pearson")

corr_results[1,5]<-picModel.gape.premax$estimate
corr_results[5,1]<-picModel.gape.premax$estimate


#################################
#      Gape Vs Maxilla
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(maxilla.length, tree.final)
picModel.gape.max <- cor.test(hPic,aPic,method="pearson")

corr_results[1,6]<-picModel.gape.max$estimate
corr_results[6,1]<-picModel.gape.max$estimate


#################################
#      Gape Vs Lower Jaw Length
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(lower.jaw, tree.final)
picModel.gape.low <- cor.test(hPic,aPic,method="pearson")

corr_results[1,7]<-picModel.gape.low$estimate
corr_results[7,1]<-picModel.gape.low$estimate

#################################
#      Gape Vs closing ma
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(close.in, tree.final)
picModel.gape.close <- cor.test(hPic,aPic,method="pearson")

corr_results[1,8]<-picModel.gape.close$estimate
corr_results[8,1]<-picModel.gape.close$estimate

#################################
#      Gape Vs opening ma
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(open.in, tree.final)
picModel.gape.open <- cor.test(hPic,aPic,method="pearson")


corr_results[1,9]<-picModel.gape.open$estimate
corr_results[9,1]<-picModel.gape.open$estimate


#################################
#      Gape Vs X pos lower jaw
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(x.pos.jaw, tree.final)
picModel.gape.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[1,10]<-picModel.gape.x.low$estimate
corr_results[10,1]<-picModel.gape.x.low$estimate

#################################
#      Gape Vs Y pos lower jaw
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.gape.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[1,11]<-picModel.gape.y.low$estimate
corr_results[11,1]<-picModel.gape.y.low$estimate

#################################
#      Gape Vs X pos pal
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.gape.x.pal <- cor.test(hPic,aPic,method="pearson")


corr_results[1,12]<-picModel.gape.x.pal$estimate
corr_results[12,1]<-picModel.gape.x.pal$estimate

#################################
#      Gape Vs Y pos pal
##################################
hPic <- pic(gape, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.gape.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[1,13]<-picModel.gape.y.pal$estimate
corr_results[13,1]<-picModel.gape.y.pal$estimate

#################################
#      Jaw Protrusion Vs Adductor mass
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.pro.add <- cor.test(hPic,aPic,method="pearson")

corr_results[2,3]<-picModel.pro.add$estimate
corr_results[3,2]<-picModel.pro.add$estimate


#################################
#      Jaw Protrusion Vs Head
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(head.height, tree.final)
picModel.pro.head <- cor.test(hPic,aPic,method="pearson")

corr_results[2,4]<-picModel.pro.head$estimate
corr_results[4,2]<-picModel.pro.head$estimate


#################################
#      Premaxilla Vs Jaw Protrusion
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(Jaw.Protrusion, tree.final)
picModel.pro.pre <- cor.test(hPic,aPic,method="pearson")

corr_results[2,5]<-picModel.pro.pre$estimate
corr_results[5,2]<-picModel.pro.pre$estimate


#################################
#      Maxilla Vs Jaw protrusion
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(Jaw.Protrusion, tree.final)
picModel.pro.max <- cor.test(hPic,aPic,method="pearson")

corr_results[2,6]<-picModel.pro.max$estimate
corr_results[6,2]<-picModel.pro.max$estimate

#################################
#      Lower Jaw Vs Jaw protrusion
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(Jaw.Protrusion, tree.final)
picModel.pro.low <- cor.test(hPic,aPic,method="pearson")

corr_results[2,7]<-picModel.pro.low$estimate
corr_results[7,2]<-picModel.pro.low$estimate

#################################
#      Jaw Protrusion Vs Close ma
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(close.in, tree.final)
picModel.pro.close <- cor.test(hPic,aPic,method="pearson")

corr_results[2,8]<-picModel.pro.close$estimate
corr_results[8,2]<-picModel.pro.close$estimate


#################################
#      Jaw Protrusion Vs open ma
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(open.in, tree.final)
picModel.pro.open <- cor.test(hPic,aPic,method="pearson")

corr_results[2,9]<-picModel.pro.open$estimate
corr_results[9,2]<-picModel.pro.open$estimate

#################################
#      Jaw Protrusion Vs X position lower jaw
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(x.pos.jaw, tree.final)
picModel.pro.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[2,10]<-picModel.pro.x.low$estimate
corr_results[10,2]<-picModel.pro.x.low$estimate


#################################
#      Jaw Protrusion Vs Y position lower jaw
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.pro.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[2,11]<-picModel.pro.y.low$estimate
corr_results[11,2]<-picModel.pro.y.low$estimate


#################################
#      Jaw Protrusion Vs X position palatine
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.pro.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[2,12]<-picModel.pro.x.pal$estimate
corr_results[12,2]<-picModel.pro.x.pal$estimate


#################################
#      Jaw Protrusion Vs Y position palatine
##################################
hPic <- pic(Jaw.Protrusion, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.pro.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[2,13]<-picModel.pro.y.pal$estimate
corr_results[13,2]<-picModel.pro.y.pal$estimate


#################################
#      Adductor Vs head
##################################
hPic <- pic(adductor.mass, tree.final)
aPic <- pic(head.height, tree.final)
picModel.add.head <- cor.test(hPic,aPic,method="pearson")

corr_results[3,4]<-picModel.add.head$estimate
corr_results[4,3]<-picModel.add.head$estimate


#################################
#      Premaxilla Vs adductor mass
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.pre <- cor.test(hPic,aPic,method="pearson")

corr_results[3,5]<-picModel.add.pre$estimate
corr_results[5,3]<-picModel.add.pre$estimate

#################################
#      Maxilla Vs adductor mass
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.max <- cor.test(hPic,aPic,method="pearson")

corr_results[3,6]<-picModel.add.max$estimate
corr_results[6,3]<-picModel.add.max$estimate

#################################
#      Lower Jaw Vs adductor mass
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.low <- cor.test(hPic,aPic,method="pearson")

corr_results[3,7]<-picModel.add.low$estimate
corr_results[7,3]<-picModel.add.low$estimate

#################################
#      Closing Vs adductor
##################################
hPic <- pic(close.in, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.close <- cor.test(hPic,aPic,method="pearson")

corr_results[3,8]<-picModel.add.close$estimate
corr_results[8,3]<-picModel.add.close$estimate

#################################
#      Adductor Vs open
##################################
hPic <- pic(adductor.mass, tree.final)
aPic <- pic(open.in, tree.final)
picModel.add.open <- cor.test(hPic,aPic,method="pearson")

corr_results[3,9]<-picModel.add.open$estimate
corr_results[9,3]<-picModel.add.open$estimate

#################################
#      X position lower jaw Vs adductor mass
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[3,10]<-picModel.add.x.low$estimate
corr_results[10,3]<-picModel.add.x.low$estimate

#################################
#      Y position lower jaw Vs adductor
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[3,11]<-picModel.add.y.low$estimate
corr_results[11,3]<-picModel.add.y.low$estimate



#################################
#      X position palatine Vs adductor
##################################
hPic <- pic(x.pal, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[3,12]<-picModel.add.x.pal$estimate
corr_results[12,3]<-picModel.add.x.pal$estimate


#################################
#      Y position palatine Vs adductor 
##################################
hPic <- pic(y.pal, tree.final)
aPic <- pic(adductor.mass, tree.final)
picModel.add.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[3,13]<-picModel.add.y.pal$estimate
corr_results[13,3]<-picModel.add.y.pal$estimate






#################################
#      Premaxilla Vs head
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(head.height, tree.final)
picModel.pre.head <- cor.test(hPic,aPic,method="pearson")

corr_results[4,5]<-picModel.pre.head$estimate
corr_results[5,4]<-picModel.pre.head$estimate

#################################
#      Maxilla Vs head
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.max <- cor.test(hPic,aPic,method="pearson")

corr_results[4,6]<-picModel.head.max$estimate
corr_results[6,4]<-picModel.head.max$estimate

#################################
#      Lower Jaw Vs head
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.low <- cor.test(hPic,aPic,method="pearson")

corr_results[4,7]<-picModel.head.low$estimate
corr_results[7,4]<-picModel.head.low$estimate

#################################
#      Closing Vs head
##################################
hPic <- pic(close.in, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.close <- cor.test(hPic,aPic,method="pearson")

corr_results[4,8]<-picModel.head.close$estimate
corr_results[8,4]<-picModel.head.close$estimate


#################################
#      head Vs open
##################################
hPic <- pic(head.height, tree.final)
aPic <- pic(open.in, tree.final)
picModel.head.open <- cor.test(hPic,aPic,method="pearson")

corr_results[4,9]<-picModel.head.open$estimate
corr_results[9,4]<-picModel.head.open$estimate


#################################
#      X position lower jaw Vs head
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[4,10]<-picModel.head.x.low$estimate
corr_results[10,4]<-picModel.head.x.low$estimate

#################################
#      Y position lower jaw Vs head
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[4,11]<-picModel.head.y.low$estimate
corr_results[11,4]<-picModel.head.y.low$estimate



#################################
#      Y position palatine Vs head
##################################
hPic <- pic(y.pal, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[4,12]<-picModel.head.x.pal$estimate
corr_results[12,4]<-picModel.head.x.pal$estimate



#################################
#      X position palatine Vs head
##################################
hPic <- pic(x.pal, tree.final)
aPic <- pic(head.height, tree.final)
picModel.head.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[4,13]<-picModel.head.y.pal$estimate
corr_results[13,4]<-picModel.head.y.pal$estimate

#################################
#      Premaxilla Vs Maxilla
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(maxilla.length, tree.final)
picModel.pre.max <- cor.test(hPic,aPic,method="pearson")

cor.test(hPic,aPic,method="pearson")
cor(hPic,aPic)

corr_results[5,6]<-picModel.pre.max$estimate
corr_results[6,5]<-picModel.pre.max$estimate

#################################
#      Premaxilla Vs Lower Jaw length
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(lower.jaw, tree.final)
picModel.pre.low <- cor.test(hPic,aPic,method="pearson")

corr_results[5,7]<-picModel.pre.low$estimate
corr_results[7,5]<-picModel.pre.low$estimate

#################################
#      Premaxilla Vs closing ma
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(close.in, tree.final)
picModel.pre.close <- cor.test(hPic,aPic,method="pearson")

corr_results[5,8]<-picModel.pre.close$estimate
corr_results[8,5]<-picModel.pre.close$estimate


#################################
#      Premaxilla Vs open ma
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(open.in, tree.final)
picModel.pre.open <- cor.test(hPic,aPic,method="pearson")

corr_results[5,9]<-picModel.pre.open$estimate
corr_results[9,5]<-picModel.pre.open$estimate

#################################
#      Premaxilla Vs X position lower jaw
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(x.pos.jaw, tree.final)
picModel.pre.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[5,10]<-picModel.pre.x.low$estimate
corr_results[10,5]<-picModel.pre.x.low$estimate

#################################
#      Premaxilla Vs Y position lower jaw
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.pre.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[5,11]<-picModel.pre.y.low$estimate
corr_results[11,5]<-picModel.pre.y.low$estimate

#################################
#      Premaxilla Vs X position palatine
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.pre.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[5,12]<-picModel.pre.x.pal$estimate
corr_results[12,5]<-picModel.pre.x.pal$estimate

#################################
#      Premaxilla Vs y position palatine
##################################
hPic <- pic(premax.length, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.pre.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[5,13]<-picModel.pre.y.pal$estimate
corr_results[13,5]<-picModel.pre.y.pal$estimate

#################################
#      Maxilla Vs Lower JAw
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(lower.jaw, tree.final)
picModel.max.low <- cor.test(hPic,aPic,method="pearson")

corr_results[6,7]<-picModel.max.low$estimate
corr_results[7,6]<-picModel.max.low$estimate

#################################
#      Maxilla Vs close in
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(close.in, tree.final)
picModel.max.close <- cor.test(hPic,aPic,method="pearson")

corr_results[6,8]<-picModel.max.close$estimate
corr_results[8,6]<-picModel.max.close$estimate


#################################
#      Maxilla Vs open in
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(open.in, tree.final)
picModel.max.open <- cor.test(hPic,aPic,method="pearson")

corr_results[6,9]<-picModel.max.open$estimate
corr_results[9,6]<-picModel.max.open$estimate

#################################
#      Maxilla Vs X position lower jaw
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(x.pos.jaw, tree.final)
picModel.max.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[6,10]<-picModel.max.x.low$estimate
corr_results[10,6]<-picModel.max.x.low$estimate


#################################
#      Maxilla Vs Y position lower jaw
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.max.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[6,11]<-picModel.max.y.low$estimate
corr_results[11,6]<-picModel.max.y.low$estimate


#################################
#      Maxilla Vs x palatine
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.max.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[6,12]<-picModel.max.x.pal$estimate
corr_results[12,6]<-picModel.max.x.pal$estimate


#################################
#      Maxilla Vs y palatine
##################################
hPic <- pic(maxilla.length, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.max.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[6,13]<-picModel.max.y.pal$estimate
corr_results[13,6]<-picModel.max.y.pal$estimate

#################################
#      Lower Jaw Vs closing in lever
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(close.in, tree.final)
picModel.low.close <- cor.test(hPic,aPic,method="pearson")

corr_results[7,8]<-picModel.low.close$estimate
corr_results[8,7]<-picModel.low.close$estimate


#################################
#      Lower Jaw Vs open in
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(open.in, tree.final)
picModel.low.open <- cor.test(hPic,aPic,method="pearson")

corr_results[7,9]<-picModel.low.open$estimate
corr_results[9,7]<-picModel.low.open$estimate

#################################
#      Lower Jaw Vs X position lower jaw
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(x.pos.jaw, tree.final)
picModel.low.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[7,10]<-picModel.low.x.low$estimate
corr_results[10,7]<-picModel.low.x.low$estimate


#################################
#      Lower Jaw Vs y position lower jaw
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.low.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[7,11]<-picModel.low.y.low$estimate
corr_results[11,7]<-picModel.low.y.low$estimate


#################################
#      Lower Jaw Vs x position palatine
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.low.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[7,12]<-picModel.low.x.pal$estimate
corr_results[12,7]<-picModel.low.x.pal$estimate


#################################
#      Lower Jaw Vs y position palatine
##################################
hPic <- pic(lower.jaw, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.low.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[7,13]<-picModel.low.y.pal$estimate
corr_results[13,7]<-picModel.low.y.pal$estimate


#################################
#      Closing Vs open
##################################
hPic <- pic(close.in, tree.final)
aPic <- pic(open.in, tree.final)
picModel.close.open <- cor.test(hPic,aPic,method="pearson")

corr_results[8,9]<-picModel.close.open$estimate
corr_results[9,8]<-picModel.close.open$estimate


#################################
#      X position lower jaw Vs closing in
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(close.in, tree.final)
picModel.close.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[8,10]<-picModel.close.x.low$estimate
corr_results[10,8]<-picModel.close.x.low$estimate


#################################
#      Y position lower jaw Vs close in
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(close.in, tree.final)
picModel.close.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[8,11]<-picModel.close.y.low$estimate
corr_results[11,8]<-picModel.close.y.low$estimate

#################################
#      X position palatine Vs close in
##################################
hPic <- pic(x.pal, tree.final)
aPic <- pic(close.in, tree.final)
picModel.close.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[8,12]<-picModel.close.x.pal$estimate
corr_results[12,8]<-picModel.close.x.pal$estimate


#################################
#      Y position palatine Vs close 
##################################
hPic <- pic(y.pal, tree.final)
aPic <- pic(close.in, tree.final)
picModel.close.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[8,13]<-picModel.close.y.pal$estimate
corr_results[13,8]<-picModel.close.y.pal$estimate

#################################
#      X position lower jaw Vs open in
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(open.in, tree.final)
picModel.open.x.low <- cor.test(hPic,aPic,method="pearson")

corr_results[9,10]<-picModel.open.x.low$estimate
corr_results[10,9]<-picModel.open.x.low$estimate


#################################
#      Y position lower jaw Vs open
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(open.in, tree.final)
picModel.open.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[9,11]<-picModel.open.y.low$estimate
corr_results[11,9]<-picModel.open.y.low$estimate

#################################
#      X position palatine Vs open 
##################################
hPic <- pic(x.pal, tree.final)
aPic <- pic(open.in, tree.final)
picModel.open.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[9,12]<-picModel.open.x.pal$estimate
corr_results[12,9]<-picModel.open.x.pal$estimate



#################################
#      Y position palatine Vs open
##################################
hPic <- pic(y.pal, tree.final)
aPic <- pic(open.in, tree.final)
picModel.open.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[9,13]<-picModel.open.y.pal$estimate
corr_results[13,9]<-picModel.open.y.pal$estimate


#################################
#      X position lower jaw Vs Y position lower jaw
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(y.pos.jaw, tree.final)
picModel.x.low.y.low <- cor.test(hPic,aPic,method="pearson")

corr_results[10,11]<-picModel.x.low.y.low$estimate
corr_results[11,10]<-picModel.x.low.y.low$estimate


#################################
#      X position lower jaw Vs X position palatine
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.x.low.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[10,12]<-picModel.x.low.x.pal$estimate
corr_results[12,10]<-picModel.x.low.x.pal$estimate


#################################
#      X position lower jaw Vs y position palatine
##################################
hPic <- pic(x.pos.jaw, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.x.low.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[10,13]<-picModel.x.low.y.pal$estimate
corr_results[13,10]<-picModel.x.low.y.pal$estimate


#################################
#      Y position lower jaw Vs x position palatine
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(x.pal, tree.final)
picModel.y.low.x.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[11,12]<-picModel.y.low.x.pal$estimate
corr_results[12,11]<-picModel.y.low.x.pal$estimate


#################################
#      Y position lower jaw Vs y position palatine
##################################
hPic <- pic(y.pos.jaw, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.y.low.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[11,13]<-picModel.y.low.y.pal$estimate
corr_results[13,11]<-picModel.y.low.y.pal$estimate



#################################
#      X position palatine Vs Y position palatine
##################################
hPic <- pic(x.pal, tree.final)
aPic <- pic(y.pal, tree.final)
picModel.x.pal.y.pal <- cor.test(hPic,aPic,method="pearson")

corr_results[12,13]<-picModel.x.pal.y.pal$estimate
corr_results[13,12]<-picModel.x.pal.y.pal$estimate














#############################################
#             Simulations
#############################################
require(phytools)

##################################
#  Simulating trait data under BM
#####################################
gape.sim<-fastBM(tree.final, nsim=1000)
prod.sim<-fastBM(tree.final, nsim=1000)
add.sim<-fastBM(tree.final, nsim=1000)
head.sim<-fastBM(tree.final, nsim=1000)
premax.sim<-fastBM(tree.final, nsim=1000)
max.sim<-fastBM(tree.final, nsim=1000)
low.jaw.sim<-fastBM(tree.final, nsim=1000)
close.out.sim<-fastBM(tree.final,nsim=1000)
open.out.sim<-fastBM(tree.final,  nsim=1000)
x.low.sim<-fastBM(tree.final,  nsim=1000)
y.low.sim<-fastBM(tree.final,  nsim=1000)
x.palatine.sim<-fastBM(tree.final, nsim=1000)
y.palatine.sim<-fastBM(tree.final,  nsim=1000)  



###########################################################
# Creating dataframe to store simulated data
############################################################

sim_results <- setNames(data.frame(matrix(ncol = 14, nrow = 110)), c("Species", "gape", "Jaw.Protrusion", "Adductor.Mass", "head.height.head.length", 
                                                                     "Premax.Length", "Maxilla.Length", "Lower.Jaw.Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                                                                     "x.palatine", "y.palatine"))
corr_results_sim <- setNames(data.frame(matrix(ncol = 78, nrow = 1000)), c("gape.vs.jaw.prod", "gape.add.mass","gape.vs.head", "gape.vs.premax", "gape.vs.max", "gape.vs.lower.lev",  
                                                                           "gape.vs.close", "gape.vs.open", "gape.vs.x.low", "gape.vs.y.low", "gape.vs.x.pal", "gape.vs.y.pal",
                                                                           "prod.vs.add", "prod.vs.head", "prod.vs.premax", "prod.vs.max", "prod.vs.low.lev", "prod.vs.close",
                                                                           "prod.vs.open", "prod.vs.x.low", "prod.vs.y.low", "prod.vs.x.pal", "prod.vs.y.pal", "add.vs.head",
                                                                           "add.vs.premax", "add.vs.max", "add.vs.lower.lev", "add.vs.close", "add.vs.open", "add.vs.x.low", 
                                                                           "add.vs.y.low", "add.vs.x.pal", "add.vs.y.pal", "head.vs.premax", "head.vs.max", "head.vs.lower.lev",  
                                                                           "head.vs.close", "head.vs.open", "head.vs.x.low", "head.vs.y.low", "head.vs.x.pal", "head.vs.y.pal",
                                                                           "premax.vs.max", "premax.vs.lower.lev", "premax.vs.close", "premax.vs.open", "premax.vs.x.low", "premax.vs.y.low",  
                                                                           "premax.vs.x.pal", "premax.vs.y.pal", "max.vs.lower.lev", "max.vs.close", "max.vs.open", "max.vs.x.low",
                                                                           "max.vs.y.low", "max.vs.x.pal", "max.vs.y.pal", "low.lev.vs.close", "low.lev.vs.open", "low.lev.vs.x.low",
                                                                           "low.lev.vs.y.low", "low.lev.vs.x.pal", "low.lev.vs.y.pal", "close.vs.open", "close.vs.x.low", "close.vs.y.low", 
                                                                           "close.vs.x.pal", "close.vs.y.pal", "open.vs.x.low", "open.vs.y.low", "open.vs.x.pal", "open.vs.y.pal", 
                                                                           "x.low.vs.y.low", "x.low.vs.x.pal", "x.low.vs.y.pal","y.low.vs.x.pal", "y.low.vs.y.pal","x.pal.vs.y.pal"))



####################################################################
#        Calculating correlation coefficient for each simulation
####################################################################

for(i in 1:1000){
  pic.gape<-pic(gape.sim[,i], tree.final)
  pic.prod<-pic(prod.sim[,i], tree.final)
  pic.add<-pic(add.sim[,i], tree.final)
  pic.head<- pic(head.sim[,i], tree.final)
  pic.premax<-pic(premax.sim[,i], tree.final)
  pic.max<-pic(max.sim[,i], tree.final)
  pic.low.jaw<-pic(low.jaw.sim[,i], tree.final)
  pic.close<-pic(close.out.sim[,i], tree.final)
  pic.open<- pic(open.out.sim[,i], tree.final)
  pic.x.low<-pic(x.low.sim[,i], tree.final)
  pic.y.low<-pic(y.low.sim[,i], tree.final)
  pic.x.pal<-pic(x.palatine.sim[,i], tree.final)
  pic.y.pal<-pic(y.palatine.sim[,i], tree.final)
  
 
  
  pic.x.pal.y.pal<-cor.test(pic.x.pal,  pic.y.pal, method="pearson")
  corr_results_sim$x.pal.vs.y.pal[i]<-pic.x.pal.y.pal$estimate[[1]]
  
  pic.y.low.x.pal<-cor.test(pic.y.low, pic.x.pal, method="pearson")
  corr_results_sim$y.low.vs.x.pal[i]<-pic.y.low.x.pal$estimate[[1]]
  
  pic.y.low.y.pal<-cor.test(pic.y.low, pic.y.pal, method="pearson")
  corr_results_sim$y.low.vs.y.pal[i]<-pic.y.low.y.pal$estimate[[1]]
  
  
  pic.x.low.y.low<-cor.test(pic.x.low, pic.y.low, method="pearson")
  corr_results_sim$x.low.vs.y.low[i]<-pic.x.low.y.low$estimate[[1]]
  
  pic.x.low.x.pal<-cor.test(pic.x.low, pic.x.pal, method="pearson")
  corr_results_sim$x.low.vs.x.pal[i]<-pic.x.low.x.pal$estimate[[1]]
  
  pic.x.low.y.pal<-cor.test(pic.x.low, pic.y.pal, method="pearson")
  corr_results_sim$x.low.vs.y.pal[i]<-pic.x.low.y.pal$estimate[[1]]
  
  
  pic.open.x.low<-cor.test(pic.open, pic.x.low, method="pearson")
  corr_results_sim$open.vs.x.low[i]<-pic.open.x.low$estimate[[1]]
  
  pic.open.y.low<-cor.test(pic.open, pic.y.low,method="pearson")
  corr_results_sim$open.vs.y.low[i]<-pic.open.y.low$estimate[[1]]
  
  pic.open.x.pal<-cor.test(pic.open, pic.x.pal, method="pearson")
  corr_results_sim$open.vs.x.pal[i]<-pic.open.x.pal$estimate[[1]]
  
  pic.open.y.pal<-cor.test(pic.open, pic.y.pal, method="pearson")
  corr_results_sim$open.vs.y.pal[i]<-pic.open.y.pal$estimate[[1]]
  
  
  pic.close.open<-cor.test(pic.close, pic.open, method="pearson")
  corr_results_sim$close.vs.open[i]<-pic.close.open$estimate[[1]]
  
  pic.close.x.low<-cor.test(pic.close, pic.x.low, method="pearson")
  corr_results_sim$close.vs.x.low[i]<-pic.close.x.low$estimate[[1]]
  
  pic.close.y.low<-cor.test(pic.close, pic.y.low, method="pearson")
  corr_results_sim$close.vs.y.low[i]<-pic.close.y.low$estimate[[1]]
  
  pic.close.x.pal<-cor.test(pic.close, pic.x.pal, method="pearson")
  corr_results_sim$close.vs.x.pal[i]<-pic.close.x.pal$estimate[[1]]
  
  pic.close.y.pal<-cor.test(pic.close, pic.y.pal, method="pearson")
  corr_results_sim$close.vs.y.pal[i]<-pic.close.y.pal$estimate[[1]]
  
  
  pic.low.jaw.close<-cor.test(pic.low.jaw, pic.close, method="pearson")
  corr_results_sim$low.lev.vs.close[i]<-pic.low.jaw.close$estimate[[1]]
  
  pic.low.jaw.open<-cor.test(pic.low.jaw, pic.open, method="pearson")
  corr_results_sim$low.lev.vs.open[i]<-pic.low.jaw.open$estimate[[1]]
  
  pic.low.jaw.x.low<-cor.test(pic.low.jaw, pic.x.low, method="pearson")
  corr_results_sim$low.lev.vs.x.low[i]<-pic.low.jaw.x.low$estimate[[1]]
  
  pic.low.jaw.y.low<-cor.test(pic.low.jaw, pic.y.low, method="pearson")
  corr_results_sim$low.lev.vs.y.low[i]<-pic.low.jaw.y.low$estimate[[1]]
  
  pic.low.jaw.x.pal<-cor.test(pic.low.jaw, pic.x.pal, method="pearson")
  corr_results_sim$low.lev.vs.x.pal[i]<-pic.low.jaw.x.pal$estimate[[1]]
  
  pic.low.jaw.y.pal<-cor.test(pic.low.jaw, pic.y.pal, method="pearson")
  corr_results_sim$low.lev.vs.y.pal[i]<-pic.low.jaw.y.pal$estimate[[1]]
  
  
  
  pic.max.low.jaw<-cor.test(pic.max, pic.low.jaw, method="pearson")
  corr_results_sim$max.vs.lower.lev[i]<-pic.max.low.jaw$estimate[[1]]
  
  pic.max.close<-cor.test(pic.max, pic.close, method="pearson")
  corr_results_sim$max.vs.close[i]<-pic.max.close$estimate[[1]]
  
  pic.max.open<-cor.test(pic.max, pic.open, method="pearson")
  corr_results_sim$max.vs.open[i]<-pic.max.open$estimate[[1]]
  
  pic.max.x.low<-cor.test(pic.max, pic.x.low, method="pearson")
  corr_results_sim$max.vs.x.low[i]<-pic.max.x.low$estimate[[1]]
  
  pic.max.y.low<-cor.test(pic.max, pic.y.low, method="pearson")
  corr_results_sim$max.vs.y.low[i]<-pic.max.y.low$estimate[[1]]
  
  pic.max.x.pal<-cor.test(pic.max, pic.x.pal, method="pearson")
  corr_results_sim$max.vs.x.pal[i]<-pic.max.x.pal$estimate[[1]]
  
  pic.max.y.pal<-cor.test(pic.max, pic.y.pal, method="pearson")
  corr_results_sim$max.vs.y.pal[i]<-pic.max.y.pal$estimate[[1]]
  
  
  pic.premax.max<-cor.test(pic.premax, pic.max, method="pearson")
  corr_results_sim$premax.vs.max[i]<-pic.premax.max$estimate[[1]]
  
  pic.premax.low.jaw<-cor.test(pic.premax, pic.low.jaw, method="pearson")
  corr_results_sim$premax.vs.lower.lev[i]<-pic.premax.low.jaw$estimate[[1]]
  
  pic.premax.close<-cor.test(pic.premax, pic.close, method="pearson")
  corr_results_sim$premax.vs.close[i]<-pic.premax.close$estimate[[1]]
  
  pic.premax.open<-cor.test(pic.premax, pic.open, method="pearson")
  corr_results_sim$premax.vs.open[i]<-pic.premax.open$estimate[[1]]
  
  pic.premax.x.low<-cor.test(pic.premax, pic.x.low, method="pearson")
  corr_results_sim$premax.vs.x.low[i]<-pic.premax.x.low$estimate[[1]]
  
  pic.premax.y.low<-cor.test(pic.premax, pic.y.low, method="pearson")
  corr_results_sim$premax.vs.y.low[i]<-pic.premax.y.low$estimate[[1]]
  
  pic.premax.x.pal<-cor.test(pic.premax, pic.x.pal, method="pearson")
  corr_results_sim$premax.vs.x.pal[i]<-pic.premax.x.pal$estimate[[1]]
  
  pic.premax.y.pal<-cor.test(pic.premax, pic.y.pal, method="pearson")
  corr_results_sim$premax.vs.y.pal[i]<-pic.premax.y.pal$estimate[[1]]
  
  
  
  pic.head.premax<-cor.test(pic.head, pic.premax, method="pearson")
  corr_results_sim$head.vs.premax[i]<-pic.head.premax$estimate[[1]]
  
  pic.head.max<-cor.test(pic.head, pic.max, method="pearson")
  corr_results_sim$head.vs.max[i]<-pic.head.max$estimate[[1]]
  
  pic.head.low.jaw<-cor.test(pic.head, pic.low.jaw, method="pearson")
  corr_results_sim$head.vs.lower.lev[i]<-pic.head.low.jaw$estimate[[1]]
  
  pic.head.close<-cor.test(pic.head, pic.close, method="pearson")
  corr_results_sim$head.vs.close[i]<-pic.head.close$estimate[[1]]
  
  pic.head.open<-cor.test(pic.head, pic.open, method="pearson")
  corr_results_sim$head.vs.open[i]<-pic.head.open$estimate[[1]]
  
  pic.head.x.low<-cor.test(pic.head, pic.x.low, method="pearson")
  corr_results_sim$head.vs.x.low[i]<-pic.head.x.low$estimate[[1]]
  
  pic.head.y.low<-cor.test(pic.head, pic.y.low, method="pearson")
  corr_results_sim$head.vs.y.low[i]<-pic.head.y.low$estimate[[1]]
  
  pic.head.x.pal<-cor.test(pic.head, pic.x.pal, method="pearson")
  corr_results_sim$head.vs.x.pal[i]<-pic.head.x.pal$estimate[[1]]
  
  pic.head.y.pal<-cor.test(pic.head, pic.y.pal, method="pearson")
  corr_results_sim$head.vs.y.pal[i]<-pic.head.y.pal$estimate[[1]]
  
  
  pic.add.head<-cor.test(pic.add, pic.head, method="pearson")
  corr_results_sim$add.vs.head[i]<-pic.add.head$estimate[[1]]
  
  pic.add.premax<-cor.test(pic.add, pic.premax, method="pearson")
  corr_results_sim$add.vs.premax[i]<-pic.add.premax$estimate[[1]]
  
  pic.add.max<-cor.test(pic.add, pic.max, method="pearson")
  corr_results_sim$add.vs.max[i]<-pic.add.max$estimate[[1]]
  
  pic.add.low.jaw<-cor.test(pic.add, pic.low.jaw, method="pearson")
  corr_results_sim$add.vs.lower.lev[i]<-pic.add.low.jaw$estimate[[1]]
  
  pic.add.close<-cor.test(pic.add, pic.close, method="pearson")
  corr_results_sim$add.vs.close[i]<-pic.add.close$estimate[[1]]
  
  pic.add.open<-cor.test(pic.add, pic.open, method="pearson")
  corr_results_sim$add.vs.open[i]<-pic.add.open$estimate[[1]]
  
  pic.add.x.low<-cor.test(pic.add, pic.x.low, method="pearson")
  corr_results_sim$add.vs.x.low[i]<-pic.add.x.low$estimate[[1]]
  
  pic.add.y.low<-cor.test(pic.add, pic.y.low, method="pearson")
  corr_results_sim$add.vs.y.low[i]<-pic.add.y.low$estimate[[1]]
  
  pic.add.x.pal<-cor.test(pic.add, pic.x.pal, method="pearson")
  corr_results_sim$add.vs.x.pal[i]<-pic.add.x.pal$estimate[[1]]
  
  pic.add.y.pal<-cor.test(pic.add, pic.y.pal, method="pearson")
  corr_results_sim$add.vs.y.pal[i]<-pic.add.y.pal$estimate[[1]]
  
  
  pic.prod.add<-cor.test(pic.prod, pic.add, method="pearson")
  corr_results_sim$prod.vs.add[i]<-pic.prod.add$estimate[[1]]
  
  pic.prod.head<-cor.test(pic.prod, pic.head, method="pearson")
  corr_results_sim$prod.vs.head[i]<-pic.prod.head$estimate[[1]]
  
  pic.prod.premax<-cor.test(pic.prod, pic.premax, method="pearson")
  corr_results_sim$prod.vs.premax[i]<-pic.prod.premax$estimate[[1]]
  
  pic.prod.max<-cor.test(pic.prod, pic.max, method="pearson")
  corr_results_sim$prod.vs.max[i]<-pic.prod.max$estimate[[1]]
  
  pic.prod.low<-cor.test(pic.prod, pic.low.jaw, method="pearson")
  corr_results_sim$prod.vs.low.lev[i]<-pic.prod.low$estimate[[1]]
  
  pic.prod.close<-cor.test(pic.prod, pic.close, method="pearson")
  corr_results_sim$prod.vs.close[i]<-pic.prod.close$estimate[[1]]
  
  pic.prod.open<-cor.test(pic.prod, pic.open, method="pearson")
  corr_results_sim$prod.vs.open[i]<-pic.prod.open$estimate[[1]]
  
  pic.prod.x.low<-cor.test(pic.prod, pic.x.low, method="pearson")
  corr_results_sim$prod.vs.x.low[i]<-pic.prod.x.low$estimate[[1]]
  
  pic.prod.y.low<-cor.test(pic.prod, pic.y.low, method="pearson")
  corr_results_sim$prod.vs.y.low[i]<-pic.prod.y.low$estimate[[1]]
  
  pic.prod.x.pal<-cor.test(pic.prod, pic.x.pal, method="pearson")
  corr_results_sim$prod.vs.x.pal[i]<-pic.prod.x.pal$estimate[[1]]
  
  pic.prod.y.pal<-cor.test(pic.prod, pic.y.pal, method="pearson")
  corr_results_sim$prod.vs.y.pal[i]<-pic.prod.y.pal$estimate[[1]]
  
  pic.gape.prod<-cor.test(pic.gape, pic.prod, method="pearson")
  corr_results_sim$gape.vs.jaw.prod[i]<-pic.gape.prod$estimate[[1]]
  
  pic.gape.add.mass<-cor.test(pic.gape, pic.add, method="pearson")
  corr_results_sim$gape.add.mass[i]<-pic.gape.add.mass$estimate[[1]]
  
  pic.gape.head<-cor.test(pic.gape, pic.head, method="pearson")
  corr_results_sim$gape.vs.head[i]<-pic.gape.head$estimate[[1]]
  
  pic.gape.premax<-cor.test(pic.gape, pic.premax, method="pearson")
  corr_results_sim$gape.vs.premax[i]<-pic.gape.premax$estimate[[1]]
  
  pic.gape.max<-cor.test(pic.gape, pic.max, method="pearson")
  corr_results_sim$gape.vs.max[i]<-pic.gape.max$estimate[[1]]
  
  pic.gape.low.jaw<-cor.test(pic.gape, pic.low.jaw, method="pearson")
  corr_results_sim$gape.vs.lower.lev[i]<-pic.gape.low.jaw$estimate[[1]]
  
  pic.gape.close<-cor.test(pic.gape, pic.close, method="pearson")
  corr_results_sim$gape.vs.close[i]<-pic.gape.close$estimate[[1]]
  
  pic.gape.open<-cor.test(pic.gape, pic.open, method="pearson")
  corr_results_sim$gape.vs.open[i]<-pic.gape.open$estimate[[1]]
  
  pic.gape.x.low<-cor.test(pic.gape, pic.x.low, method="pearson")
  corr_results_sim$gape.vs.x.low[i]<-pic.gape.x.low$estimate[[1]]
  
  pic.gape.y.low<-cor.test(pic.gape, pic.y.low, method="pearson")
  corr_results_sim$gape.vs.y.low[i]<-pic.gape.y.low$estimate[[1]]
  
  pic.gape.x.pal<-cor.test(pic.gape, pic.x.pal, method="pearson")
  corr_results_sim$gape.vs.x.pal[i]<-pic.gape.x.pal$estimate[[1]]
  
  pic.gape.y.pal<-cor.test(pic.gape, pic.y.pal, method="pearson")
  corr_results_sim$gape.vs.y.pal[i]<-pic.gape.y.pal$estimate[[1]]
  
  
  
  
  
  
}

#######################################################
#      Calculate max simulated correlation per trait
#######################################################

colMax <- function(data) sapply(data, max, na.rm = TRUE)

corr_results_sim_max<-colMax(corr_results_sim)


#######################################################################
# Remove empirical correlations that are not larger than simulated max
#######################################################################


# set up the output correlation matrix
corr_res_sig <- corr_results

# take just the lower part of the correlation matrix
lwr <- lower.tri(corr_results, diag = F)
# find the indices of the TRUE because those will be compared
lwrIdx <- which(lwr, arr.ind = TRUE)
# loop through lower half of matrix and put zeros if sim < empirical
for (f in 1:length(corr_results_sim_max)){
  simMax <- corr_results_sim_max[f]
  emp <- corr_results[lwrIdx[f,1], lwrIdx[f,2]]
  if (emp < simMax){
    corr_res_sig[lwrIdx[f,1], lwrIdx[f,2]] = 0
  }
}

##################################
#   Correlation heatmap
##################################
require(corrplot)



corrplot::corrplot(as.matrix(corr_res_sig), method = 'circle', type = 'upper', order = "hclust", 
                   tl.col = "black", tl.srt = 45, is.corr=FALSE)




################################
#        Bivariate plots
################################

#pdf(file="Plot_lower jaw vs gape_contrasts.pdf", width=8, height=8)
plot(hPic ~ aPic, pch=16, xlab="Lower jaw length contrasts", ylab="Mouth gape contrasts")
abline(a = 0, b = picModel.gape.low$estimate)
#dev.off()


#pdf(file="Plot_premaxilla vs gape_contrasts.pdf", width=8, height=8)
plot(hPic ~ aPic, pch=16, xlab="Premaxilla length contrasts", ylab="Mouth gape contrasts")
abline(a = 0, b = picModel.gape.premax$estimate)
#dev.off()


#pdf(file="Plot_addcutor vs open ma_contrasts.pdf", width=8, height=8)
plot(hPic ~ aPic, pch=16, xlab="Opening mechanical advantage contrasts", ylab="Adductor muscle mass contrasts")
abline(a = 0, b = picModel.add.open$estimate)
#dev.off()


#pdf(file="Plot_addcutor vs y palatine_contrasts.pdf", width=8, height=8)
plot(hPic ~ aPic, pch=16, xlab="Adductor muscle mass contrasts", ylab="Y position of the palatine contrasts")
abline(a = 0, b = cor.test(hPic,aPic,method="pearson"))
#dev.off()































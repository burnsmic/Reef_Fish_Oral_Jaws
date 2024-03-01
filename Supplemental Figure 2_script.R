
######################################################
#      Morphological disparity of each trait for diet
#######################################################

disparity_results <- setNames(data.frame(matrix(ncol = 8, nrow = 13)), c("trait","MIB_disparity", "GC_disparity", "HD_disparity", "MIS_disparity",
                                                                         "OM_disparity", "PK_disparity", "SI_disparity"))
   
diet<-dat.rat.na.log$Diet

names(diet)<-dat.rat.na.log$Species 

i=1

for(i in 1:13){
  print(paste("Now starting trait:",trait_vector[i]))
  feeding <- dat.rat.na$Diet
  now_trait <- as.data.frame(dat.rat.na.log[,(i+4)])
  row.names(now_trait) <- dat.rat.na.log$Species
  feeding_disparity <- morphol.disparity(as.matrix(now_trait) ~ diet, groups= ~ diet, iter = 10000, print.progress = FALSE)
  disparity_results$trait[i] <- trait_vector[i]
  disparity_results$MIB_disparity[i] <- feeding_disparity$Procrustes.var[1]
  disparity_results$GC_disparity[i] <- feeding_disparity$Procrustes.var[2]
  disparity_results$HD_disparity[i] <- feeding_disparity$Procrustes.var[3]
  disparity_results$MIS_disparity[i] <- feeding_disparity$Procrustes.var[4]
  disparity_results$OM_disparity[i] <- feeding_disparity$Procrustes.var[5]
  disparity_results$PK_disparity[i] <- feeding_disparity$Procrustes.var[6]
  disparity_results$SI_disparity[i] <- feeding_disparity$Procrustes.var[7]
}



dr_pca = setNames(data.frame(names(disparity_results)[-1], 
                             t(unname(disparity_results[,-1]))), 
                  c('Diet', disparity_results[,1]))



# bit by bit
# pull out the column names, but skip the first one (we don't need it)
col_names = names(disparity_results)[-1]
# trim the data frame to remove the first column and remove the column headers
disparity_results_simp = unname(disparity_results[,-1])
# transpose the numbers only
disparity_results_t = t(disparity_results_simp)
# combine the original column names and the transposed data frame into a new 
# data frame with the column names as the first column
disparity_results_pca = data.frame(col_names, disparity_results_t)
# assign new column names based on the original data frames first row
disparity_results_pca = setNames(disparity_results_pca, c('Diet', disparity_results[,1]))



###############
#  PCA
###############

disparity.pca<-prcomp(disparity_results_pca[,2:14], scale=TRUE)



disparity_results_pca$Color<-disparity_results_pca$Diet
disparity_results_pca$Color[disparity_results_pca$Color == "GC_disparity"]<-"#F9A21E"
disparity_results_pca$Color[disparity_results_pca$Color == "OM_disparity"]<-"#F06A95"
disparity_results_pca$Color[disparity_results_pca$Color == "HD_disparity"]<-"#11619B"
disparity_results_pca$Color[disparity_results_pca$Color == "MIS_disparity"]<-"#27B79B"
disparity_results_pca$Color[disparity_results_pca$Color == "SI_disparity"]<-"#47ADE1"
disparity_results_pca$Color[disparity_results_pca$Color == "PK_disparity"]<-"#9B519E" 
disparity_results_pca$Color[disparity_results_pca$Color == "MIB_disparity"]<-"#E21C4D" 


plot(disparity.pca$x[,1], disparity.pca$x[,2], pch=16, col=disparity_results_pca$Color,
     xlab="PC1", ylab="PC2", main="PC1 vs. PC2")

########################
# Supplemental Table 3
########################

disparity.pca$rotation

##############
# Dotplot
##############

dotchart(log(disparity_results[,2]), pch = 21, labels = disparity_results$trait, bg = "#E21C4D", xlab="Morphological variance (log transformed)",
         pt.cex = 1.5, xlim=c(-8,0))
points(log(disparity_results[,3]), 1:nrow(dat.dot), col = "#F9A21E", pch = 19, cex = 1.5)
points(log(disparity_results[,4]), 1:nrow(dat.dot), col = "#11619B", pch = 19, cex = 1.5)
points(log(disparity_results[,5]), 1:nrow(dat.dot), col = "#27B79B", pch = 19, cex = 1.5)
points(log(disparity_results[,6]), 1:nrow(dat.dot), col = "#F06A95", pch = 19, cex = 1.5)
points(log(disparity_results[,7]), 1:nrow(dat.dot), col = "#9B519E", pch = 19, cex = 1.5)
points(log(disparity_results[,8]), 1:nrow(dat.dot), col = "#47ADE1", pch = 19, cex = 1.5)




##################################################
#      Morphological rate of each trait for diet
##################################################

rate_results <- setNames(data.frame(matrix(ncol = 8, nrow = 13)), c("trait","MIB_disparity", "GC_disparity", "HD_disparity", "MIS_disparity",
                                                                         "OM_disparity", "PK_disparity", "SI_disparity"))

diet<-dat.rat.na.log$Diet

names(diet)<-dat.rat.na.log$Species

i=1

for(i in 1:13){
  print(paste("Now starting trait:",trait_vector[i]))
  feeding <- dat.rat.na$Diet
  now_trait <- as.data.frame(dat.rat.na.log[,(i+4)])
  row.names(now_trait) <- dat.rat.na.log$Species
  feeding_rate <- compare.evol.rates(A=as.matrix(now_trait), phy=tree.final, gp=diet, iter = 9999, seed=NULL, method=c("simulation"), print.progress = FALSE)
  rate_results$trait[i] <- trait_vector[i]
  rate_results$MIB_disparity[i] <- feeding_rate$sigma.d.gp[1]
  rate_results$GC_disparity[i] <- feeding_rate$sigma.d.gp[2]
  rate_results$HD_disparity[i] <- feeding_rate$sigma.d.gp[3]
  rate_results$MIS_disparity[i] <- feeding_rate$sigma.d.gp[4]
  rate_results$OM_disparity[i] <- feeding_rate$sigma.d.gp[5]
  rate_results$PK_disparity[i] <- feeding_rate$sigma.d.gp[6]
  rate_results$SI_disparity[i] <- feeding_rate$sigma.d.gp[7]
}



rate_pca = setNames(data.frame(names(rate_results)[-1], 
                             t(unname(rate_results[,-1]))), 
                  c('Diet', rate_results[,1]))



# bit by bit
# pull out the column names, but skip the first one (we don't need it)
col_names = names(rate_results)[-1]
# trim the data frame to remove the first column and remove the column headers
rate_results_simp = unname(rate_results[,-1])
# transpose the numbers only
rate_results_t = t(rate_results_simp)
# combine the original column names and the transposed data frame into a new 
# data frame with the column names as the first column
rate_results_pca = data.frame(col_names, rate_results_t)
# assign new column names based on the original data frames first row
rate_results_pca = setNames(rate_results_pca, c('Diet', rate_results[,1]))



###############
#  PCA
###############

rate.pca<-prcomp(rate_results_pca[,2:14], scale=TRUE)


rate_results_pca$Color<-rate_results_pca$Diet
rate_results_pca$Color[rate_results_pca$Color == "GC_disparity"]<-"#F9A21E"
rate_results_pca$Color[rate_results_pca$Color == "OM_disparity"]<-"#F06A95"
rate_results_pca$Color[rate_results_pca$Color == "HD_disparity"]<-"#11619B"
rate_results_pca$Color[rate_results_pca$Color == "MIS_disparity"]<-"#27B79B"
rate_results_pca$Color[rate_results_pca$Color == "SI_disparity"]<-"#47ADE1"
rate_results_pca$Color[rate_results_pca$Color == "PK_disparity"]<-"#9B519E" 
rate_results_pca$Color[rate_results_pca$Color == "MIB_disparity"]<-"#E21C4D" 


plot(rate.pca$x[,1], rate.pca$x[,2], pch=16, col=rate_results_pca$Color,
     xlab="PC1", ylab="PC2", main="PC1 vs. PC2")

########################
# Supplemental Table 3
########################

rate.pca$rotation

##############
# Dotplot
##############

dotchart(log(rate_results[,2]), pch = 21, labels = rate_results$trait, bg = "#E21C4D", xlab="Rate of evolution (log transformed)",
         pt.cex = 1.5, xlim=c(-11,-3))
points(log(rate_results[,3]), 1:nrow(dat.dot), col = "#F9A21E", pch = 19, cex = 1.5)
points(log(rate_results[,4]), 1:nrow(dat.dot), col = "#11619B", pch = 19, cex = 1.5)
points(log(rate_results[,5]), 1:nrow(dat.dot), col = "#27B79B", pch = 19, cex = 1.5)
points(log(rate_results[,6]), 1:nrow(dat.dot), col = "#F06A95", pch = 19, cex = 1.5)
points(log(rate_results[,7]), 1:nrow(dat.dot), col = "#9B519E", pch = 19, cex = 1.5)
points(log(rate_results[,8]), 1:nrow(dat.dot), col = "#47ADE1", pch = 19, cex = 1.5)






################################################
#            Calculating Multivariate Disparity- functional feeding mode
################################################

require(geomorph)

trop.func<-dat.rat.na.log$Trophic
names(trop.func)<-dat.rat.na.log$Species

traits <- as.data.frame(dat.rat.na.log[,c(5:17)])
rownames(traits)<-dat.rat.na.log$Species

morph_disparity <- morphol.disparity(as.matrix(traits) ~ trop.func, groups= ~ trop.func, iter = 99999, print.progress = TRUE)
summary(morph_disparity)

##########################################
#               Multivariate Rate analysis
##########################################

morph_rates <- compare.evol.rates(
  A = as.matrix(traits),
  phy = tree.final,
  gp = trop.func,
  iter = 9999,
  seed = NULL,
  method = c("simulation"),
  print.progress = TRUE
)

summary(morph_rates)

###########################################
#               MANOVA
###########################################

morph_manova <- procD.pgls(as.matrix(traits) ~ trop.func, phy = tree.final, iter = 9999)

summary(morph_manova)

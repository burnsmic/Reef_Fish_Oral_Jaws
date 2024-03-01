##########################################
#         Pairwise Disparity
##########################################

trait_vector <- c("gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", 
                  "Premax Length", "Maxilla Length", "Lower Jaw Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")

disparities_list <- as.data.frame(trait_vector)

disparities_results <- setNames(data.frame(matrix(ncol = 4, nrow = 13)), c("trait","biting_disparity", "suction_disparity", "pval"))
i=3

for(i in 1:13){
  #print(paste("Now starting trait:",trait_vector[i]))
  feeding <- dat.rat.na.log$Trophic
  now_trait <- as.data.frame(dat.rat.na.log[,(i+4)])
  row.names(now_trait) <- dat.rat.na.log$Species
  feeding_disparity <- morphol.disparity(as.matrix(now_trait) ~ trop.func, groups= ~ trop.func, iter = 10000, print.progress = FALSE)
  disparities_results$trait[i] <- trait_vector[i]
  disparities_results$biting_disparity[i] <- feeding_disparity$Procrustes.var[1]
  disparities_results$suction_disparity[i] <- feeding_disparity$Procrustes.var[2]
  disparities_results$pval[i]<-feeding_disparity$PV.dist.Pval[2]
}


disparities_results


###########
# Dotplot
###########

#pdf(file="DotPlot_disparity_function.pdf", width=12, height=8)
dotchart(log(disparities_results[,2]), pch = 21, labels = disparities_results$trait, bg = "black", xlab="Trait values",
         pt.cex = 1.5, xlim=c(-7,0.5))
points(log(disparities_results[,3]), 1:nrow(dat.dot), col = "red", pch = 19, cex = 1.5)
#dev.off()


########################################################
#          Rate analysis - individual traits
#########################################################

trait_vector <- c("gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", 
                  "Premax Length", "Maxilla Length", "Lower Jaw Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")

rate_list <- as.data.frame(trait_vector)

rate_results <- setNames(data.frame(matrix(ncol = 4, nrow = 13)), c("trait","biting_rate", "suction_rate", "pval"))
i=1

for(i in 1:13){
  #print(paste("Now starting trait:",trait_vector[i]))
  feeding <- dat.rat.na$Diet
  now_trait <- as.data.frame(dat.rat.na.log[,(i+4)])
  row.names(now_trait) <- dat.rat.na$Species
  feeding_rate <- compare.evol.rates(as.matrix(now_trait), phy=tree.final, gp= trop.func, iter = 9999, print.progress = FALSE)
  rate_results$trait[i] <- trait_vector[i]
  rate_results$biting_rate[i] <- feeding_rate$sigma.d.gp[1]
  rate_results$suction_rate[i] <- feeding_rate$sigma.d.gp[2]
  rate_results$pval[i]<-feeding_rate$P.value
}

rate_results

###########
# Dotplot
###########

#pdf(file="DotPlot_rate_function.pdf", width=12, height=8)
dotchart(log(rate_results[,2]), pch = 21, labels = rate_results$trait, bg = "black", xlab="Trait values",
         pt.cex = 1.5, xlim=c(-11,-4))
points(log(rate_results[,3]), 1:nrow(dat.dot), col = "red", pch = 19, cex = 1.5)
#dev.off()


##################################
#       Mean analysis
##################################

traits <- as.data.frame(dat.rat.na.log[,c(5:17)])
rownames(traits)<-dat.rat.na.log$Species

trait_troph<-traits

trait_troph$Trophic<-dat.rat.na.log$Trophic

aovs_table <- setNames(data.frame(matrix(ncol = 4, nrow = 13)), c("trait", "Biter_Mean", "Suction_Mean","p.val"))
i=1
anovas <- for(i in 1:13){
  print(paste("Now starting trait:",trait_vector[i]))
  trait.ind<-traits[,i]
  names(trait.ind)<-rownames(traits)
  morph_anova <- procD.pgls(trait.ind ~ trop.func, phy = tree.final, iter = 9999)
  #aovs_list[i,2] <- list(morph_anova)
  aovs_table$trait[i] <- trait_vector[i]
  aovs_table$Biter_Mean[i] <- mean(trait_troph[,i][which(dat.rat.na$Trophic == "biter")]) 
  aovs_table$Suction_Mean[i] <- mean(trait_troph[,i][which(dat.rat.na$Trophic == "suction")])
}

aovs_table


###########
# Dotplot
###########

#pdf(file="DotPlot_mean_function.pdf", width=12, height=8)
dotchart(aovs_table[,2], pch = 21, labels = aovs_table$trait, bg = "black", xlab="Trait values",
         pt.cex = 1.5, xlim=c(-6,1))
points(aovs_table[,3], 1:nrow(dat.dot), col = "red", pch = 19, cex = 1.5)
#dev.off()

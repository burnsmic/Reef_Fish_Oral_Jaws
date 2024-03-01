library("tidyverse")
require(phytools)

######################################################
#         Calculating simmap of diet
#######################################################
trop.diet<-dat.rat.na$Diet
names(trop.diet)<-dat.rat.na$Species

func_simmaps <- make.simmap(tree = tree.final, x=trop.diet, Q="empirical", pi = "estimated", nsim = 1000)



##################################################
#   Making a dataframe for OUwie analysis
##################################################

#change df order to [,1] species, [,2] current selective regime, and [,3] continuous trait
trait_vector <- c("gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", 
                  "Premax Length", "Maxilla Length", "Lower Jaw Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")
#colnames(final.trait.data)[1:8] <- trait_vector 
#OUwie.df <- rownames_to_column(final.trait.data, var = "Species") %>% relocate(Habitat, .before = Standard_length)
#save(OUwie.df, file="OUwie.df.RData")

OUwie.df<-setNames(data.frame(matrix(ncol = 15, nrow = 110)), c("Species","Diet","gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", "Premax Length", 
                                                                "Maxilla Length", "Lower Jaw Lever", "close.out", "open.out", "x.low.jaw",
                                                              "y.low.jaw", "x.palatine", "y.palatine"))

OUwie.df[,1]<-dat.rat.na.log$Species
OUwie.df[,2]<-dat.rat.na.log$Diet
OUwie.df[,3]<-dat.rat.na.log$log.Gape.head.length
OUwie.df[,4]<-dat.rat.na.log$log.prod.head.length
OUwie.df[,5]<-dat.rat.na.log$log.add.mass.body.mass
OUwie.df[,6]<-dat.rat.na.log$log.head.height.head.length
OUwie.df[,7]<-dat.rat.na.log$log.premax.head.length
OUwie.df[,8]<-dat.rat.na.log$log.max.head.length
OUwie.df[,9]<-dat.rat.na.log$log.low.jaw.lever
OUwie.df[,10]<-dat.rat.na.log$log.close.out
OUwie.df[,11]<-dat.rat.na.log$log.open.out
OUwie.df[,12]<-dat.rat.na.log$log.x.low.jaw.head.length
OUwie.df[,13]<-dat.rat.na.log$log.y.low.jaw.head.length
OUwie.df[,14]<-dat.rat.na.log$log.x.palatine.head.length
OUwie.df[,15]<-dat.rat.na.log$log.y.palatine.head.length


###################  Set up data frame for results #############################

best_model_data <- tibble(.rows = 13000) # 13 traits * 1 models * 1000 simmaps
best_model_data$trait <- rep(NA,13000)
best_model_data$simmap_count <- rep(NA,13000)
best_model_data$model <- rep(NA,13000)
best_model_data$loglik <- rep(NA,13000)
best_model_data$aicc <- rep(NA,13000)
best_model_data$eigval <- rep(NA,13000)
best_model_data$sigma.sq_d <- rep(NA, 13000)
best_model_data$sigma.sq_gc <- rep(NA, 13000)
best_model_data$sigma.sq_HD <- rep(NA, 13000)
best_model_data$sigma.sq_MI <- rep(NA, 13000)
best_model_data$sigma.sq_OM <- rep(NA, 13000)
best_model_data$sigma.sq_PK <- rep(NA, 13000)
best_model_data$sigma.sq_SI <- rep(NA, 13000)
best_model_data$alpha_d <- rep(NA,13000)
best_model_data$alpha_gc <- rep(NA,13000)
best_model_data$alpha_HD <- rep(NA,13000)
best_model_data$alpha_MI <- rep(NA,13000)
best_model_data$alpha_OM <- rep(NA,13000)
best_model_data$alpha_PK <- rep(NA,13000)
best_model_data$alpha_SI <- rep(NA,13000)
best_model_data$theta_d <- rep(NA,13000)
best_model_data$theta_gc <- rep(NA,13000)
best_model_data$theta_HD <- rep(NA,13000)
best_model_data$theta_MI <- rep(NA,13000)
best_model_data$theta_OM <- rep(NA,13000)
best_model_data$theta_PK <- rep(NA,13000)
best_model_data$theta_SI <- rep(NA,13000)
best_model_data$theta_d_se <- rep(NA,13000)
best_model_data$theta_gc_se <- rep(NA,13000)
best_model_data$theta_HD_se <- rep(NA,13000)
best_model_data$theta_MI_se <- rep(NA,13000)
best_model_data$theta_OM_se <- rep(NA,13000)
best_model_data$theta_PK_se <- rep(NA,13000)
best_model_data$theta_SI_se <- rep(NA,13000)
best_model_data$saddle <- rep(NA,13000)

#########################  Run OUwie  ##########################################

library(OUwie)

# define function
OUwie.model <- function(model, phy, data) {
  print(paste("Now starting model: ", model))
  return(OUwie(phy, data, model, simmap.tree=T,algorithm="invert", diagn=T)) 
}

models <- c("OUM")

trait_vector <- c("gape", "Jaw.Protrusion", "Adductor.Mass", "head.height.head.length", 
                  "Premax.Length", "Maxilla.Length", "Lower.Jaw.Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")

j = 1 #trait
i = 1 #simmap
r = 1 
q = 1 #model

# subset data frame for each trait
for (j in 1:13){
  print(paste("Now starting trait: ", trait_vector[j]))
  OUwie_data <- data.frame(OUwie.df[ ,c(1,2,j+2)]) 
  i = 1
  
  # extract simmap to use
  for (i in 1:1000){
    print(paste("Now starting SIMMAP: ", i))
    tree <- func_simmaps[[i]]
    
    # apply function to each model
    results <- lapply(models, OUwie.model, phy=tree, data=OUwie_data) 
    q = 1
    
    for (q in 1:1) {
      each_model <- results[[q]] 
      
      # write results into data frame
      best_model_data$trait[r] <- trait_vector[j]
      best_model_data$simmap_count[r] <- i
      best_model_data$model[r] <- each_model$model
      best_model_data$loglik[r] <- each_model$loglik
      best_model_data$aicc[r] <- each_model$AICc
      best_model_data$eigval[r] <- list(each_model$eigval)
      best_model_data$sigma.sq_d[r] <- each_model$solution[2,1]
      best_model_data$sigma.sq_gc[r] <- each_model$solution[2,2]
      best_model_data$sigma.sq_HD[r] <- each_model$solution[2,3]
      best_model_data$sigma.sq_MI[r] <- each_model$solution[2,4]
      best_model_data$sigma.sq_OM[r] <- each_model$solution[2,5]
      best_model_data$sigma.sq_PK[r] <- each_model$solution[2,6]
      best_model_data$sigma.sq_SI[r] <- each_model$solution[2,7]
      best_model_data$alpha_d[r] <- each_model$solution[1,1]
      best_model_data$alpha_gc[r] <- each_model$solution[1,2]
      best_model_data$alpha_HD[r] <- each_model$solution[1,3]
      best_model_data$alpha_MI[r] <- each_model$solution[1,4]
      best_model_data$alpha_OM[r] <- each_model$solution[1,5]
      best_model_data$alpha_PK[r] <- each_model$solution[1,6]
      best_model_data$alpha_SI[r] <- each_model$solution[1,7]
      best_model_data$theta_d[r] <- each_model$theta[1,1]
      best_model_data$theta_gc[r] <- each_model$theta[2,1]
      best_model_data$theta_HD[r] <- each_model$theta[3,1]
      best_model_data$theta_MI[r] <- each_model$theta[4,1]
      best_model_data$theta_OM[r] <- each_model$theta[5,1]
      best_model_data$theta_PK[r] <- each_model$theta[6,1]
      best_model_data$theta_SI[r] <- each_model$theta[7,1]
      best_model_data$theta_d_se[r] <- each_model$theta[1,2]
      best_model_data$theta_d_se[r] <- each_model$theta[2,2]
      best_model_data$theta_d_se[r] <- each_model$theta[3,2]
      best_model_data$theta_d_se[r] <- each_model$theta[4,2]
      best_model_data$theta_d_se[r] <- each_model$theta[5,2]
      best_model_data$theta_d_se[r] <- each_model$theta[6,2]
      best_model_data$theta_d_se[r] <- each_model$theta[7,2]
      best_model_data$saddle[r] <- any(each_model$eigval < 0) #all should be over to 0 if well-fit healthy OU model/good run
      
      q <- q + 1
      r <- r + 1
    }
    i <- i + 1
  }
  j <- j + 1
}


OUM_diet<-best_model_data



###########################################################
#   Pulling out average theta values from OUwie analysis
###########################################################

trait_vector <- c("gape", "Jaw.Protrusion", "Adductor.Mass", "head.height.head.length", 
                  "Premax.Length", "Maxilla.Length", "Lower.Jaw.Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")


theta_results <- setNames(data.frame(matrix(ncol = 8, nrow = 13)), c("trait", "MIB", "GC", "MIS", "HD", "SI", "OM", "PK")) 


for(i in 1:length(trait_vector)){
  theta_results$trait[i]<-trait_vector[i]
  OUM<-filter(OUM_diet, trait == trait_vector[i], model=="OUM")
  theta_results$MIB[i]<- mean(OUM$theta_d)
  theta_results$GC[i]<- mean(OUM$theta_gc)
  theta_results$HD[i]<- mean(OUM$theta_HD)
  theta_results$MIS[i]<- mean(OUM$theta_MI)
  theta_results$OM[i]<- mean(OUM$theta_OM)
  theta_results$PK[i]<- mean(OUM$theta_PK)
  theta_results$SI[i]<- mean(OUM$theta_SI)
}



###############################################
#       Calculating mean trait values
###############################################

trait_troph<-dat.rat.na.log[,c(1,3, 5:17)]

trait_vector <- c("gape", "Jaw Protrusion", "Adductor Mass", "head.height.head.length", 
                  "Premax Length", "Maxilla Length", "Lower Jaw Lever", "close.out", "open.out", "x.low.jaw", "y.low.jaw",
                  "x.palatine", "y.palatine")

aovs_table <- setNames(data.frame(matrix(ncol = 8, nrow = 13)), c("trait", "GC", "MIB", "HD", "MIS", "OM", "PK", "SI" 
))

anovas <- for(i in 1:13){
  print(paste("Now starting trait:",trait_vector[i]))
  #morph_anova <- procD.pgls(trait ~ trop.func, phy = tree.final, iter = 9999)
  #aovs_list[i,2] <- list(morph_anova)
  aovs_table$trait[i] <- trait_vector[i]
  aovs_table$GC[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "GC")]) 
  aovs_table$MIB[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "D")])
  aovs_table$HD[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "HD")])
  aovs_table$MIS[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "MI")])
  aovs_table$OM[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "OM")])
  aovs_table$PK[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "PK")])
  aovs_table$SI[i] <- mean(trait_troph[,i+2][which(trait_troph$Diet == "SI")])
  #aovs_table$p.val[i] <- morph_anova$aov.table$`Pr(>F)`[[1]]
}




########################################################################
#       Replacing thetas that could not converge with mean trait values
#########################################################################

theta_results_mean<-theta_results
theta_results_mean[2,]<-aovs_table[2,]
theta_results_mean[4,]<-aovs_table[4,]
theta_results_mean[10,]<-aovs_table[10,]
theta_results_mean[11,]<-aovs_table[11,]
theta_results_mean[12,]<-aovs_table[12,]


##########################################
#      Transposing dataframe for PCA
##########################################

theta_pca = setNames(data.frame(names(theta_results_mean)[-1], 
                             t(unname(theta_results_mean[,-1]))), 
                  c('Diet', theta_results_mean[,1]))



# bit by bit
# pull out the column names, but skip the first one (we don't need it)
col_names = names(theta_results_mean)[-1]
# trim the data frame to remove the first column and remove the column headers
theta_results_simp = unname(theta_results_mean[,-1])
# transpose the numbers only
theta_results_t = t(theta_results_simp)
# combine the original column names and the transposed data frame into a new 
# data frame with the column names as the first column
theta_results_pca = data.frame(col_names, theta_results_t)
# assign new column names based on the original data frames first row
theta_results_pca = setNames(theta_results_pca, c('Diet', theta_results_mean[,1]))



####################################
#     Figure 5A: PCA
###################################


pca.theta<-prcomp(theta_results_pca[,2:13], scale=TRUE)

theta_results_pca$Color<-theta_results_pca$Diet
theta_results_pca$Color[theta_results_pca$Color == "GC"]<-"#F9A21E"
theta_results_pca$Color[theta_results_pca$Color == "OM"]<-"#F06A95"
theta_results_pca$Color[theta_results_pca$Color == "HD"]<-"#11619B"
theta_results_pca$Color[theta_results_pca$Color == "MIS"]<-"#27B79B"
theta_results_pca$Color[theta_results_pca$Color == "SI"]<-"#47ADE1"
theta_results_pca$Color[theta_results_pca$Color == "PK"]<-"#9B519E" 
theta_results_pca$Color[theta_results_pca$Color == "MIB"]<-"#E21C4D" 

plot(pca.theta$x[,1], pca.theta$x[,2], pch=16, col=theta_results_pca$Color,
     xlab="PC1", ylab="PC2", main="PC1 vs. PC2")

####################
#    Figure 5B: Dot plot
####################

dotchart(theta_results_mean[,2], pch = 21, labels = theta_results_mean$trait, bg = "#E21C4D", xlab="Trait Value",
         pt.cex = 1.5, xlim=c(-4,1))
points(theta_results_mean[,3], 1:nrow(dat.dot), col = "#F9A21E", pch = 19, cex = 1.5)
points(theta_results_mean[,4], 1:nrow(dat.dot), col = "#27B79B", pch = 19, cex = 1.5)
points(theta_results_mean[,5], 1:nrow(dat.dot), col = "#11619B", pch = 19, cex = 1.5)
points(theta_results_mean[,6], 1:nrow(dat.dot), col = "#47ADE1", pch = 19, cex = 1.5)
points(theta_results_mean[,7], 1:nrow(dat.dot), col = "#F06A95", pch = 19, cex = 1.5)
points(theta_results_mean[,8], 1:nrow(dat.dot), col = "#9B519E", pch = 19, cex = 1.5)








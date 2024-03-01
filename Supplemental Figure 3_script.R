##################################
#  Spearman rank correlation analysis for 
dat.rank.rate<-read.csv("Trait_rates_diet_rank order.csv")

cor.rat<-cor(dat.rank.rate[,2:14], method="spearman")

pdf(file="Heatmap_correlations_rate_half.pdf", width=10, height=10)
corrplot(cor.rat, methods="square", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, is.corr=FALSE)
dev.off()


dat.rank.disp<-read.csv("Trait_disparity_diet_rotated.csv")

cor.disp<-cor(dat.rank.disp[,2:14], method="spearman")

pdf(file="Heatmap_correlations_disp_half.pdf", width=10, height=10)
corrplot(cor.disp, methods="square", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, is.corr=FALSE)
dev.off()


#################################
#       Data Prep
#################################


dat.rank<-read.csv("Theta diet rank values.csv")

dat.rank$Color<-dat.rank$trait
dat.rank$Color[dat.rank$Color == "Benthic.invertivores"]<-"#E21C4D" 
dat.rank$Color[dat.rank$Color == "General.carnivore"]<-"#F9A21E"
dat.rank$Color[dat.rank$Color == "Herbivore.detritivore"]<-"#11619B"
dat.rank$Color[dat.rank$Color == "Mobile.invertivore"]<-"#27B79B"
dat.rank$Color[dat.rank$Color == "Omnivore"]<-"#F06A95"
dat.rank$Color[dat.rank$Color == "Planktivore"]<-"#9B519E" 
dat.rank$Color[dat.rank$Color == "Sessile.invertivore"]<-"#47ADE1"


#################################
# Perform Spearman correlation
#################################

require(corrplot)

cor.trait<-cor(dat.rank[,2:9],method="spearman")

pdf(file="Heatmap_correlations_optima_half.pdf", width=10, height=10)
corrplot(cor.trait, methods="square", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, is.corr=FALSE)
dev.off()


#############################################
#         Bivariate plots
#############################################

# Mouth gape vs Premaxilla length
corr <- cor.test(x=dat$gape, y=dat$Premax.Length, method = 'spearman')
pdf(file="Correlation_plot_gape_premax.pdf", width=8, height=8)
plot(dat.rank$gape, dat.rank$Premax.Length, col=dat.rank$Color, pch=16)
abline(a = 1, b = corr$estimate)
dev.off()

# Mouth gape vs Lower Jaw length
corr <- cor.test(x=dat$gape, y=dat$Lower.Jaw.Lever, method = 'spearman')
pdf(file="Correlation_plot_gape_lower jaw.pdf", width=8, height=8)
plot(dat.rank$gape, dat.rank$Lower.Jaw.Lever, col=dat.rank$Color, pch=16)
abline(a = 1, b = corr$estimate)
dev.off()

# Mouth gape vs Closing mechanical advantage
corr <- cor.test(x=dat$gape, y=dat$close.out, method = 'spearman')
pdf(file="Correlation_plot_gape_close.pdf", width=8, height=8)
plot(dat.rank$gape, dat.rank$close.out, col=dat.rank$Color, pch=16)
abline(a = 7, b = corr$estimate)
dev.off()

# Maxilla Length vs Closing mechanical advantage
corr <- cor.test(x=dat$Maxilla.Length, y=dat$close.out, method = 'spearman')
pdf(file="Correlation_plot_max_close.pdf", width=8, height=8)
plot(dat.rank$Maxilla.Length, dat.rank$close.out, col=dat.rank$Color, pch=16)
abline(a = 7, b = corr$estimate)
dev.off()
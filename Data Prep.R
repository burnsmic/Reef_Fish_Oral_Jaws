#####################################
#      Read in Raw Data
#####################################

dat<-read.csv("Combined_dataset_log.csv")

############################
#    Species averages
# 
# Lines 13-23: Is used to create a dataframe where the species average is calculated for each morphological trait
############################

dat.species<-aggregate(x=dat[c("head_length", "head_height", "Dentigerous_premax", "maxilla", "out_lever_lower_jaw", "closing_in_lever", "opening_in_lever",
                               "Xpos_lower_jaw_joint", "Ypos_lower_jaw_joint_inv", "Xpos_maxilla_joint_with_palatine", "Ypos_maxilla_joint_with_palatine_inv", "SL",
                               "Gape", "Protrusion_diff_abs", "Body_Mass", "Body_Depth", "Body_Width", "Head_Width", "Adductor_Mass", "log_SL" )], by=dat[c("Tree_Name")], FUN=mean, na.rm=TRUE)

dat.eco<-aggregate(x=dat["Trophic"], by=dat["Tree_Name"], FUN=unique)
dat.troph<-aggregate(x=dat["TG"], by=dat["Tree_Name"], FUN=unique)
dat.function<-aggregate(x=dat["Function"], by=dat["Tree_Name"], FUN=unique)

dat.species$Trophic<-dat.eco$Trophic
dat.species$Diet<-dat.troph$TG
dat.species$Function<-dat.function$Function

#############################
# Lines 32-145: Code used to take ratios of traits to create size corrected variables
#############################


########################
#   Gape/head_length Species
########################

gape.head.length<-as.data.frame(dat.species$Gape/dat.species$head_length)
gape.head.length$Species<-dat.species$Tree_Name
gape.head.length.na<-na.omit(gape.head.length)

##############################
#    Protrusion/head_length
##############################

prod.head.length <-as.data.frame(dat.species$Protrusion_diff_abs/dat.species$head_length)
prod.head.length$Species<-dat.species$Tree_Name
prod.head.length.na<-na.omit(prod.head.length)

#################################
#     Max Body Depth/SL
#################################

body.depth.sl <-as.data.frame(dat.species$Body_Depth/dat.species$SL)
body.depth.sl$Species<-dat.species$Tree_Name
body.depth.sl.na<-na.omit(body.depth.sl)

#################################
#    Max Body Width/SL
##################################
body.width.sl <-as.data.frame(dat.species$Body_Width/dat.species$SL)
body.width.sl$Species<-dat.species$Tree_Name
body.width.sl.na<-na.omit(body.width.sl)

#################################
#    Adductor Mass/Body Mass
##################################

add.mass.body.mass<-as.data.frame(dat.species$Adductor_Mass/dat.species$Body_Mass)
add.mass.body.mass$Species<-dat.species$Tree_Name
add.mass.body.mass.na<-na.omit(add.mass.body.mass)

###################################
#      Head_height/head_length
###################################

head.height.head.length<-as.data.frame(dat.species$head_height/dat.species$head_length)
head.height.head.length$Species<-dat.species$Tree_Name
head.height.head.length.na<-na.omit(head.height.head.length)

#######################################
#      Dentigerous_premax/head_length
#########################################

premax.head.length<-as.data.frame(dat.species$Dentigerous_premax/dat.species$head_length)
premax.head.length$Species<-dat.species$Tree_Name
premax.head.length.na<-na.omit(premax.head.length)

#############################################
#       Maxilla/head_length
#############################################

max.head.length<-as.data.frame(dat.species$maxilla/dat.species$head_length)
max.head.length$Species<-dat.species$Tree_Name
max.head.length.na<-na.omit(max.head.length)

############################################
#       Lower Jaw out lever/Head Length
############################################

low.jaw.lever<-as.data.frame(dat.species$out_lever_lower_jaw/dat.species$head_length)
low.jaw.lever$Species<-dat.species$Tree_Name
low.jaw.lever.na<-na.omit(low.jaw.lever)

############################################
#        Close_in_lever/out_lever_lower_jaw
#############################################

close.out<-as.data.frame(dat.species$closing_in_lever/dat.species$out_lever_lower_jaw)
close.out$Species<-dat.species$Tree_Name
close.out.na<-na.omit(close.out)

############################################
#        Opening_in_lever/out_lever_lower_jaw
#############################################

open.out<-as.data.frame(dat.species$opening_in_lever/dat.species$out_lever_lower_jaw)
open.out$Species<-dat.species$Tree_Name
open.out.na<-na.omit(open.out)

############################################
#        Xpos_lower_jaw_joint/head_length
#############################################

x.low.jaw.head.length<-as.data.frame(dat.species$Xpos_lower_jaw_joint/dat.species$head_length +1) # need to add one to make the values positive for log transformation
x.low.jaw.head.length$Species<-dat.species$Tree_Name
x.low.jaw.head.length.na<-na.omit(x.low.jaw.head.length)

############################################
#        Ypos_lower_jaw_joint/head_length
#############################################

y.low.jaw.head.length<-as.data.frame(dat.species$Ypos_lower_jaw_joint_inv/dat.species$head_length +1) # need to add one to make the values positive for log transformation
y.low.jaw.head.length$Species<-dat.species$Tree_Name
y.low.jaw.head.length.na<-na.omit(y.low.jaw.head.length)

############################################
# Xpos_maxilla_joint_with_palatine/head_length
#############################################

x.palatine.head.length<-as.data.frame(dat.species$Xpos_maxilla_joint_with_palatine/dat.species$head_length + 1)
x.palatine.head.length$Species<-dat.species$Tree_Name
x.palatine.head.length.na<-na.omit(x.palatine.head.length)

############################################
# Ypos_maxilla_joint_with_palatine/head_length
#############################################

y.palatine.head.length<-as.data.frame(dat.species$Ypos_maxilla_joint_with_palatine_inv/dat.species$head_length +1)
y.palatine.head.length$Species<-dat.species$Tree_Name
y.palatine.head.length.na<-na.omit(y.palatine.head.length)


###################
# Lines 153-176: Creates a dataframe with log ratios of traits, species names, functional feeding mode, and diet.
###############

require(tidyverse)
dat.list<-list(gape.head.length,prod.head.length,  add.mass.body.mass, head.height.head.length,
               premax.head.length, max.head.length, low.jaw.lever, close.out, open.out, x.low.jaw.head.length, y.low.jaw.head.length,
               x.palatine.head.length, y.palatine.head.length)

dat.rat<-dat.list %>% reduce(full_join, by='Species')
dat.rat$Trophic<-dat.species$Trophic
dat.rat$Diet<-dat.species$Diet
dat.rat$Function<-dat.species$Function
dat.rat.na<-na.omit(dat.rat)

dat.rat.na$log.Gape.head.length<-log(dat.rat.na$`dat.species$Gape/dat.species$head_length`)
dat.rat.na$log.prod.head.length<-log(dat.rat.na$`dat.species$Protrusion_diff_abs/dat.species$head_length`)
dat.rat.na$log.add.mass.body.mass<-log(dat.rat.na$`dat.species$Adductor_Mass/dat.species$Body_Mass`)
dat.rat.na$log.head.height.head.length<-log(dat.rat.na$`dat.species$head_height/dat.species$head_length`)
dat.rat.na$log.premax.head.length<-log(dat.rat.na$`dat.species$Dentigerous_premax/dat.species$head_length`)
dat.rat.na$log.max.head.length<-log(dat.rat.na$`dat.species$maxilla/dat.species$head_length`)
dat.rat.na$log.low.jaw.lever<-log(dat.rat.na$`dat.species$out_lever_lower_jaw/dat.species$head_length`)
dat.rat.na$log.close.out<-log(dat.rat.na$`dat.species$closing_in_lever/dat.species$out_lever_lower_jaw`)
dat.rat.na$log.open.out<-log(dat.rat.na$`dat.species$opening_in_lever/dat.species$out_lever_lower_jaw`)
dat.rat.na$log.x.low.jaw.head.length<-log(dat.rat.na$`dat.species$Xpos_lower_jaw_joint/dat.species$head_length`)
dat.rat.na$log.y.low.jaw.head.length<-log(dat.rat.na$`dat.species$Ypos_lower_jaw_joint_inv/dat.species$head_length`)
dat.rat.na$log.x.palatine.head.length<-log(dat.rat.na$`dat.species$Xpos_maxilla_joint_with_palatine/dat.species$head_length`)
dat.rat.na$log.y.palatine.head.length<-log(dat.rat.na$`dat.species$Ypos_maxilla_joint_with_palatine_inv/dat.species$head_length`)

dat.rat.na.clean<-dat.rat.na[,c(2,15:17,18:30)]

#######################################################
#       Read in phylogeny and match tips to tree
#######################################################
dat.rat.na<-dat.rat.na.clean

require(fishtree)
require(geiger)

tree<-fishtree_phylogeny()

dd<-name.check(tree, dat.rat.na, data.names=dat.rat.na$Species)
tree.final<-drop.tip(tree,dd$tree_not_data)
#save(tree.final, file="Reef_fish_tree.RData")

###########################################################


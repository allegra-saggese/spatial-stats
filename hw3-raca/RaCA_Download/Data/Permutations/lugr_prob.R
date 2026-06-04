rm(list = ls() ) #clear environment

############################################################
#Calculation of LUGR pixel count weights
#########################################################

#load packages
library(plyr) 
#########################################################
#folder locations of source and output data - change to local path that contains input data
#sub folder paths are then updated


#raca = "F:/RaCA-newphase/RaCA_Cloudvault/"
raca = 'C:/Users/skyew/Documents/data/'

wt <- read.csv(paste0(raca, "LUGR_pixelcount.csv"))
wt$LU <- as.factor(substr(wt$LUGR,1,1))
names(wt)

#sum total pixel counts for the entire area, each MO and each LU
#while inclusion propobabilities are not precisely equal within LUGR, we are considering this sample draw random

wt$TOTAL <- sum(as.numeric(wt$group_count))

LU_TOTAL <- ddply(wt, "LU", summarise, 
                       LU_TOTAL= sum(as.numeric(group_count))
                       )

MO_TOTAL <- ddply(wt, "MO", summarise, 
                       MO_TOTAL = sum(as.numeric(group_count))
)

LUbyMO_TOTAL <- ddply(wt, .(LU, MO), summarise, 
                  LUbyMO_TOTAL = sum(as.numeric(group_count))
)


#join to pixel count file
Wt <- join(wt, LU_TOTAL, by='LU')
wT <- join(Wt, MO_TOTAL, by='MO')
WT <- join(wT, LUbyMO_TOTAL, by=c('MO', 'LU'))


WT$prob.oa <- as.numeric(WT$group_count)/WT$TOTAL
WT$prob.LU <- as.numeric(WT$group_count)/WT$LU_TOTAL
WT$prob.MO <- as.numeric(WT$group_count)/WT$MO_TOTAL
WT$prob.LUbyMO  <- as.numeric(WT$group_count)/WT$LUbyMO_TOTAL

names(WT)

write.csv(WT, paste0(raca, "LUGR_prob.csv"))


#check that probabilites sum to zero
sum(WT$prob.oa)

check_sum <- ddply(WT, .(LU), summarise, 
                   l.sum = sum(prob.LU)
                   )

check_sum <- ddply(WT, .(MO), summarise, 
                   l.sum = sum(prob.LU)
)

check_sum <- ddply(WT, .(MO, LU), summarise, 
                   l.sum = sum(prob.LUbyMO)
)

summary(WT$prob.oa)





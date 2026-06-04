
###################################################
#LUGR and LU summary for mapping##
####################################################
#weights used for regional summary
##################################
install.packages(c("soilDB", "shrpshootR", "aqp", "lattice", "plyr", "reshape2", "ggplot2", "lme4", "stringr", "Hmisc", "SDMTools"))

rm(list = ls() ) #clear environment

#load packages
library(soilDB)
library(sharpshootR)
library(aqp)
library(lattice)
library(plyr)
library(reshape2)
library(ggplot2)
library(lme4)
library(stringr)
library(Hmisc)
library(SDMTools)

#folder locations of source and output data - change to local path that contains input data

raca = "F:/RaCA-newphase/RaCA_Cloudvault/RaCA_Download/Data/"
#raca = 'C:/Users/skyew/Documents/data/'


LUGR_wt <- read.csv(paste0(raca, "LUGR_pixelcount.csv"))

SOCpedons <- read.csv(paste0(raca, "RaCA_SOC_pedons.csv"))


#create subfolders to save output
summ = paste0(raca, "Summary/")
ifelse(!dir.exists(file.path(summ)), dir.create(summ), FALSE)

wt_avg_folder = paste0(summ, "wt_avg/")
ifelse(!dir.exists(file.path(wt_avg_folder)), dir.create(wt_avg_folder), FALSE)


#create subfolders to save output
####DOEs not save directly to cloudvault#####
LUGR = "F:/RaCA-newphase/RaCA_Cloudvault/LUGR_summaries/"
ifelse(!dir.exists(file.path(LUGR)), dir.create(LUGR), FALSE)




#combine land use - land cover class and  soil group - LUGR

SOCpedons$LUGR <- as.factor(paste0(as.character(SOCpedons$LU),substr(SOCpedons$rcasiteid,2,5)))


#Count by MO and LU
count <- ddply(SOCpedons, .(MO,LU), summarise, N = sum(!is.na(upedon)))
Count.t <- reshape(count, idvar="MO", timevar = "LU", direction = "wide" )
write.table(Count.t, paste0(summ, "RaCA_Table_count.csv"), sep=",", row.names=F)


#LUGR Means

LUGR <- recast(SOCpedons, LUGR + MO + LU  ~ variable, mean
               , id.var = c("LUGR", "MO", "LU"), measure.var = 9:11, na.rm=T)

LUGRsd <- recast(SOCpedons, LUGR + MO + LU  ~ variable, sd
               , id.var = c("LUGR", "MO", "LU"), measure.var = 9:11, na.rm=T)

names(LUGRsd)[4] <- "SOCsd5"
names(LUGRsd)[5] <- "SOCsd30"
names(LUGRsd)[6] <- "SOCsd100"

LUGRj <- cbind(LUGR,LUGRsd)

LUbyMO <- recast(SOCpedons,  MO + LU ~ variable, mean
                 , id.var = c( "MO", "LU"), measure.var = 9:11, na.rm=T)

LUbyMOsd <- recast(SOCpedons,  MO + LU ~ variable, sd
                 , id.var = c( "MO", "LU"), measure.var = 9:11, na.rm=T)
names(LUbyMOsd)[3] <- "SOCsd5"
names(LUbyMOsd)[4] <- "SOCsd30"
names(LUbyMOsd)[5] <- "SOCsd100"

LUbyMO$LUGR <- paste0(as.character(LUbyMO$LU), sprintf("%02d", LUbyMO$MO), "00")
LUbyMOsd$LUGR <- paste0(as.character(LUbyMOsd$LU), sprintf("%02d", LUbyMOsd$MO), "00")

LUbyMOj <- cbind(LUbyMO, LUbyMOsd)

MOGR <- recast(SOCpedons,  MO + MOGr ~ variable, mean
                 , id.var = c("MOGr", "MO"), measure.var = 9:11, na.rm=T)
MOGRsd <- recast(SOCpedons,  MO + MOGr ~ variable, sd
               , id.var = c("MOGr", "MO"), measure.var = 9:11, na.rm=T)
names(MOGRsd)[3] <- "SOCsd5"
names(MOGRsd)[4] <- "SOCsd30"
names(MOGRsd)[5] <- "SOCsd100"

MOGR$LUGR <- paste(sprintf("%04d", MOGR$MOGr))
MOGR$LU <- "Unknown"
MOGRsd$LUGR <- paste(sprintf("%04d", MOGR$MOGr))
MOGRsd$LU <- "Unknown"

MOGRj <- cbind2(MOGR, MOGRsd)

MO <- recast(SOCpedons,  MO ~ variable, mean
             , id.var = c( "MO"), measure.var = 9:11, na.rm=T)

MOsd <- recast(SOCpedons,  MO ~ variable, sd
             , id.var = c( "MO"), measure.var = 9:11, na.rm=T)
names(MOsd)[2] <- "SOCsd5"
names(MOsd)[3] <- "SOCsd30"
names(MOsd)[4] <- "SOCsd100"

MO$LUGR <- paste0(sprintf("%02d", MO$MO), "00")
MO$LU <- "Unknown"

MOsd$LUGR <- paste0(sprintf("%02d", MO$MO), "00")
MOsd$LU <- "Unknown"


MOj <- cbind2(MO, MOsd)


LUGR_join <- rbind.fill(LUGRj, LUbyMOj, MOj, MOGRj)


write.table(LUGR_join, paste0(wt_avg_folder , "LUGR_join.csv"), sep=",", row.names=F)


write.table(LUGR_join, paste0(summ, "RaCA_LUGR_stocks.csv"), sep=",", row.names=F)


#######################################################################
#CALCULATE LUGR means
#transform each pedon value to better approximate normality
#then back transform for export

#LU and MO averages are used to fill pixels with incomplete LUGR information
##################################################


#previously loaded pedon stocks
SOCpedons <- SOCpedons[SOCpedons$MOGrLU!="",]


#calculate natural log of each depth stock of interest
lnLUGr5 <- aggregate(log(SOCstock5)~MO+LU+LUGR+MOGrLU, data=SOCpedons, FUN=mean)
names(lnLUGr5)[5] <- "lnLUGr5"
lnLUGr30 <- aggregate(log(SOCstock30)~LUGR, data=SOCpedons, FUN=mean)
names(lnLUGr30)[2] <- "lnLUGr30"
lnLUGr100 <- aggregate(log(SOCstock100)~LUGR, data=SOCpedons, FUN=mean)
names(lnLUGr100)[2] <- "lnLUGr100"

LUGrsd5 <- aggregate(log(SOCstock5)~LUGR, data=SOCpedons, sd)
names(LUGrsd5)[2] <- "lnLUGrsd5"
LUGrsd30 <- aggregate(log(SOCstock30)~LUGR, data=SOCpedons, sd)
names(LUGrsd30)[2] <- "lnLUGrsd30"
LUGrsd100 <- aggregate(log(SOCstock100)~LUGR, data=SOCpedons, sd)
names(LUGrsd100)[2] <- "lnLUGrsd100"


LUGrse5 <- aggregate(log(SOCstock5)~LUGR, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUGrse5)[2] <- "lnLUGrse5"
LUGrse30 <- aggregate(log(SOCstock30)~LUGR, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUGrse30)[2] <- "lnLUGrse30"
LUGrse100 <- aggregate(log(SOCstock100)~LUGR, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUGrse100)[2] <- "lnLUGrse100"

#bind columns of mean and se for each depth
LUGR5 <- cbind(lnLUGr5 , LUGrsd5[2], LUGrse5[2])
LUGR30 <- cbind(lnLUGr30, LUGrsd30[2], LUGrse30[2])
LUGR100 <- cbind(lnLUGr100,LUGrsd100[2], LUGrse100[2])

LUGR530<- join(LUGR5, LUGR30, by = "LUGR")
LUGR <- join(LUGR530, LUGR100, by = "LUGR")

#take out replicate columns and rename

names(LUGR)<- c("MO", "LU", "LUGR", "MOGRLU", "lnSOC5",  "lnSOCsd5","lnSOCse5", "lnSOC30", "lnSOCsd30", "lnSOCse30", "lnSOC100","lnSOCsd100","lnSOCse100")


LUGR$lnCIl5 <- with(LUGR, lnSOC5-lnSOCse5*1.96)
LUGR$lnCIu5 <- with(LUGR, lnSOC5+lnSOCse5*1.96)

LUGR$lnCIl30 <- with(LUGR, lnSOC30-lnSOCse30*1.96)
LUGR$lnCIu30 <- with(LUGR, lnSOC30+lnSOCse30*1.96)

LUGR$lnCIl100 <- with(LUGR, lnSOC100-lnSOCse100*1.96)
LUGR$lnCIu100 <- with(LUGR, lnSOC100+lnSOCse100*1.96)

write.table(LUGR, paste0(wt_avg_folder,"LUGR_ln.csv"), sep=",", row.names=F)

bt_LUGR <- cbind(LUGR[1:4], exp(LUGR[5:length(LUGR)]))
names(bt_LUGR)[5:length(bt_LUGR)] <- gsub('ln', 'bt_', names(bt_LUGR[5:length(bt_LUGR)]))

write.table(bt_LUGR, paste0(summ,"bt_LUGR_wCL.csv"), sep=",", row.names=F)

################################
#For MO and LU
#calculate natural log of each depth stock of interest
lnLUMO5 <- aggregate(log(SOCstock5)~MO+LU, data=SOCpedons, FUN=mean)
names(lnLUMO5)[3] <- "lnLUMO5"
lnLUMO30 <- aggregate(log(SOCstock30)~MO +LU, data=SOCpedons, FUN=mean)
names(lnLUMO30)[3] <- "lnLUMO30"
lnLUMO100 <- aggregate(log(SOCstock100)~MO + LU, data=SOCpedons, FUN=mean)
names(lnLUMO100)[3] <- "lnLUMO100"

LUMOsd5 <- aggregate(log(SOCstock5)~MO + LU, data=SOCpedons, sd)
names(LUMOsd5)[3] <- "lnLUMOsd5"
LUMOsd30 <- aggregate(log(SOCstock30)~MO + LU, data=SOCpedons, sd)
names(LUMOsd30)[3] <- "lnLUMOsd30"
LUMOsd100 <- aggregate(log(SOCstock100)~MO + LU, data=SOCpedons, sd)
names(LUMOsd100)[3] <- "lnLUMOsd100"


LUMOse5 <- aggregate(log(SOCstock5)~LU + MO, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUMOse5)[3] <- "lnLUMOse5"
LUMOse30 <- aggregate(log(SOCstock30)~LU + MO, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUMOse30)[3] <- "lnLUMOse30"
LUMOse100 <- aggregate(log(SOCstock100)~MO + LU, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(LUMOse100)[3] <- "lnLUMOse100"

#bind columns of mean, sd, and se for each depth
LUMO5 <- cbind(lnLUMO5 , LUMOsd5[3], LUMOse5[3])
LUMO5$LUGR <- paste0(as.character(LUMO5$LU), sprintf("%02d", LUMO5$MO), "00")

LUMO30 <- cbind(lnLUMO30, LUMOsd30[3], LUMOse30[3])
LUMO30$LUGR <- paste0(as.character(LUMO30$LU), sprintf("%02d", LUMO30$MO), "00")

LUMO100 <- cbind(lnLUMO100,LUMOsd100[3], LUMOse100[3])
LUMO100$LUGR <- paste0(as.character(LUMO30$LU), sprintf("%02d", LUMO100$MO), "00")

LUMO530<- join(LUMO5, LUMO30[,3:6], by = "LUGR")
LUMO <- join(LUMO530, LUMO100[,3:6], by = "LUGR")

#take out replicate columns and rename
names(LUMO)<- c("MO", "LU",  "lnSOC5",  "lnSOCsd5","lnSOCse5", "LUGR", "lnSOC30", "lnSOCsd30", "lnSOCse30", "lnSOC100","lnSOCsd100","lnSOCse100")

LUMO$lnCIl5 <- with(LUMO, lnSOC5-lnSOCse5*1.96)
LUMO$lnCIu5 <- with(LUMO, lnSOC5+lnSOCse5*1.96)

LUMO$lnCIl30 <- with(LUMO, lnSOC30-lnSOCse30*1.96)
LUMO$lnCIu30 <- with(LUMO, lnSOC30+lnSOCse30*1.96)

LUMO$lnCIl100 <- with(LUMO, lnSOC100-lnSOCse100*1.96)
LUMO$lnCIu100 <- with(LUMO, lnSOC100+lnSOCse100*1.96)


################################
#For MO
#calculate natural log of each depth stock of interest
lnMO5 <- aggregate(log(SOCstock5)~MO, data=SOCpedons, FUN=mean)
names(lnMO5)[2] <- "lnMO5"
lnMO30 <- aggregate(log(SOCstock30)~MO, data=SOCpedons, FUN=mean)
names(lnMO30)[2] <- "lnMO30"
lnMO100 <- aggregate(log(SOCstock100)~MO, data=SOCpedons, FUN=mean)
names(lnMO100)[2] <- "lnMO100"


MOsd5 <- aggregate(log(SOCstock5)~MO, data=SOCpedons, sd)
names(MOsd5)[2] <- "lnMOsd5"
MOsd30 <- aggregate(log(SOCstock30)~MO, data=SOCpedons, sd)
names(MOsd30)[2] <- "lnMOsd30"
MOsd100 <- aggregate(log(SOCstock100)~MO , data=SOCpedons, sd)
names(MOsd100)[2] <- "lnMOsd100"


MOse5 <- aggregate(log(SOCstock5)~MO, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOse5)[2] <- "lnMOse5"
MOse30 <- aggregate(log(SOCstock30)~MO, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOse30)[2] <- "lnMOse30"
MOse100 <- aggregate(log(SOCstock100)~MO, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOse100)[] <- "lnMOse100"

#bind columns of mean, sd, and se for each depth
MO5 <- cbind(lnMO5 , MOsd5[2], MOse5[2])
MO5$LUGR <- paste0(sprintf("%02d", MO5$MO), "00")
MO5$LU <- "Unknown"

MO30 <- cbind(lnMO30, MOsd30[2], MOse30[2])
MO30$LUGR <- paste0(sprintf("%02d", MO30$MO), "00")
MO30$LU <- "Unknown"

MO100 <- cbind(lnMO100,MOsd100[2], MOse100[2])
MO100$LUGR <- paste0(sprintf("%02d", MO100$MO), "00")
MO100$LU <- "Unknown"


MO530<- join(MO5, MO30[,2:5], by = "LUGR")
MO <- join(MO530, MO100[,2:5], by = "LUGR")

#take out replicate columns and rename
names(MO)<- c("MO",  "lnSOC5",  "lnSOCsd5", "lnSOCse5","LUGR", "LU", "lnSOC30", "lnSOCsd30", "lnSOCse30"
              , "lnSOC100","lnSOCsd100","lnSOCse100")

MO$lnCIl5 <- with(MO, lnSOC5-lnSOCse5*1.96)
MO$lnCIu5 <- with(MO, lnSOC5+lnSOCse5*1.96)

MO$lnCIl30 <- with(MO, lnSOC30-lnSOCse30*1.96)
MO$lnCIu30 <- with(MO, lnSOC30+lnSOCse30*1.96)

MO$lnCIl100 <- with(MO, lnSOC100-lnSOCse100*1.96)
MO$lnCIu100 <- with(MO, lnSOC100+lnSOCse100*1.96)

#####################
#calculate average by MO and Group
##################
#calculate natural log of each depth stock of interest

lnMOGr5 <- aggregate(log(SOCstock5)~MOGr, data=SOCpedons, FUN=mean)
names(lnMOGr5)[2] <- "lnMOGr5"
lnMOGr30 <- aggregate(log(SOCstock30)~MOGr, data=SOCpedons, FUN=mean)
names(lnMOGr30)[2] <- "lnMOGr30"
lnMOGr100 <- aggregate(log(SOCstock100)~MOGr, data=SOCpedons, FUN=mean)
names(lnMOGr100)[2] <- "lnMOGr100"

MOGrsd5 <- aggregate(log(SOCstock5)~MOGr, data=SOCpedons, sd)
names(MOGrsd5)[2] <- "lnMOGrsd5"
MOGrsd30 <- aggregate(log(SOCstock30)~MOGr, data=SOCpedons, sd)
names(MOGrsd30)[2] <- "lnMOGrsd30"
MOGrsd100 <- aggregate(log(SOCstock100)~MOGr , data=SOCpedons, sd)
names(MOGrsd100)[2] <- "lnMOGrsd100"


MOGrse5 <- aggregate(log(SOCstock5)~MOGr, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOGrse5)[2] <- "lnMOGrse5"
MOGrse30 <- aggregate(log(SOCstock30)~MOGr, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOGrse30)[2] <- "lnMOGrse30"
MOGrse100 <- aggregate(log(SOCstock100)~MOGr, data=SOCpedons, function(x) c(SE = sd(x)/sqrt(length(x))))
names(MOGrse100)[] <- "lnMOGrse100"

#bind columns of mean, sd, and se for each depth
MOGr5 <- cbind(lnMOGr5 , MOGrsd5[2], MOGrse5[2])
MOGr5$LUGR <- paste0(sprintf("%04d", MOGr5$MOGr))
MOGr5$LU <- "Unknown"

MOGr30 <- cbind(lnMOGr30, MOGrsd30[2], MOGrse30[2])
MOGr30$LUGR <- paste0(sprintf("%04d", MOGr30$MOGr))
MOGr30$LU <- "Unknown"

MOGr100 <- cbind(lnMOGr100,MOGrsd100[2], MOGrse100[2])
MOGr100$LUGR <- paste0(sprintf("%04d", MOGr100$MOGr))
MOGr100$LU <- "Unknown"



MOGr530<- join(MOGr5, MOGr30[,2:5], by = "LUGR")
MOGr <- join(MOGr530, MOGr100[,2:5], by = "LUGR")

#take out replicate columns and rename
names(MOGr)<- c("MO",  "lnSOC5",  "lnSOCsd5", "lnSOCse5","LUGR", "LU", "lnSOC30", "lnSOCsd30", "lnSOCse30"
              , "lnSOC100","lnSOCsd100","lnSOCse100")

MOGr$lnCIl5 <- with(MOGr, lnSOC5-lnSOCse5*1.96)
MOGr$lnCIu5 <- with(MOGr, lnSOC5+lnSOCse5*1.96)

MOGr$lnCIl30 <- with(MOGr, lnSOC30-lnSOCse30*1.96)
MOGr$lnCIu30 <- with(MOGr, lnSOC30+lnSOCse30*1.96)

MOGr$lnCIl100 <- with(MOGr, lnSOC100-lnSOCse100*1.96)
MOGr$lnCIu100 <- with(MOGr, lnSOC100+lnSOCse100*1.96)


#Join all LUGR and replacements together 
################
ln_LUGRall_join <- rbind.fill(LUGR[,-(4)], LUMO, MO, MOGr)

#output
###################
write.table(ln_LUGRall_join, paste0(wt_avg_folder,"LUGRall_ln.csv"), sep=",", row.names=F)

bt_LUGRall <- cbind(ln_LUGRall_join[1:3], exp(ln_LUGRall_join [4:length(ln_LUGRall_join )]))

names(bt_LUGRall)[4:length(bt_LUGRall)] <- gsub('ln', 'bt_', names(bt_LUGRall[4:length(bt_LUGRall)]))


#for GIS mapping

write.table(bt_LUGRall, paste0(summ,"bt_RaCA_LUGRall_stocks.csv"), sep=",", row.names=F)


################################################################################################
#################################################################################################

#some unweighted means
NOTwt_avgLM <- ddply(SOCpedons, .(LU, as.factor(MO)), summarise, 
                   m5 = mean(log(SOCstock5), na.rm=T),
                   m30 = mean(log(SOCstock30),na.rm=T),
                   m100 = mean(log(SOCstock100), na.rm=T)
)
names(NOTwt_avgLM)[2] <- "MO"

NOTwt_avg <- ddply(SOCpedons, .(LU), summarise, 
                   m5 = mean(log(SOCstock5), na.rm=T),
                   m30 = mean(log(SOCstock30),na.rm=T),
                   m100 = mean(log(SOCstock100), na.rm=T)
)

NOTwt_avg$MO <- "Overall"

NOT <- rbind(NOTwt_avg, NOTwt_avgLM)

No_wt_avg <- cbind(NOT[,c(1,5)], exp(NOT[,2:4]))

#write.table(No_wt_avg, paste0(wt_avg_folder,"UNwt_avg.csv"), sep =",", row.names = F)

#############

#should I do some group means and sd from ln pedons??????????????????????????????

#########################################################
#Weights
#################################
# use pixel counts to weight LUGR averages for each MO 
#########################################################

#create new folder for new output
#if the folder already exsists it returns FALSE

wt <- join(LUGR_wt, LUGR, by ="MOGRLU")
wt$MO <- factor(wt$MO)

#weight sd

#calculate weighted average overall
ln_wt_avg <- ddply(wt, .(), summarise, 
                   wmn5 = weighted.mean(lnSOC5, group_count, na.rm=T),
                   wmn30 = weighted.mean(lnSOC30, group_count,na.rm=T),
                   wmn100 = weighted.mean(lnSOC100, group_count, na.rm=T),
                  wsd5 = sqrt(wtd.var(lnSOC5, w = as.numeric(group_count), na.rm=T)),
                   wsd30 = sqrt(wtd.var(lnSOC30, as.numeric(group_count), na.rm=T)),
                   wsd100 = sqrt(wtd.var(lnSOC100, as.numeric(group_count), na.rm=T))
                 )



ln_wt_avg[1,1] <- "ALL"
names(ln_wt_avg)[1] <- "Average For"


ln_wt_avg_MO <- ddply(wt, "MO", summarise,
                      N = sum(!is.na(lnSOC5)),    
                   wmn5 = weighted.mean(lnSOC5, group_count, na.rm=T),
                  wmn30 = weighted.mean(lnSOC30, group_count,na.rm=T),
                  wmn100 = weighted.mean(lnSOC100, group_count, na.rm=T),
                  t.wm.100 = wtd.mean(lnSOC100, as.numeric(group_count), na.rm=T),
                  wsd5 = sqrt(wtd.var(lnSOC5, as.numeric(group_count), na.rm=T)),
                  wsd30 = sqrt(wtd.var(lnSOC30, as.numeric(group_count), na.rm=T)),
                  wsd100 = sqrt(wtd.var(lnSOC100, as.numeric(group_count), na.rm=T)),
                  se100 = wsd100/sqrt(N),
                  CI100_l = wmn100 - 1.96*se100,
                  CI100_u = wmn100 +  1.96*se100              
                  )


 
ln_wt_avg_LU <- ddply(wt, "LU", summarise,
                      N = sum(!is.na(lnSOC5)),   
                      wmn5 = weighted.mean(lnSOC5, group_count, na.rm=T),
                   wmn30 = weighted.mean(lnSOC30, group_count,na.rm=T),
                   wmn100 = weighted.mean(lnSOC100, group_count, na.rm=T),
                   wsd5 = sqrt(wtd.var(lnSOC5, as.numeric(group_count), na.rm=T)),
                   wsd30 = sqrt(wtd.var(lnSOC30, as.numeric(group_count), na.rm=T)),
                   wsd100 = sqrt(wtd.var(lnSOC100, as.numeric(group_count), na.rm=T)),
                   se100 = wsd100/sqrt(N),
                   CI100_l = wmn100 - 1.96*se100,
                   CI100_u = wmn100 +  1.96*se100
                   ) 


ln_wt_avg_LUbyMO <- ddply(wt, .(MO,LU), summarise, 
                          N = sum(!is.na(lnSOC5)),   
                          wmn5 = weighted.mean(lnSOC5, group_count, na.rm=T),
                   wmn30 = weighted.mean(lnSOC30, group_count,na.rm=T),
                   wmn100 = weighted.mean(lnSOC100, group_count, na.rm=T),
                   wsd5 = sqrt(wtd.var(lnSOC5, as.numeric(group_count), na.rm=T)),
                   wsd30 = sqrt(wtd.var(lnSOC100, as.numeric(group_count), na.rm=T)),
                   wsd100 = sqrt(wtd.var(lnSOC30, as.numeric(group_count), na.rm=T)),
                   se100 = wsd100/sqrt(N),
                   CI100_l = wmn100 - 1.96*se100,
                   CI100_u = wmn100 +  1.96*se100
) 


#back transform values
wt_avg <- cbind(ln_wt_avg[,1], exp(ln_wt_avg[,2:length(ln_wt_avg)]))
names(wt_avg)[1] <- "Average For"
wt_avg_MO <- cbind(ln_wt_avg_MO[,1], exp(ln_wt_avg_MO[,2:length(ln_wt_avg_MO)]))
names(wt_avg_MO)[1] <- "Average For MO"
wt_avg_LU <- cbind(ln_wt_avg_LU[,1:2], exp(ln_wt_avg_LU[,3:length(ln_wt_avg_LU)]))
names(wt_avg_LU)[1] <- "Average For LU Class"
wt_avg_LUbyMO <- cbind(ln_wt_avg_LUbyMO[,1:2], exp(ln_wt_avg_LUbyMO[,3:length(ln_wt_avg_LUbyMO)]))


#export tables  - putting them into one file for convienence

write.table(wt_avg, paste0(wt_avg_folder,"ALL_wt_avg_wSD.csv"), sep =",", row.names = F)
write.table(wt_avg_MO, paste0(wt_avg_folder,"ALL_wt_avg_wSD.csv"), sep =",", row.names = F, append=T)
write.table(wt_avg_LU, paste0(wt_avg_folder,"ALL_wt_avg_wSD.csv"), sep =",", row.names = F, append=T)
write.table(wt_avg_LUbyMO, paste0(wt_avg_folder,"ALL_wt_avg_wSD.csv"), sep =",", row.names = F, append=T)

write.table(wt_avg_MO, paste0(wt_avg_folder,"MO_wm_all_wSD.csv"), sep =",", row.names = F)



############
#make graphs for LU's

library(ggthemes)

graphs = paste0(raca, "Graphs/")
ifelse(!dir.exists(file.path(graphs)), dir.create(graphs), FALSE)

LU <- wt_avg_LU[!is.na(wt_avg_LU[,1]), ]
names(LU)[1] <- "LU"

cplot <-   ggplot(data=LU, aes(LU, y=wmn100))+ geom_bar(stat = "identity", fill = "blue") + 
  geom_errorbar(aes(ymin=CI100_l, ymax=CI100_u)) + scale_y_continuous(limits=c(0,500), expand=c(0,0))+ scale_x_discrete(drop=FALSE) +
  ylab("Geometric Mean Mg SOC per ha to 1m") + xlab("Land Use/Cover Class") + 
  ggtitle("All MOs") 
print(cplot)

#################
#Tables - info for report
###################
#SOC stocks to 100cm
SOC100 <- ddply(bt_LUGR, .(MO, LU), summarise, Mean = mean(bt_SOC100, na.rm=T))
SOC100r <- reshape(SOC100, idvar = "MO", timevar= "LU", direction = "wide")
write.table(SOC100r, paste0(summ, "RaCA_Table_bt100_means.csv"), sep=",", row.names=F)

#SOC weighted stocks
SOC100w <- reshape(wt_avg_LUbyMO[!is.na(wt_avg_LUbyMO$wmn100),c(1:2,6)], idvar = "MO", timevar= "LU", direction = "wide")
write.table(SOC100w, paste0(summ, "RaCA_Table_bt100_WTmeans.csv"), sep=",", row.names=F)



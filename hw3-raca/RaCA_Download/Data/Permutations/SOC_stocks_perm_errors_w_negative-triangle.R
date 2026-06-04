######################################################################################################
#Permutations of measurement error for each pedon
#this file allow negative values to center error at 0
#this file uses the triangle distribution to approximate a normal distribution with a low and high limit
#####################################################################################################

rm(list = ls() ) #clear environment

#a large number of replications will cause memory issues with some machines
memory.size(max = FALSE)
memory.limit(size = NA)

############################################################
#SOC stock calculation - incorporating errors
#########################################################

#load packages
library(plyr)
library(SDMTools)
library(Hmisc)
library(triangle)
library(data.table)
#########################################################
#folder locations of source and output data - change to local path that contains input data
#sub folder paths are then updated


#raca = 'Y:/RaCA_Download/RaCA_Download/Data/' #direct path to cloudvault

#raca = "F:/RaCA-newphase/RaCA_Cloudvault/RaCA_Download/Data/"
raca = 'C:/Users/skyew/ownCloud/RaCA_Download/RaCA_Download/Data/'
#raca = 'C:/Users/skyew/Documents/data/Cloudvault/RaCA_Download/Data/'

summ = paste0(raca, "summary/")
ifelse(!dir.exists(file.path(summ)), dir.create(summ), FALSE)

perm = paste0(raca, "Permutations/")
ifelse(!dir.exists(file.path(perm)), dir.create(perm), FALSE)


perm_out = paste0(raca, "Permutations/perm_out/")
ifelse(!dir.exists(file.path(perm_out)), dir.create(perm_out), FALSE)

perm_summ = paste0(perm_out, "summaries/")
ifelse(!dir.exists(file.path(perm_summ)), dir.create(perm_summ), FALSE)



##########
#Data
###########################################
#complete sample list with lab data and modeled BD (bulk density) - July 2016
samp.s <- read.csv(paste0(raca, "RaCA_samples.csv"), na.strings = c("NA", "NULL"))

#limit analysis to central pedons
samp.s <- samp.s[samp.s$pedon_no ==1,]

#label organic (L horizons labeled as O's) and mineral horizons
samp.s$M<- as.factor(ifelse(samp.s$M== "L", "O", as.character(samp.s$M)))
samp.s$Horizon <- ifelse(samp.s$M== "O"|samp.s$M== "L", "Organic", "Mineral")

#pedon info to screen for use by depth

#remove stocks based on Use field 
SOCpedons <- read.csv(paste0(raca, "RaCA_SOC_pedons.csv"))

samp.s <- join(samp.s, SOCpedons[, c("upedon", "USE")], by = "upedon", type="left", match = "first")



#########
#Errors
#########
#import error table
me <- read.csv(paste0(perm, "Horizon_propSD.csv"))
BD_eval<- read.csv(paste0(perm, "BD_resid_eval.csv"))

#join errors to sample by Horizon type
samp <- join(samp.s, me, by="Horizon", type="left", match = "first")

samp <- join(samp, BD_eval[,c("M", "sd")], by = "M", type="left", match = "first")
names(samp)[length(samp)] <- "BDmodel.SD"

#define coarsefragment error
samp$CF.SD <- ifelse(samp$fragvolc ==0, 0, ifelse(samp$fragvolc >30, 5, 2.5))

#single bulk density record - used measured if available
samp$BD <- as.numeric(ifelse(!is.na(samp$Measure_BD), samp$Measure_BD, samp$Model_BD))

#if BD measured is used set model error to zero
samp$BDmodel.SD <- as.numeric(ifelse(!is.na(samp$Measure_BD), 0, samp$BDmodel.SD))

#remove2 samples with no lab information
samp <- samp[!is.na(samp$c_tot_ncs),]
samp <- samp[!is.na(samp$BDmodel.SD),]


#sample top, bottom, and thickness to be used in pedon stock calculaitons
samp$top <-samp$TOP
samp$bottom <- samp$BOT
samp$thick <- samp$bottom - samp$top



########################################################
#Permutations of SOC density and Stock calculations
#include SD based error
#########################################################

#samp.compare <- samp # for tracking errors caused by looping

########################################################
#create set of commands that run on depth(s) of interest
set.seed<- 337766
for (i in c(5, 30, 100)){

##########################################################################
#############################################################################  
  cdepth=i
 
  #rename data for manipulation
   #limit samples to those with complete samples for appropriate depth (use SOC stock calculation field)
   sampi <- samp[!is.na(samp$USE),]
   sampi <- sampi[as.numeric(sampi$USE)>=cdepth,]   
   sampi$upedon <- droplevels(sampi$upedon)

   
   r = 500 #number of replication to create distribution
   
  #puts a number on the amount (cm) of each layer used in calculation
  sampi$cthick <- ifelse(sampi$bottom <= cdepth, sampi$thick,
                               ifelse(sampi$top < cdepth, cdepth-(sampi$top),0))
                                                 
  #remove zeros from samples to permute
  sampi <- sampi[!is.na(sampi$cthick) & sampi$cthick!=0, ]
  upedonid <- as.character(sampi$upedon) #a label to id pedons
                                                 

#replications for distributions
  id <- c(rep(as.character(sampi$samp), each=r)) #unique sample id 
  repno <- c(rep(1:r, length(sampi$samp)))
               
  pedon <- factor(rep(sampi$upedon, each=r)) #unique pedon id 
                                                 
##errors
  TotCmeasureerror <- rep(sampi$TotC.sd, each=r)
  esptotc <- rnorm(length(id), mean = 0, sd = TotCmeasureerror)
  tri <- runif(length(TotCmeasureerror),-2*TotCmeasureerror , 2*TotCmeasureerror)
  tcf <- ifelse(tri <TotCmeasureerror & tri >TotCmeasureerror, rtriangle(1,-2*TotCmeasureerror, 2*TotCmeasureerror, 0), 0)
  tritotc <- tri+tcf
  #histogram(tritotc)

  CARBmeasureerror <- rep(sampi$CCE.sd, each=r)
  espcce <- rnorm(length(id), mean = 0, sd = CARBmeasureerror)
  ti <- runif(length(CARBmeasureerror),-2*CARBmeasureerror, 2*CARBmeasureerror)
  tcce <- ifelse(ti <CARBmeasureerror & ti >CARBmeasureerror, rtriangle(1,-2*CARBmeasureerror, 2*CARBmeasureerror, 0), 0)
  tricce <- ti+tcce
  #histogram(tricce)
  
  BDerror <- rep(sampi$bd.sd, each=r)
  espbd<- rnorm(length(id), mean = 0, sd=BDerror)
  BDmodelerror <- rep(sampi$BDmodel.SD, each=r) #add additional error for modeled BD, set to zero for measured data
  espbdm<- rnorm(length(id), mean = 0, sd=BDmodelerror) + espbd # term that combines measurement and modeling error
                                                 
  CFmodelerror <- rep(ifelse(is.na(sampi$CF.SD), 1, sampi$CF.SD +.1), each=r) #CF error calculation  
  espcf<- rnorm(length(id), mean = 0, sd=CFmodelerror)
  t <- runif(length(CFmodelerror), -2*CFmodelerror, 2*CFmodelerror)
  tcf <- ifelse(t <CFmodelerror & t >CFmodelerror, rtriangle(1,-2*CFmodelerror, 2*CFmodelerror, 0), 0)
  tricf <- t + tcf
  #histogram(tricf)
  
#input values for each sampile, replicated r times
  TotC <- c(rep(sampi$c_tot_ncs, each=r)) #Total carbon measured with combustion %
  CCE <-  c(rep(ifelse(is.na(sampi$caco3),0,sampi$caco3), each=r)) # Calcium Carbonate equivalent, pressure calcimeter, %
  BD <- c(rep(sampi$BD, each=r))	#BD measured or modeled
  CF <- c(rep(sampi$fragvolc, each = r))  #coarse fragments estimated by volume at sampling
  Cthick <- c(rep(sampi$cthick, each=r)) #thickness of layer used in each depth of stock calculation
                                                                                      
#SOC in each individual sample (includes density, coarse fragments and thickness of layer)

  #calculate SOC = total C - inorganic C, 12% of calcium carbonate is carbon (CCE*0.12)
    #add boundary of 0 for minimum carbon content
  IC <- .12*(CCE+espcce)
  TC <- TotC+esptotc
  SOCc <- TC - IC
  SOC <- SOCc
  
  ICt <- ifelse(.12*(CCE+tricce)<0,0, .12*(CCE+tricce))
  TCt <- ifelse(TotC+tritotc<0, 0, TotC+tritotc)
  SOCt <- ifelse(TotC - (.12*CCE) <0, 0, TCt - ICt)

  cf <- CF + espcf 
  CFv <- cf/100 #calculate as a ratio to whole soil
  
  cft <- ifelse(CF + tricf <0, 0, CF +tricf)
  CFt <- ifelse(cft<=0,0, cft/100) #calculate as a ratio to whole soil

  bd <- ifelse(BD+espbdm <=0, 0.01, BD+espbdm) #create a minimum possible value for bulk density
  

#samp$SOCden <- with(samp,SOC*BD*(1 - ifelse(fragvolc == 0, 0, fragvolc/100)))

SOCsdT <- SOCt* bd * (1 - CFt) * Cthick 
      
      #check calculations; used when altering calculation
      #checkcalc<-cbind.data.frame(id,repno, SOCsdT, ICt, TCt, SOCt, bd , 1-CFt, Cthick) 
                                #repno,CF, CFmodelerror, 
                                #espcf, cf, CFv, tricf, CFt)   
        #SOCsdt, TotC, TC, esptotc, CCE, espcce, IC, SOC,   BD, espbd, espbdm, bd, CF, espcf, cf, CFv, Cthick)
      ##########################################################################################################          

#remove unused data
  rm(SOCt, TCt, ICt, BDerror, BDmodelerror, CFmodelerror
     , TotCmeasureerror, CARBmeasureerror, espcce, esptotc, espcf, espbd
     , espbdm, TotC, IC, TC, SOCc, CCE, BD, CF, cf, CFv, CFt, cft
     , Cthick, SOC, id, bd)                                                                                  
                                                                                      
#########################################                                                                                      
#summarize by pedon
################
## create new dataset to run stats on for horizons  ## need to include rcasiteid then produce means and sd for each pedon
upedonrep <- paste(pedon, repno, sep= " ")

rm(pedon, repno)#remove unsed data

sumpedon<-as.data.frame.table(tapply(SOCsdT, list(upedonrep), sum))
colnames(sumpedon)<-c("upedonrep","SOC_pedon") ### post distributions for each pedon

rm(upedonrep)

##########################################################  import sumpedon
ff<-strsplit(as.character(sumpedon$upedonrep), " ")
sumpedon$pedon<-as.factor(sapply(ff, "[", 1))
sumpedon$repno<-as.integer(sapply(ff, "[", 2))
sumpedon$LU <- as.factor(substr(sumpedon$pedon, 6, 6))
sumpedon$region <- as.factor(substr(sumpedon$pedon, 2, 3))
sumpedon$LUGR <- as.factor(paste0(sumpedon$LU,substr(sumpedon$pedon, 2,5)))

upedonmean<-aggregate(SOC_pedon ~ pedon, data=sumpedon, mean)  # gC/cm2 to depth of cdepth
upedonmean$LU <- as.factor(substr(upedonmean$pedon, 6, 6))
upedonmean$region <- as.factor(substr(upedonmean$pedon, 2, 3))
upedonmean$LUGR <- as.factor(paste0(upedonmean$LU,substr(upedonmean$pedon, 2,5)))

upedonmean$SOCsd<-aggregate(SOC_pedon ~ pedon, data=sumpedon, sd)[,2]
upedonmean$SOC_q5<-aggregate(SOC_pedon ~ pedon, data=sumpedon, FUN = function(x) quantile(x, probs = 0.05))[,2]
upedonmean$SOC_q25<-aggregate(SOC_pedon ~ pedon, data=sumpedon, FUN = function(x) quantile(x, probs = 0.25))[,2]
upedonmean$SOC_q50<-aggregate(SOC_pedon ~ pedon, data=sumpedon, FUN = function(x) quantile(x, probs = 0.50))[,2]
upedonmean$SOC_q75<-aggregate(SOC_pedon ~ pedon, data=sumpedon, FUN = function(x) quantile(x, probs = 0.75))[,2]
upedonmean$SOC_q95<-aggregate(SOC_pedon ~ pedon, data=sumpedon, FUN = function(x) quantile(x, probs = 0.95))[,2]
upedonmean$calc_rep <- r
upedonmean$calc_depth <- i
title<-paste("RaCA pedonstat ",i, sep=" ")
Sys.time()
write.csv(upedonmean, file=paste0(perm_out, "NEG-mc_", i, "_", "stocks.csv"))
write.csv(sumpedon, file=paste0(perm_out, "NEG-REPS_", i, "_", "stocks.csv"))


##############
#Take permutation reps as outputs and summarize them the same way as the initial summary
###################################################################
sp <- data.table(sumpedon)

LUGR_rep1 <- sp[ , j=list(lugr.m = mean(SOC_pedon, na.rm=T)
                          , lugr.sd = sd(SOC_pedon, na.rm=T)
                          , lugr.q5 = quantile(SOC_pedon, probs = 0.05, na.rm=T)
                          ,lugr.q25 = quantile(SOC_pedon, probs = 0.25, na.rm=T)
                          ,lugr.q50 = quantile(SOC_pedon, probs = 0.50, na.rm=T)
                          ,lugr.q75 = quantile(SOC_pedon, probs = 0.75, na.rm=T)
                          ,lugr.q95 = quantile(SOC_pedon, probs = 0.95, na.rm=T))
                          , by=list(LUGR,repno)
                          ]

                        

LUGR_reps <-  as.data.frame(LUGR_rep1)

LUGR_range <- ddply(LUGR_reps, .(LUGR), summarise,
                    o.mean = mean(lugr.m, na.rm=T),
                    sd = sd(lugr.m, na.rm=T),
                    m.sd = mean(lugr.sd, na.rm=T),
                    m.q5 = mean(lugr.q5, na.rm=T),
                    m.q25 = mean(lugr.q25, na.rm=T),
                    m.q50 = mean(lugr.q50, na.rm=T),
                    m.q75 = mean(lugr.q75 , na.rm=T),
                    m.q95 = mean(lugr.q95, na.rm=T),
                    q5.wm =  quantile(lugr.m, probs = 0.05, na.rm=T),
                    q25.wm =  quantile(lugr.m, probs = 0.25, na.rm=T),
                    q5.wm =  quantile(lugr.m, probs = 0.5, na.rm=T),
                    q75.wm =  quantile(lugr.m, probs = 0.75, na.rm=T),
                    q95.wm =  quantile(lugr.m, probs = 0.95, na.rm=T)
)
#
#
write.csv(LUGR_reps, file=paste0(perm_out, "NEG-mc_", r, "_reps_", i, "_LUGR_", "meanstocks_REPS.csv"))
write.csv(LUGR_range, file=paste0(perm_out, "NEG-mc_", r, "_reps_", i, "_LUGR_", "meanstocks_RANGE.csv"))


#Add weights for LU and MO by LUGR (for selection probabilities)
prob <- read.csv(paste0(perm, "LUGR_prob.csv"))

P <- join(LUGR_reps, prob[,c(2:3,7, 12:15)], by="LUGR")
#remove X (na for weights) - as it doesn't have spatial meaning
P <- P[!is.na(P$LU),]


oa_reps <-  ddply(P, .(repno), summarise,
                  s.m = mean(lugr.m, na.rm=T),
                  w.m = weighted.mean(lugr.m, w = prob.oa),
                  sd= sd(lugr.m, na.rm=T),
                  #wt.sd = sqrt(wtd.var(lugr.m, w = prob.oa, na.rm=T)),
                  q5 = quantile(lugr.m, probs = 0.05, na.rm=T),
                  q25 = quantile(lugr.m, probs = 0.25, na.rm=T),
                  q50 = quantile(lugr.m, probs = 0.50, na.rm=T),
                  q75 = quantile(lugr.m, probs = 0.75, na.rm=T),
                  q95 = quantile(lugr.m, probs = 0.95, na.rm=T)
)
#
oa_summ <- ddply(oa_reps, .(), summarise,
                 std.m = mean(s.m),
                 wt.m = mean(w.m),
                 sd.wm = sd(w.m),
                 m.sd = mean(sd),
                 #wt.sd = mean(wt.sd),
                 m.q5 = mean(q5),
                 m.q25 = mean(q25),
                 m.q50 = mean(q50),
                 m.q75 = mean(q75),
                 m.q95 = mean(q95),
                 q05.wm =  quantile(w.m, probs = 0.05, na.rm=T),
                 q25.wm =  quantile(w.m, probs = 0.25, na.rm=T),
                 q50.wm =  quantile(w.m, probs = 0.5, na.rm=T),
                 q75.wm =  quantile(w.m, probs = 0.75, na.rm=T),
                 q95.wm =  quantile(w.m, probs = 0.95, na.rm=T)
)

LU_reps <-  ddply(P, .(repno, LU), summarise,
                  s.m = mean(lugr.m, na.rm=T),
                  w.m = weighted.mean(lugr.m, w = prob.LU),
                  sd= sd(lugr.m, na.rm=T),
                  wt.sd = sqrt(wt.var(lugr.m, w = prob.LU)),
                  q5 = quantile(lugr.m, probs= 0.05, na.rm=T),
                  q25 = quantile(lugr.m, probs = 0.25, na.rm=T),
                  q50 = quantile(lugr.m, probs = 0.50, na.rm=T),
                  q75 = quantile(lugr.m, probs = 0.75, na.rm=T),
                  q95 = quantile(lugr.m, probs = 0.95, na.rm=T),
                  q05.wm =  quantile(w.m, probs = 0.05, na.rm=T),
                  q25.wm =  quantile(w.m, probs = 0.25, na.rm=T),
                  q50.wm =  quantile(w.m, probs = 0.5, na.rm=T),
                  q75.wm =  quantile(w.m, probs = 0.75, na.rm=T),
                  q95.wm =  quantile(w.m, probs = 0.95, na.rm=T)
)

#
LU_summ <-  ddply(LU_reps, .(LU), summarise,
                  m = mean(w.m, na.rm=T),
                  m.sm = mean(s.m),
                  m.wsd = mean(wt.sd),
                  m.sd = mean(sd),
                  sd.wm = sd(w.m, na.rm=T),
                  m.q5 = mean(q5),
                  m.q25 = mean(q25),
                  m.q50 = mean(q50),
                  m.q75 = mean(q75),
                  m.q95 = mean(q95),
                  q05.wm =  quantile(w.m, probs = 0.05, na.rm=T),
                  q25.wm =  quantile(w.m, probs = 0.25, na.rm=T),
                  q50.wm =  quantile(w.m, probs = 0.5, na.rm=T),
                  q75.wm =  quantile(w.m, probs = 0.75, na.rm=T),
                  q95.wm =  quantile(w.m, probs = 0.95, na.rm=T)
)

MO_reps <-  ddply(P, .(repno, MO), summarise,
                  s.m = mean(lugr.m, na.rm=T),
                  w.m = weighted.mean(lugr.m, w = prob.MO),
                  sd = sd(lugr.m, na.rm=T),
                  wt.sd = wt.sd(lugr.m, w = prob.MO),
                  q05 = quantile(lugr.m, probs = 0.05, na.rm=T),
                  q25 = quantile(lugr.m, probs = 0.25, na.rm=T),
                  q50 = quantile(lugr.m, probs = 0.50, na.rm=T),
                  q75 = quantile(lugr.m, probs = 0.75, na.rm=T),
                  q95 = quantile(lugr.m, probs = 0.95, na.rm=T)
)

MO_summ <-  ddply(MO_reps, .(MO), summarise,
                  m = mean(w.m, na.rm=T),
                  s.m = mean(s.m),
                  m.sd = mean(wt.sd),
                  sd.m = sd(w.m, na.rm=T),
                  m.q5 = mean(q05),
                  m.q25 = mean(q25),
                  m.q50 = mean(q50),
                  m.q75 = mean(q75),
                  m.q95 = mean(q95),
                  q05.wm =  quantile(w.m, probs = 0.05, na.rm=T),
                  q25.wm =  quantile(w.m, probs = 0.25, na.rm=T),
                  q50.wm =  quantile(w.m, probs = 0.5, na.rm=T),
                  q75.wm =  quantile(w.m, probs = 0.75, na.rm=T),
                  q95.wm =  quantile(w.m, probs = 0.95, na.rm=T)
                  
                  
)


LUbyMO_reps <-  ddply(P, .(repno, LU, MO), summarise,
                      s.m = mean(lugr.m, na.rm=T),
                      w.m = weighted.mean(lugr.m, w = prob.LUbyMO),
                      sd.sd = mean(lugr.sd, na.rm=T),
                      sd.mean = sd(lugr.m, na.rm=T),
                      wt.sd = wt.sd(lugr.m, w = prob.LUbyMO),
                      q5 = quantile(lugr.m, probs = 0.05, na.rm=T),
                      q25 = quantile(lugr.m, probs = 0.25, na.rm=T),
                      q50 = quantile(lugr.m, probs = 0.50, na.rm=T),
                      q75 = quantile(lugr.m, probs = 0.75, na.rm=T),
                      q95 = quantile(lugr.m, probs = 0.95, na.rm=T)
)

LUbyMO_summ <-  ddply(LUbyMO_reps, .(MO), summarise,
                      m = mean(w.m, na.rm=T),
                      s.m = mean(s.m),
                      m.sd = mean(sd.sd),
                      m.wsd = mean(wt.sd),
                      sd.m = sd(w.m, na.rm=T),
                      m.q5 = mean(q5),
                      m.q25 = mean(q25),
                      m.q50 = mean(q50),
                      m.q75 = mean(q75),
                      m.q95 = mean(q95),
                      q05.wm =  quantile(w.m, probs = 0.05, na.rm=T),
                      q25.wm =  quantile(w.m, probs = 0.25, na.rm=T),
                      q50.wm =  quantile(w.m, probs = 0.5, na.rm=T),
                      q75.wm =  quantile(w.m, probs = 0.75, na.rm=T),
                      q95.wm =  quantile(w.m, probs = 0.95, na.rm=T)
                      
)



write.csv(MO_reps, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "MO_reps.csv"))
write.csv(LU_reps, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "LU_REPS.csv"))
write.csv(oa_reps, file=paste0(perm_summ, "NEG-mc_", r, "reps_",i, "cm_LUGR_", "OA_reps.csv"))
write.csv(MO_summ, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "MO_summ.csv"))
write.csv(LU_summ, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "LU_summ.csv"))
write.csv(oa_summ, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "OA_summ.csv"))
write.csv(LUbyMO_reps, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "LUbyMO_reps.csv"))
write.csv(LUbyMO_summ, file=paste0(perm_summ, "NEG-mc_",r, "reps_", i, "cm_LUGR_", "LUbyMO_summ.csv"))


#rm(list=ls(all=TRUE))
rm(list=setdiff(ls(), c("raca","perm", "perm_out", "perm_summ","samp", "prob","SU", "summ", "LUGR_range", "LUGR_reps",'MO_reps', 'LU_pedreps','LU_reps', "MO_reps", "OA_reps",
                        "OA_summ", "MO_summ", "LU_summ")))
}

###################################################################                                                                 


#example                                                                                        
#Count.t <- reshape(count, idvar="MO", timevar = "LU", direction = "wide" )
#write.table(Count.t, paste0(raca, "RaCA_Table_count.csv"), sep=",", row.names=F)
# reorganize tables for output


#save sumpedons to run summaries

                                                                                        
                                                                                        
                                                                                        
                                                                                      
rm(list = ls() ) #clear environment
############################################################
#SOC stock calculation - using RaCA sample values (from central pedon)
#######################################################3

###--- Install and Load Libraries as needed
#instal packages
#install.packages(c("ggplot2","soilDB","plyr", "reshape2", "ggthems"))

#load packages
library(soilDB)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggthemes)



# #currently download data from USDA Cloudvault then save in local folder
# "https://www.cloudvault.usda.gov/index.php/s/nZDwWAdlmKUW1g4"


##############################################
#set folders and subfolders for data and  output

#raca = "F:/RaCA-newphase/RaCA_Cloudvault/RaCA_Download/Data/"
#raca = 'C:/Users/skyew/Documents/data/'
raca = 'C:/Users/skyew/ownCloud/RaCA_Download/RaCA_Download/Data/'

summ= paste0(raca, "Summary/")
ifelse(!dir.exists(file.path(summ)), dir.create(summ), FALSE) 

#complete sample list with lab data - July 2016
samp <- read.csv(paste0(raca, "RaCA_samples.csv"), na.strings = c("NA", "NULL"))
str(samp)

#limit calculations to central pedons
samp <- samp[samp$pedon_no == 1,]

#count sites
agg <- aggregate(data=samp, rcasiteid ~ MO + LU , function(x) length(unique(x)))
Count.r <- reshape(agg, idvar="MO", timevar = "LU", direction = "wide" )
write.csv(Count.r, paste0(raca, "RaCA_site_count.csv"))

#simplify variables and calculate SOC density for each sample
samp$CALC_SOC <- ifelse(!is.na(samp$caco3),
                        ifelse(samp$caco3>0,samp$c_tot_ncs-(samp$caco3*0.12),samp$c_tot_ncs), samp$c_tot_ncs) 
samp$SOC <- ifelse(samp$CALC_SOC<0, 0, samp$CALC_SOC)
samp$BD <- ifelse(!is.na(samp$Measure_BD), samp$Measure_BD, samp$Model_BD)
samp$SOCd <- with(samp,SOC*BD)
samp$SOCden <- with(samp,SOC*BD*(1 - ifelse(fragvolc == 0, 0, fragvolc/100)))

samp$top <-samp$TOP
samp$bottom <- samp$BOT

samp$thick <- samp$bottom - samp$top

#assign how layers are counted in stock calculations
#puts a number on the amount (cm) of each layer used in calculation
samp$SOC_5 <- ifelse(samp$bottom<=5, samp$thick,
                      ifelse(samp$top<5, 5-samp$top,0))

samp$SOC_30 <- ifelse(samp$bottom<=30, samp$thick,
                    ifelse(samp$top<30, 30-samp$top,0))

samp$SOC_100 <- ifelse(samp$bottom<=100, samp$thick,
                    ifelse(samp$top<100, 100-samp$top,0))

#calculate layer stocks - assigned depth * SOCden 
#native units are Mg/ha
#chack math against pedon example excel and Kabindra script below
samp$SOCstock5 <- with(samp, SOCden*SOC_5)
samp$SOCstock30 <- with(samp, SOCden*SOC_30)
samp$SOCstock100 <- with(samp, SOCden*SOC_100)

#sum SOC by pedon
SOC_5 <- aggregate(SOCstock5~MO + MOGr + LU + MOGrLU + rcasiteid + upedonid + upedon, data=samp, FUN=sum)
SOC_30 <- aggregate(SOCstock30~upedonid, data=samp, FUN=sum)
SOC_100 <-aggregate(SOCstock100~upedonid, data=samp, FUN=sum) 

#pedon id info
ped <- aggregate(sample.id ~MO + MOGr + LU + MOGrLU + rcasiteid + upedonid + upedon, data=samp, FUN=length)
names(ped)[8] <- "Sample_count"

#pedon information
pedon_thick <- aggregate(thick  ~ upedonid, data=samp, FUN=sum)
names(pedon_thick)[2] <- "total_thickness"
ped_bott <- aggregate(bottom~upedonid, data=samp, FUN=max)
lab_no <- aggregate(natural_key~upedonid, data=samp, FUN=length)
names(lab_no)[2] <- "Lab_count"

sample_notR <- aggregate(sample.id~upedonid, data=samp[grep("R", samp$model_desg, invert=T, ignore.case=T),] , FUN=length)
names(sample_notR)[2] <- "Non-R_SampleCount"

SOC_no <- aggregate(SOC~upedonid, data=samp[!is.na(samp$SOC), ], FUN=length)
names(SOC_no)[2] <- "SOC_count"

SOC_thick <- aggregate(thick~upedonid, data=samp[!is.na(samp$SOC), ], FUN=sum)
names(SOC_thick)[2] <- "SOC_thickness"

top <- aggregate(top~upedonid, data=samp, FUN=min)
names(top)[2] <- "Pedon_top"

#depth of R or other restrictive layer (assumed to have negligble SOC)

hard <- c("Bkm", "Bkmm", "Bqm", "R", "Cr", "m")
depthR <- aggregate(top~upedonid, data=samp[grepl(paste(hard, collapse='|'),  samp$model_desg, ignore.case=T),], FUN=min)
#only uses R, added Cr, massive and other above
#depthR <- aggregate(top~upedonid, data=samp[grepl("R", samp$model_desg, ignore.case=T),], FUN=min)
names(depthR)[2] <- "Depth_to_R"

#combine pedon information
### this method didn't work..........m <- merge_all(SOC_5,SOC_30, SOC_100, by=upedonid)
SOC5 <- join(ped, SOC_5, by="upedonid", type="full", match="first")
SOC30 <- join(SOC5, SOC_30, by="upedonid", type="full", match="first")
SOC100 <- join(SOC30,SOC_100, by="upedonid", type="full", match="first")
SOCp <- join(SOC100, pedon_thick, by="upedonid", type="full", match="first")
SOCn <- join(SOCp, SOC_no, by="upedonid", type="full", match="first")
SOCnn <- join(SOCn, lab_no, by="upedonid", type="full", match="first")
SOCpn <- join(SOCnn, depthR, by="upedonid", type="full", match="first")
SOCpp <- join(SOCpn, sample_notR, by="upedonid", type="full", match="first")
SOCpedons <- join(SOCpp, SOC_thick, by="upedonid", type="full", match="first")
                                                      
SOCpedons$USE <-ifelse(SOCpedons$SOCstock5 == 0, 
                                                ifelse(SOCpedons$SOC_count<SOCpedons$Lab_count,"MISS_anal", "MISS_labsamp"), 
                       ifelse(is.na(SOCpedons$total_thickness), 'NA',
                              ifelse(SOCpedons$SOC_thickness >= SOCpedons$total_thickness, '100',
                                     ifelse(!is.na(SOCpedons$Depth_to_R),  
                                            ifelse(SOCpedons$Depth_to_R <= SOCpedons$SOC_thickness , '100', 
                                                   ifelse(SOCpedons$SOC_thickness >= 100, '100',
                                                          ifelse(SOCpedons$SOC_thickness >= 30, '30',
                                                                 ifelse(SOCpedons$SOC_thickness >= 5, '5', 'NA')))),
                                            ifelse(SOCpedons$SOC_thickness >= 100, '100',
                                                   ifelse(SOCpedons$SOC_thickness >= 30, '30',
                                                          ifelse(SOCpedons$SOC_thickness >= 5, '5', 'NA'))))
                              )))
                        

# check NAs and missing samples
table(SOCpedons$MO,SOCpedons$USE, useNA= "ifany")

#keep ones with values
SOCpedons <- SOCpedons[!is.na(SOCpedons$SOCstock5),]
SOCpedons <- SOCpedons[SOCpedons$SOCstock5 != 0,]

#remove stocks based on Use field 
SOCpedons$SOCstock5 <- as.numeric(ifelse(SOCpedons$USE %in% c(5, 30, 100), SOCpedons$SOCstock5, ""))
SOCpedons$SOCstock30 <- as.numeric(ifelse(SOCpedons$USE %in% c(30, 100), SOCpedons$SOCstock30, ""))
SOCpedons$SOCstock100 <- as.numeric(ifelse(SOCpedons$USE %in% c(100), SOCpedons$SOCstock100, ""))


summary(SOCpedons$SOCstock30)
summary(SOCpedons$SOCstock100)


#export
#############
write.table(SOCpedons, paste0(raca, "RaCA_SOC_pedons.csv"), sep=",", row.names=F)



############
#make graphs
#####################
# To save graphs in a new folder


graphs = paste0(summ, "Graphs/")
ifelse(!dir.exists(file.path(graphs)), dir.create(graphs), FALSE)


#add data if not currently available
#SOCpedons <- read.csv(paste0(raca, "RaCA_SOC_pedons.csv"))

#  change labels LU with Land use - land cover classes
levels(SOCpedons$LU)
#reorder land use
SOCpedons$LU <- factor(SOCpedons$LU, levels(SOCpedons$LU)[c(5,2,3,4,1,6)])
levels(SOCpedons$LU)
levels(SOCpedons$LU)<- c('Wetland', 'Forest land', 'Pastureland',  "Rangeland", 'Cropland', "CRP")


#create graph using ggplot
C1 <- ggplot(SOCpedons, aes(SOCstock100, color = LU)) + stat_ecdf( size=2) 
Ct <- C1 +scale_colour_few() + labs(color =  "Land Use - Land Cover Class", y = "ECDF - F(x)") + theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=10)) 

#view graph
Ct

#save graph
ggsave('Fig6c.Cumm_SOC100.png', plot = Ct, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


#desntiy plots


D1 <- ggplot(SOCpedons, aes(SOCstock100)) + geom_density(aes(fill=LU), alpha= .75) 
Dt <- D1+scale_fill_few()  + labs(fill ="Land Use - Land Cover Class", x="SOC stock to 100cm") +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),legend.text=element_text(size=12)) 

Dt

ggsave('Fig.6a.Density_SOC100.png', plot = Dt, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)

Dl <- D1 + labs(fill ="Land Use - Land Cover Class", x="Log SOC stock to 100cm") + scale_fill_few() + scale_x_log10() +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),legend.text=element_text(size=12)) 


Dl

ggsave('Fig.6b.Density_LOGSOC100.png', plot = Dl, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


#boxplots
#use only those with values to 100cm
S <- SOCpedons[!is.na(SOCpedons$SOCstock100),]

#check levels of LULC; change and relabel if needed
levels(S$LU)

#100cm
s <- na.omit(S[ ,c("LU","MO","SOCstock100")]) #limit data used to those with relevant values

B <-  ggplot(s, aes(y =SOCstock100, x = LU)) + geom_boxplot()

Bl <- B + labs(x ="Land Use - Land Cover Class", y="SOC stock to 100cm") + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

Bl

ggsave('Fig.8c2.DBoxplot_SOC100.png', plot = Bl, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)

#30cm

s30 <- na.omit(S[,c("LU","MO","SOCstock30")]) 

B30 <-  ggplot(s30, aes(y =SOCstock30, x = LU)) + geom_boxplot()
+ scale_y_log10() 

B30l <- B30 + labs(x ="Land Use - Land Cover Class", y="SOC stock to 30cm") + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

B30l

ggsave('Fig.8b2.DBoxplot_SOC30.png', plot = B30l, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


#5cm

s5 <- na.omit(S[,c("LU","MO","SOCstock30")]) 

B5 <-  ggplot(s5, aes(y =SOCstock30, x = LU)) + geom_boxplot()


B5l <- B5 + labs(x ="Land Use - Land Cover Class", y="SOC stock to 5cm") + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

B5l

ggsave('Fig.8a2.DBoxplot_SOC5.png', plot = B5l, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


#Change scales to log for each depth increment
B5L <- B5 + labs(x ="Land Use - Land Cover Class", y="Log SOC stock to 5cm") + scale_y_log10() + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

B5L

ggsave('Fig.8a1.DBoxplot_LogSOC5.png', plot = B5L, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)

B30l <- B30 + labs(x ="Land Use - Land Cover Class", y="Log SOC stock to 30cm") + scale_y_log10() +  
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

B30l

ggsave('Fig.8b1.DBoxplot_LogSOC30.png', plot = B30l, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


Bl <- B + labs(x ="Land Use - Land Cover Class", y="Log SOC stock to 100cm") + scale_y_log10() +  
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 

Bl

ggsave('Fig.8c1.DBoxplot_logSOC100.png', plot = Bl, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)

################################################
#explore samples
#sample graphs

#reorder LU
levels(samp$LU)
samp$LU <- factor(samp$LU, levels(samp$LU)[c(5,2,3,4,1,6)])
levels(samp$LU)
levels(samp$LU)<- c('Wetland', 'Forest land', 'Pastureland',  "Rangeland", 'Cropland', "CRP")


s <- ggplot(data=samp, aes(SOC)) + geom_density(aes(fill=LU), alpha= .75) + scale_x_log10()


s1 <- s +scale_fill_few() + labs(fill ="Land Use - Land Cover Class", x="SOC concentration - Log Scale") + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=12)) 


s1

ggsave('Fig.3.Density_Log_SOC.png', plot = s1, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)



m <- ggplot(samp[samp$M %in% c("A", "E", "B", "C", "O","L"),], aes(M, SOC)) + geom_boxplot(aes(fill=M), show.legend = F) + scale_y_log10()

m1 <- m + labs(x ="Horizon Master", y="SOC concentration - Log Scale") +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),legend.text=element_text(size=12)) 


m1

# ggsave('BoxbyM_LogSOCconc.png', plot = m1, device = "png", path = graphs, scale = 1, width = 6,   
#        height = 3, units = "in", dpi = 600, limitsize = TRUE)


t <- ggplot(samp, aes(texture, SOC)) + geom_boxplot(aes(fill=texture), show.legend = F) + scale_y_log10()
t

t1 <- t + labs(x ="texture class", y="SOC concentration - Log Scale") +
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=8),legend.text=element_text(size=12)) 

t1

#surface SOC
r <- ggplot(samp[samp$TOP==0,], aes(LU, SOC)) + geom_boxplot() + scale_y_log10() 

cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# scale_fill_manual(values=cbPalette)

r1 <- r + labs(x ="RaCA Region", y="SOC concentration - Log Scale") +scale_fill_manual(values=cbPalette)+
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=8),legend.text=element_text(size=12)) 

r1


# ggsave('BoxbyM_LogSOCconc.png', plot = m1, device = "png", path = graphs, scale = 1, width = 6,   
#        height = 3, units = "in", dpi = 600, limitsize = TRUE)


#Master Horizon by LU
levels(samp$M)
samp$M <- factor(samp$M, levels(samp$M)[c(6,5,1,4,2,3,7,8,9)])
levels(samp$M)


#create a function for a contigency/mosaic plot

ggMMplot <- function(var1, var2){
  require(ggplot2)
  levVar1 <- length(levels(var1))
  levVar2 <- length(levels(var2))
  
  jointTable <- prop.table(table(var1, var2))
  plotData <- as.data.frame(jointTable)
  plotData$marginVar1 <- prop.table(table(var1))
  plotData$var2Height <- plotData$Freq / plotData$marginVar1
  plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) +
    plotData$marginVar1 / 2
  
  ggplot(plotData, aes(var1Center, var2Height)) +
    geom_bar(stat = "identity", aes(width = marginVar1, fill = var2), col = "Black") +
    geom_text(aes(label = as.character(var1), x = var1Center, y = 1.05)) 
}

h <- samp[samp$M %in% c("A", "E", "B", "C", "O", "R"),]

h$M <- droplevels(h$M)

levels(h$M)
#h$M <- factor(h$M, levels(h$M)[c(5,1, 4, 2,3,6)])
#levels(h$M)


ml <- ggMMplot(h$M,h$LU)
# 
# cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# scale_fill_manual(values=cbPalette)

ml1 <- ml + labs(x ="Master Horizon Proportion", y="LULC Proportion", fill = "Land Use - Land Cover") +  scale_fill_brewer(palette="Set1")  +
  theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.text.x=element_text(size=10),legend.text=element_text(size=10)) 

ml1

ggsave('Fig.4.Mosaic_LUbyM.png', plot = ml1, device = "png", path = graphs, scale = 1, width = 6,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


####################################################
# Depth Plots for SOC concentration and Density
#using aqp
####################################################
library(aqp)
library(lattice)

x <- samp

# promote to SPC
depths(x) <- upedon  ~ TOP + BOT
# move some site-level data to site slot
site(x) <- ~ rcasiteid + MO + LU

# depth-wise quantiles, by LU
a <- slab(x, LU ~ SOC + BD + SOCden)

levels(a$LU)
# fix LU factor levels - if not done previously, check first
levels(a$LU)
#a$LU <- factor(a$LU, levels=c('C', 'F', 'P', 'R', 'W', 'X'), labels=c('Cropland', 'Forestland', 'Pastureland', 'Rangeland', 'Wetland', 'Other'))

# custom colors
tps <- list(superpose.line=list(col=c('DarkRed', 'DarkGreen', 'grey', "RoyalBlue", 'black', 'Yellow'), lwd=2))


p.SOC <- xyplot(top ~ p.q50, groups=LU, data=a, ylab='Depth (cm)',
                xlab='median bounded by 25th and 75th percentiles',
                lower=a$p.q25, upper=a$p.q75, ylim=c(40,-1),
                panel=panel.depth_function, alpha=0.25, sync.colors=TRUE,
                prepanel=prepanel.depth_function,
                strip=strip.custom(bg=grey(0.85)),
                scales=list(x=list(alternating=1)),
                par.settings=tps, subset=variable == 'SOC',
                auto.key=list(columns=2, lines=TRUE, points=FALSE)
)

p.BD <- xyplot(top ~ p.q50, groups=LU, data=a, ylab='Depth (cm)',
               xlab='median bounded by 25th and 75th percentiles',
               lower=a$p.q25, upper=a$p.q75, ylim=c(40,-1),
               panel=panel.depth_function, alpha=0.25, sync.colors=TRUE,
               prepanel=prepanel.depth_function,
               strip=strip.custom(bg=grey(0.85)),
               scales=list(x=list(alternating=1)),
               par.settings=tps, subset=variable == 'BD',
               auto.key=list(columns=2, lines=TRUE, points=FALSE)
)

p.SOCden <- xyplot(top ~ p.q50, groups=LU, data=a, ylab='Depth (cm)',
                   xlab='median bounded by 25th and 75th percentiles',
                   lower=a$p.q25, upper=a$p.q75, ylim=c(40,-1),
                   panel=panel.depth_function, alpha=0.25, sync.colors=TRUE,
                   prepanel=prepanel.depth_function,
                   strip=strip.custom(bg=grey(0.85)),
                   scales=list(x=list(alternating=1)),
                   par.settings=tps, subset=variable == 'SOCden',
                   auto.key=list(columns=2, lines=TRUE, points=FALSE)
)

png(file = paste0(graphs,'Fig.5a.SOC.png'), width = 4, height = 4, units = 'in', res = 600)
print(p.SOC) # Make plot
dev.off()

png(file = paste0(graphs,'Fig.5b.png'), width = 4, height = 4, units = 'in', res = 600)
print(p.BD) # Make plot
dev.off()

png(file = paste0(graphs,'Fig.5c.SOCden.png'), width = 4, height = 4, units = 'in', res = 600)
print(p.SOCden) # Make plot
dev.off()


#######################
#Create stacked bar graphs
#to use SOCpedons from calculation

O <- SOCpedons

# check LU order
levels(O$LU)
O$LU <- factor(O$LU, levels(O$LU)[c(5,2,3,4,1,6)])
levels(O$LU)
levels(O$LU)<- c('Wetland', 'Forest land', 'Pastureland',  "Rangeland", 'Cropland', "CRP")


#create factors for stacked bar plots
O$zero_to_five <- O$SOCstock5
O$five_to_thirty <- O$SOCstock30-O$SOCstock5
O$thirty_to_hundred <- O$SOCstock100- O$SOCstock30 - O$SOCstock5

#change orientation
t <- melt(O, id = names(O[,c("upedonid", "MO", "LU", "LUGR")]), 
          measure.vars= names(O[,c("zero_to_five", "five_to_thirty", "thirty_to_hundred")]))
names(t)
names(t)[5] <- "Depth"
names(t)[6] <- "SOCstock"

#need to calculate metrics by depth first, I think
M <- ddply(t, .(LU, Depth), summarise,
           N = length(SOCstock),
           SOC = mean(SOCstock, na.rm=T),
           sd = sd(SOCstock, na.rm=T),
           se = sd/sqrt(N),
           CL_l = SOC-(se*1.96),
           CL_u = SOC+ (se*1.96)
)


#prep levels for plotting
levels(M$Depth) <- c("0 to 5 cm","5 to 30 cm","30 to 100cm")

St <- ggplot(M, aes( LU, SOC)) + geom_bar(stat = "identity", aes(fill=Depth)) +   scale_y_reverse() 

St1 <- St + labs(x ="Land Use - Land Cover Class", y="SOC stock (Mg/ha)", fill = "Depth Increment") + 
  theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=10),legend.text=element_text(size=10)) 


St1

ggsave('Fig.9.STACK_SOC.png', plot = St1, device = "png", path = graphs, scale = 1, width = 8,   
       height = 3, units = "in", dpi = 600, limitsize = TRUE)


Ste <- St1 +  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))

Ste

# ggsave('STACK_SOC_alt.png', plot = Ste, device = "png", path = graphs, scale = 1, width = 8,   
#        height = 3, units = "in", dpi = 600, limitsize = TRUE)

#######LIBRARIES##########
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(plyr)
library(reshape)
library(foreign)
library(scales)
library(vegan)
library(pastecs)
library(data.table)
library(lme4)
library(lqmm)
library(bbmle)
library(quantreg)
library(hexbin)
library(viridis)
library(gam)
library(mgcv)
library(nortest)
library(png)
library(grid)
library(car)
library(boot)

rootdir <- "F:/Chapter3_ecoevo/" #UPDATE!
datadir <- file.path(rootdir, 'data')
resdir <- file.path(rootdir, 'results')
figdir <- file.path(rootdir, 'doc/Manuscript/Figures')

##########################################################################################################################################################
################################################################################## 1. FORMAT DATA #####################################################
############## Format CPUE ########
CPUEOR <- read.csv(file.path(datadir,"Fieldworkdata/CPUE_OR.csv"))

#Format CPUE for mapping
CPUEOR_kickmean <- ddply(CPUEOR[CPUEOR$Method == "Kick",], c("Site", "Method"), summarise, kickmean = mean(CPUE/Effort), kickstd=sd(CPUE/Effort))
CPUEOR_dipmean <- ddply(CPUEOR[CPUEOR$Method == "Dipnetting",], c("Site", "Method"), summarise, dipmean = mean(CPUE/Effort),  dipstd=sd(CPUE/Effort))
CPUEOR_mean <- ddply(CPUEOR, c("Site", "Method"), summarise, mean = mean(CPUE/Effort), std=sd(CPUE/Effort))
CPUEOR_melt <- melt(CPUEOR_mean, id.var=c('Site','Method'))
CPUEOR_cast <- cast(CPUEOR_melt, Site~Method+variable)
CPUEOR_cast <- CPUEOR_cast[,-which(colnames(CPUEOR_cast) %in% c("Dipnetting_std","Snork/Dipnetting_mean","Snork/Dipnetting_std","Snorkeling_std", "Snorkeling/Dipnetting_mean","Snorkeling/Dipnetting_std"))]
#After joining table to sites, make sure to create new field and replace NA by NULL if needed
# write.dbf(CPUEOR_mean, "F:/Chapter2_HexSim_Crayfish/data/Field_work_Data/CPUE_OR_stat.dbf")
# write.dbf(CPUEOR_kickmean, "F:/Chapter2_HexSim_Crayfish/data/Field_work_Data/CPUEkick_OR_stat.dbf")
# write.dbf(CPUEOR_dipmean, "F:/Chapter2_HexSim_Crayfish/data/Field_work_Data/CPUEdip_OR_stat.dbf")

################################################# FORMAT SITE-LEVEL HABITAT DATA #############################
#Data on distance between site and closest downstream confluence
site_streamdist <- read.csv(file.path(datadir,"Fieldworkdata/Sampled_sites_gps_streamdist_edit.csv"))
#Habitat data (Velocity, depth, temperature,etc.)
habdat <- read.csv(file.path(datadir,"Fieldworkdata/Site_habdat3.csv"))
#Sampling info for site
siteinfo <- read.csv(file.path(datadir,"Fieldworkdata/Site_info.csv"))
#Macroinvertebrate data
macroinv <- read.csv(file.path(datadir,"Bug_Labwork_spreadsheet.csv"))
macroinv$DW_net <- macroinv$DW - macroinv$Dish_weight  
macroinv$AFDW <- macroinv$DW_net - (macroinv$AW/1000)

site_streamdist$Site <- substr(site_streamdist$Name, 4, 10)
#Summarize habitat data
habdat_stat <- ddply(habdat, .(Site),summarize, 
                     tempmean = mean(Temperature), 
                     tempsd = sd(Temperature),
                     depthmean = mean(Depth),
                     depthsd = sd(Depth),
                     velmean = mean(Velocity),
                     velsd = sd(Velocity),
                     greenmean = mean(Green, na.rm = T), 
                     greensd = sd(Green, na.rm = T),
                     cyanomean = mean(Cyano, na.rm = T),
                     cyanosd = sd(Cyano, na.rm = T),
                     diatommean = mean(Diatom, na.rm = T),
                     diatomsd = sd(Diatom, na.rm = T))
#Summarize macroinvertebrate biomass data
#ggplot(macroinv, aes(x=Site, y=AFDW, color=factor(Site))) + geom_point(size=4, alpha=1/2) + geom_text(aes(label=Site))
#Take out outlier from Holliday Clyde with 111g of trichoptera for analysis
macroinv <- macroinv[macroinv$AFDW<2,]
macroinv_stat <- ddply(macroinv[!is.na(macroinv$DW_net),], .(Site),summarize, 
                     DWmean = mean(DW_net, na.rm =T), 
                     DWsd = sd(DW_net, na.rm = T),
                     AFDWmean = mean(AFDW, na.rm=T),
                     AFDWsd = sd(AFDW, na.rm=T),
                     macn = length(DW_net))


habdatdist <- merge(habdat_stat, site_streamdist, by = "Site", all.x=T,all.y=T)
habdatdist <- merge(habdatdist, macroinv_stat, by = "Site", all.x=T,all.y=T)
habdatdistinfo <- merge(siteinfo, habdatdist, by = "Site", all.x=T, all.y=T)
habdatdistinfo <- merge(habdatdistinfo, CPUEOR_cast, by='Site', all.x=T)

############## Format distance from confluence into distance from introduction point ###################
#Needed to edit names of tributaries in edges_18 to get correct distances (Bridge Creek, Cottonwood Creek, Canyon Creek, change direction of Beech Creek)
#Introduction point
introdist <- 379179

#Distance between two confluences (JD is JD mouth)
JDBeech <- 377903
JDBridge <- 213770
JDCanyon <- 391271
JDCottonwood <- 326853
BeechEastForkBeech <- 17503
EastForkBeechLake <- 10725
SFMurderers <- 26159
JDNF <- 292428
JDRock <- 322403
JDSF <- 335231

habdatdistinfo[habdatdistinfo$River_Tributary == 'Beech Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Beech Creek', 'interval'] + introdist-JDBeech
habdatdistinfo[habdatdistinfo$River_Tributary == 'Bridge Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Bridge Creek', 'interval'] + abs(JDBridge-introdist)
habdatdistinfo[habdatdistinfo$River_Tributary == 'Canyon Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Canyon Creek', 'interval'] + abs(JDCanyon-introdist)
#NHD confluence of Cottonwood Creek and mainstem is innacurate. So point ends up sitting on John Day River
habdatdistinfo[habdatdistinfo$River_Tributary == 'Cottonwood Creek', 'Spread_dist'] <- 250 + abs(introdist - JDCottonwood)
#Do not add Lake Creek as it seems like the source population is Magone Lake
#habdatdistinfo[habdatdistinfo$River_Tributary == 'Lake Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Lake Creek', 'interval'] + introdist-JDBeech + BeechEastForkBeech
habdatdistinfo[habdatdistinfo$River_Tributary == 'Lower mainstem' | 
                 habdatdistinfo$River_Tributary == 'Upper mainstem', 'Spread_dist'] <- abs(habdatdistinfo[habdatdistinfo$River_Tributary == 'Lower mainstem' |
                                                                                                            habdatdistinfo$River_Tributary == 'Upper mainstem', 'interval']-introdist)
habdatdistinfo[habdatdistinfo$River_Tributary == 'Murderers Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Murderers Creek', 'interval'] + SFMurderers + abs(JDSF-introdist)
habdatdistinfo[habdatdistinfo$River_Tributary == 'North Fork', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'North Fork', 'interval'] + abs(JDNF-introdist)
habdatdistinfo[habdatdistinfo$River_Tributary == 'Rock Creek', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'Rock Creek', 'interval'] + abs(JDRock-introdist)
habdatdistinfo[habdatdistinfo$River_Tributary == 'South Fork', 'Spread_dist'] <- habdatdistinfo[habdatdistinfo$River_Tributary == 'South Fork', 'interval'] + abs(JDSF-introdist)

############## Get network COMID (see end of script Crayfish_model_visualization for joining of sites with network) #######
sitesnetjoin <- read.dbf(file.path(datadir,'Sampled_sites_notes_CPUE_kick_edges_18_join_1.dbf'))
#Original path: "F:\Chapter2_HexSim_Crayfish\src\Crayfish_model\Parameters_and_outputs\Network_23\Network23_test16_2025\Sampled_sites_notes_CPUE_kick_edges_18_join_1.dbf"
sitesnetjoin <- sitesnetjoin[,c('Site_ID', 'COMID','rid')]
colnames(sitesnetjoin) <- c('Site','COMID','rid')
############## Get modeled stream temperature ######################
temp <- read.csv(file.path(datadir,'tempdat_9920_mean.csv'))
#Original path: "F:\Chapter2_HexSim_Crayfish\data\ISEMP\HexSim_ready\tempdat_9920_mean.csv"
colnames(temp)
habdatdistinfo <- merge(habdatdistinfo, sitesnetjoin, by='Site', all.x=T)
habdatdistinfo[habdatdistinfo$Site == 43, 'rid'] <- 676 #Nudge Murderer's Creek site a little downstream to get invasion time estimate from model
habdatdistinfo <- merge(habdatdistinfo, temp[,c(1,196:207)], by='rid', all.x=T, all.y=F)
#Degree days: sum over months of matrix multiplication of mean temperature per month by vector of number of days per month
habdatdistinfo$degdays <- rowSums(t(t(habdatdistinfo[,grepl('X20', colnames(habdatdistinfo))])*c(31,30,31,30,31,31,29,31,30,31,30,31)))

temp[temp$rid == 1383 | temp$rid == 1704,]

############## Get NHD monthly estimated discharge data ############
networkflow <- read.dbf(file.path(datadir,'NHDplus_streams_HU6clip_format_UTM_positiveflow_perennial_McNysetpointjoin_correct_withfields_o025_nodry.dbf'))
#Original path: "F:\Chapter2_HexSim_Crayfish\data\River_network\NHDplus_streams_HU6clip_format_UTM_positiveflow_perennial_McNysetpointjoin_correct_withfields_o025_nodry.dbf"
habdatdistinfo <- merge(habdatdistinfo, networkflow[,c('COMID',colnames(networkflow)[grepl('Q0001E',colnames(networkflow))])], by='COMID', all.x=T,all.y=F)
# which(networkflow$COMID != networkflow$COMID_1)
# which(networkflow$COMID != networkflow$COMID_12) #Don't use COMID_12


############## Get average modeled first invasion date #####################
#HexSim output directory
root_dir <- file.path(datadir, 'Network23_test16_2025')
#Original path: "F:\Chapter2_HexSim_Crayfish\src\Crayfish_model\Network_23\Results\Network23_test16_2025"
#Wanted directory to write results
root_dir_out <- file.path(resdir, 'Network23_test16_2025')
#Original path: "F:\Chapter2_HexSim_Crayfish\src\Crayfish_model\Parameters_and_outputs\Network_23\Network23_test16_2025"
#List all directories in HexSim output directory, each corresponding to 
scenario_reps_list <- list.files(root_dir)


#List all directories in HexSim output directory, each corresponding to 
scenario_reps_list <- list.files(root_dir)
scenario_reps_list
#Compute first year that each site was estimated to be colonized by rusty crayfish
mininvyear <- lapply(scenario_reps_list, function(x) {
  scenario_rep <- paste(root_dir, x, "Generated Properties", sep = "/")
  scenario_reps_res_list <- as.vector(list.files(scenario_rep))
  scenario_reps_res<-paste(scenario_rep, scenario_reps_res_list[pmatch("Dens",scenario_reps_res_list)], sep = "/")
  temp_dens <- read.csv(scenario_reps_res)
  temp_densdat <- temp_dens[temp_dens$Reach %in% sitesnetjoin$rid,]
  mininvyr <- adply(temp_densdat[,2:314], 1, function(x) {
    minyr <- min(as.numeric(substr(colnames(x[which(x>0)]),6,10)))
    minyr[minyr==Inf] <- NA
    minyr
  })
  cbind(temp_densdat$Reach, mininvyr$V1)
})
#Convert list of dataframes to a single data frame
mininvyear_df <- ldply(mininvyear, data.frame)
#Compute statistics on the estimated first year of colonization of reach based on multiple replicates of model
mininvyear_df <- ddply(mininvyear_df[!is.na(mininvyear_df$X2),], .(X1), summarize, meaninvyr = (mean(X2, na.rm=T)/12), maxinvyr= 1999.5+(max(X2,na.rm=T)/12), mininvyr =1999.5+(min(X2,na.rm=T)/12), sdinvyr=sd(X2,na.rm=T), lengthinvyr=length(X2))
colnames(mininvyear_df)[1] <- 'rid'

#Get towns' location on the network
JD_towns <-data.frame(town=c('John day','Mount Vernon','Dayville','Spray'),town_dist=c(391271,378420,335231,268244))
habdatdistinfo[habdatdistinfo$Site >= 106 & habdatdistinfo$River_Tributary == 'Upper mainstem', 'meaninvyr'] <- 9



################################################# PUT CRAYFISH DATA TOGETHER ###############################################################
#Import crayfish individual data
craydat <- read.csv(file.path(datadir,"Fieldworkdata/Crayfish_4.csv"))

#Correct some quirks in Crayfish data (fill in NAs)
craydat <- craydat[!is.na(craydat$Cray_ID),]
craydat[craydat$Meso_type=='RuorRi' & !is.na(craydat$Meso_type),'Meso_type'] <- 'Ri-Ru'
craydat[craydat$Meso_type=='Po/Ru'& !is.na(craydat$Meso_type),'Meso_type'] <- 'Ru-Ro'
craydat[craydat$Meso_type=='Ru or Pool-edge'& !is.na(craydat$Meso_type),'Meso_type'] <- 'Pool-edge'
craydat[craydat$Site==27 & is.na(craydat$Meso_type),'Meso_type'] <- 'Ru'
craydat[craydat$Site==20 & is.na(craydat$Meso_type),'Meso_type'] <- 'Ru'
craydat[craydat$Site==8 & is.na(craydat$Meso_type),'Meso_type'] <- 'Ru'
craydat[craydat$Site==10 & is.na(craydat$Meso_type),'Meso_type'] <- 'Ru'


#Join with RNADNA data
RNADNA_dat <- read.csv(file.path(resdir,"RNADNA_analysis/RNADNA_datformat.csv"))
RNADNA_dat[RNADNA_dat$Site == 20 & RNADNA_dat$Cray_ID == 21, 'Cray_ID'] <- 20
RNADNA_params_format <- read.csv(file.path(resdir,"RNADNA_analysis/RNADNA_paramsformat.csv"))
craydat <-merge(craydat, RNADNA_dat, by=c('Site','Cray_ID'), all.x=T, all.y=T)
#Missing ones: 2/17 (digestion did not work), 113/8 (digestion did not work),2/6,27/27

#Join with isotope data
isodat <- read.csv(file.path(datadir,'Isotope_analysis/Messager JDR-Cray2016Trays1-4 0217_format.csv'))
isodat$Sample.ID <- as.character(isodat$Sample.ID)
isodat[isodat$Sample.ID == '35-D-O','Sample.ID'] <- '35-D-17'
isodat[isodat$Sample.ID == '15F-6','Sample.ID'] <- '15-F-6'
isodat[isodat$Sample.ID == '25-MI','Sample.ID'] <- '25-MI-3'
isodat[isodat$Sample.ID == '15-O-O','Sample.ID'] <- '15-O-18'
isodat[isodat$Sample.ID == '6-O-O-2','Sample.ID'] <- '6-O-17'
isodat[isodat$Sample.ID == '11-O-15','Sample.ID'] <- '10-O-16'
isodat[isodat$Sample.ID == '27-E-19','Sample.ID'] <- '107-E-19'
isodat[isodat$Sample.ID == '31-E-4','Sample.ID'] <- '31-E-41'
isodat[isodat$Sample.ID == '119-O-13','Sample.ID'] <- '113-O-13'
isodat[isodat$Sample.ID == '41-F-29','Sample.ID'] <- '41-F-27'

isodat$Site <- substr(isodat$Sample.ID, 1,3)
isodat[grepl('-', isodat$Site),"Site"] <- substr(isodat[grepl('-', isodat$Site),"Site"],1,2)
isodat[grepl('-', isodat$Site),"Site"] <- substr(isodat[grepl('-', isodat$Site),"Site"],1,1)
isodat$Site <- as.numeric(isodat$Site)

isodat$Cray_ID <- substr(isodat$Sample.ID, nchar(as.character(isodat$Sample.ID))-1,nchar(as.character(isodat$Sample.ID)))
isodat[grepl('-', isodat$Cray_ID),"Cray_ID"] <- substr(isodat[grepl('-', isodat$Cray_ID),"Cray_ID"],2,2)
isodat[grepl('AND', isodat$Sample.ID),"Cray_ID"] <- 1

craydat_testclean <-merge(craydat, isodat[isodat$Type.of.Material == 'Crayfish muscle',], by=c('Site','Cray_ID'), all.x=T, all.y=T)
#Not even written in oven list: 27-11, 8-12
#Confirmed missing from oven:41-11 (never made it), 31-35 (insufficient material),31-40 (insufficient material),107-8 (insufficient material),107-10 (insufficient material)
#Confirmed missing from tinning:41-16,4-12
#Check around plate 3 well E-6
#Check around plate 1 well F-2

#Check macroinvertebrate isotopic data
ggplot(isodat[isodat$Type.of.Material == 'Whole macroinvert',], aes(x=Site, y=d15N, color=factor(Site))) + geom_point()
ggplot(isodat, aes(x=Site, y=d15N, color=factor(Type.of.Material))) + geom_point()
#Check which sites don't have macroinv ref data
unique(isodat$Site)[which(!(unique(isodat$Site) %in% unique(isodat[isodat$Type.of.Material == 'Whole macroinvert',"Site"])))]
#MI isotope for 107 missing as did not have enough material, others were never collected in the field

#Compute statistics for d15N of macroinvertebrates
isodat_MI <- ddply(isodat[isodat$Type.of.Material == 'Whole macroinvert',], .(Site), summarize, 
      MId15N_mean=mean(d15N),MId15N_max=max(d15N),MId15N_min=min(d15N),MId15N_max=max(d15N))
#Merge with crayfish d15N data
isodatcray <- merge(isodat[isodat$Type.of.Material=='Crayfish muscle',], isodat_MI, by='Site', all.x=T,all.y=T)
#Merge with site info to understand potential drivers of base d15N
habdatdistinfo_MI <- merge(habdatdistinfo, isodat_MI, by='Site', all.x=T)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=greenmean, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=cyanomean, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=diatommean, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=greenmean+cyanomean+diatommean, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=Q0001E_MA, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=degdays, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=degdays, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=X2016.07.01, y=MId15N_mean, color=factor(River_Tributary))) + geom_point()+ geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ggplot(habdatdistinfo_MI[!is.na(habdatdistinfo_MI$MId15N_mean),], aes(x=interval, y=MId15N_mean, color=factor(River_Tributary))) + geom_point() + geom_text(aes(label = Site, size = 3), check_overlap = TRUE)
# ###Could not identify obvious drivers of baseline d15N so:

#Use site 35 to set baseline d15N for site 107 (as nearest site on mainstem)
isodatcray[isodatcray$Site == 107,'MId15N_mean'] <-  isodatcray[isodatcray$Site == 35,'MId15N_mean'][1]
#Use average of site 4 and 8 to set baseline d15N for site 6
isodatcray[isodatcray$Site == 6,'MId15N_mean'] <-  mean(c(isodatcray[isodatcray$Site == 4,'MId15N_mean'][1], isodatcray[isodatcray$Site == 8,'MId15N_mean'][1]))
#Use median baseline d15N for site 46, 49, 52
qplot(isodat_MI$MId15N_mean) + geom_vline(xintercept=median(isodat_MI$MId15N_mean), color='red') + geom_vline(xintercept=mean(isodat_MI$MId15N_mean), color='blue')
isodatcray[isodatcray$Site == 46|isodatcray$Site == 49 |isodatcray$Site == 52,'MId15N_mean'] <- median(isodat_MI$MId15N_mean)

#Compute crayfish trophic level (based on discrimination of 2.54 from Glon et al. 2.54)
isodatcray$trophiclevel <- 2+((isodatcray$d15N-isodatcray$MId15N_mean)/2.54)

craydat <-merge(craydat, isodatcray[,c('Site','Cray_ID','d15N','Namount_ug','MId15N_mean','trophiclevel')], 
                by=c('Site','Cray_ID'), all.x=T, all.y=T)

str(craydat)
craydat$Weight <- as.numeric(as.character(craydat$Weight))
craydat$CL <- as.numeric(as.character(craydat$CL))
craydat$Chelae_L <- as.numeric(as.character(craydat$Chelae_L))
craydat <- craydat[craydat$Species=='OR',]
craydat <- craydat[!is.na(craydat$CL),]
craydat <- craydat[craydat$Meso_type!='Lake',]
craydat <- merge(craydat, habdatdistinfo, by="Site", all.x=T)
craydat[craydat$greenmean == 0 & !is.na(craydat$greenmean),'greenmean'] <- 0.05 #To allow the log of greenmean as very log-normally distributed
craydat <- craydat[craydat$Site!=47 & craydat$Site!=48,] #Take out sites from East Beech Creek Fork as these are crayfish that dispersed from Magone Lake -- different population
craydat$Spread_dist <- craydat$Spread_dist/1000
craydat[!is.na(craydat$Miss_App),'missap_binary'] <- 'Y'
craydat[is.na(craydat$Miss_App),'missap_binary'] <- 'N'
#colnames(craydat)
#summary(craydat)

#GENERAL DESCRIPTION OF VARIABLES 
#stat.desc(craydat)
#hist.plots(craydat[,sapply(craydat, is.numeric)])
#Log: CL, Chelae_L,Weight, RNADNAratio, sqrt(greenmean), sqrt(cyanomean)?,sqrt(diatommean),log(AFDWmean),log(Kick_mean),log(Q0001E_MA),
#qqnorm.plots(craydat[,sapply(craydat, is.numeric)])

################################################# COMPUTE RESIDUALS FOR WEIGHT AND CHELA LENGTH #######################
###############Estimate relationship between CL and Weight
#By linear regression
Weight_CL_lm <-lm(log(Weight)~log(CL), data=craydat[!is.na(craydat$Weight) & is.na(craydat$Miss_App),])
summary(Weight_CL_lm)
mean(abs(craydat[!is.na(craydat$Weight) & is.na(craydat$Miss_App),'Weight'] - exp(predict(Weight_CL_lm))))
par(mfrow=c(2,2))
plot(Weight_CL_lm)
par(mfrow=c(1,1))
qplot(craydat$Weight_CL_lm_res)
craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),'CL_weight_lm_pred'] <- predict(Weight_CL_lm) 
craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),'CL_weight_lm_res'] <- residuals(Weight_CL_lm)
CL_Weight_loglog <- ggplot(craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),], aes(x = CL, y = Weight)) +
  geom_point() + 
  #geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  #geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) +
  geom_line(aes(y=exp(CL_weight_lm_pred)), color = "#225ea8", size = 1.25) +
  scale_x_log10(breaks=c(5,10,20,30,40,50)) +
  scale_y_log10(breaks=c(0.1,1,10,25,50))+
  theme_bw()
ggsave('cray_Weight0_CL_loglog.png',CL_Weight_loglog)

#By MLE
CL_Weight_mod <- function(a, b, sigma, data = craydat[!is.na(craydat$Weight) & is.na(craydat$Miss_App),]) {
  CL <- data$CL
  W <- data$Weight
  # Fit the equation with the parameters
  model.predR <- a*CL^b
  # Calculates the negative log likelihood 
  NLL <- -sum(dnorm(x= W,mean=model.predR,sd=sigma,log=TRUE))
  return(NLL)
}
#Estimate the parameters that will minimize the negative log-likeligood of the model fit 
MLE.res.R <- mle2(minuslogl=CL_Weight_mod, 
                  start=list(a=0.2, b=1.5, sigma=5), 
                  method="Nelder-Mead",
                  control = list(maxit = 10000))
MLE.res.R 

#By nls
Weight_CL_nls <- nls(Weight~b*CL^c, data=craydat[!is.na(craydat$Weight) & is.na(craydat$Miss_App),], start=c(b=0.1,c=3), trace=T)
summary(Weight_CL_nls)
Weight_CL_nlsfun <- function(CL) {log10(0.0002243*CL^3.11812)}
craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),'CL_weight_pred'] <- predict(Weight_CL_nls) 
craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),'CL_weight_nls_res'] <- residuals(Weight_CL_nls)
qplot(log10(residuals(Weight_CL_nls)))
plot(Weight_CL_nls)

CL_Weight_nls <- ggplot(craydat[!is.na(craydat$Weight)& is.na(craydat$Miss_App),], aes(x = CL, y = Weight)) +
  geom_point() + 
  #geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  #geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) +
  stat_function(fun = Weight_CL_nlsfun, color = "#225ea8", size = 1.25) +
  scale_x_log10(breaks=c(5,10,20,30,40,50)) +
  scale_y_log10(breaks=c(0.1,1,10,25,50))+
  theme_bw()
ggsave('cray_Weight0_CL_nls.png',CL_Weight_nls)

#################Estimate relationship between CL and chela length
#With log-log model
ChelaL_CL_lm <-lm(log(Chelae_L)~log(CL), data=craydat[!is.na(craydat$Chelae_L),])
summary(ChelaL_CL_lm)
plot(ChelaL_CL_lm)

craydat[!is.na(craydat$Chelae_L),'CL_ChelaL_lm_pred'] <- predict(ChelaL_CL_lm) 
craydat[!is.na(craydat$Chelae_L),'CL_ChelaL_lm_res'] <- residuals(ChelaL_CL_lm)
qplot(residuals(ChelaL_CL_lm))

chelaratio_CL_trib_loglog <- ggplot(craydat[!is.na(craydat$Chelae_L)& is.na(craydat$Miss_App),], aes(x = CL, y = Chelae_L)) +
  geom_jitter(size=3, alpha=1/2,aes(color=factor(River_Tributary))) +
  #geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  #geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) +
  geom_line(aes(y=CL_ChelaL_lm_res), color = "#225ea8", size = 1.25) + 
  scale_x_log10(breaks=c(5,10,25,50), labels=c(5,10,25,50))+
  scale_y_log10(breaks=c(5,10,25,50))+
  theme_bw()
chelaratio_CL_trib_loglog
ggsave('cray_chelaratio0CL_loglog_trib.png', plot=chelaratio_CL_trib_loglog)

#With Non-linear least-square model
ChelaL_CL_nls <- nls(Chelae_L~b*CL^c, data=craydat[!is.na(craydat$Chelae_L),], start=c(b=0.2,c=1.5), trace=T)
summary(ChelaL_CL_nls)
plot(ChelaL_CL_nls)
mean(abs(craydat[!is.na(craydat$Chelae_L),'Chelae_L'] - predict(ChelaL_CL_nls)))

ChelaL_CL_nlsfun <- function(CL) {log10(0.1733*CL^1.4653)}
craydat[!is.na(craydat$Chelae_L),'CL_ChelaL_nls_pred'] <- predict(ChelaL_CL_nls) 
craydat[!is.na(craydat$Chelae_L),'CL_ChelaL_nls_res'] <- residuals(ChelaL_CL_nls)
qplot(residuals(ChelaL_CL_nls))

chelaratio_CL_trib_nls <- ggplot(craydat, aes(x = CL, y = Chelae_L)) +
  geom_jitter(size=3, alpha=1/2,aes(color=factor(River_Tributary))) +
  #geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  #geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) +
  stat_function(fun = ChelaL_CL_nlsfun, color = "#225ea8", size = 1.25) + 
  scale_x_log10(breaks=c(5,10,25,50), labels=c(5,10,25,50))+
  scale_y_log10(breaks=c(5,10,25,50))+
  theme_bw()
chelaratio_CL_trib_nls
ggsave('cray_chelaratio0CL_nls_trib.png', plot=chelaratio_CL_trib_nls)

#######
craydatOR <- craydat[craydat$Species=='OR' & (craydat$Method=='Kick' | craydat$Method=='Snorkel/Dipnet' | craydat$Method == 'Dipnet'),]
unique(craydatOR$Site)
########
################################################# FORMAT BINNED LENGTH DATA ############################################
CPUEOR$CRUSHED.MISSED <- as.numeric(as.character(CPUEOR$CRUSHED.MISSED))
CPUEOR[is.na(CPUEOR$CRUSHED.MISSED),'CRUSHED.MISSED'] <- 0
CPUEOR$u15_rel <- CPUEOR$u15/(CPUEOR$CPUE-CPUEOR$CRUSHED.MISSED)
CPUEOR$o15u25_rel <- CPUEOR$o15u25/(CPUEOR$CPUE-CPUEOR$CRUSHED.MISSED)
CPUEOR$o25_rel <- CPUEOR$o25/(CPUEOR$CPUE-CPUEOR$CRUSHED.MISSED)
CPUEOR$reltest <- CPUEOR$u15_rel+CPUEOR$o15u25_rel+CPUEOR$o25_rel
CPUEOR[is.nan(CPUEOR$u15_rel),'u15_rel'] <- 0
CPUEOR[is.nan(CPUEOR$o15u25_rel),'o15u25_rel'] <- 0
CPUEOR[is.nan(CPUEOR$o25_rel),'o25_rel'] <- 0

CPUEOR_stat <- CPUEOR[CPUEOR$CPUE!=0 & CPUEOR$CPUE!=CPUEOR$CRUSHED.MISSED,]
CPUEOR_statfig <- ddply(CPUEOR_stat, .(Method, Site), summarize, 
                        u15_relmean =mean(u15_rel), o15u25_relmean=mean(o15u25_rel), o25_relmean=mean(o25_rel), totalcatch = sum(CPUE))
CPUEOR_habdat <- merge(CPUEOR_statfig, habdatdistinfo[,c('Site', 'meaninvyr','Spread_dist', 'River_Tributary')], by='Site', all.x=T)
CPUEOR_habdat <- melt(CPUEOR_habdat[CPUEOR_habdat$Method=='Kick',], 
                      id.vars=c('Site', 'meaninvyr','Spread_dist', 'River_Tributary','totalcatch'), 
                      measure.vars=c('u15_relmean','o15u25_relmean','o25_relmean'))

CPUEOR_statmod <- ddply(CPUEOR_stat, .(Site), summarize, 
                        u15_rel =sum(u15)/(sum(CPUE)-sum(CRUSHED.MISSED)), 
                        o15u25_rel=sum(o15u25)/(sum(CPUE)-sum(CRUSHED.MISSED)), 
                        o25_rel=sum(o25)/(sum(CPUE)-sum(CRUSHED.MISSED)))

################################################# FORMAT SITE-LEVEL CRAYFISH DATA ###############################################################
craydat_stat <- ddply(craydatOR, .(Site), function(x) {
  ncray <- nrow(x)
  sexratio <- nrow(x[x$Sex=='M',])/(nrow(x[x$Sex=='M',])+nrow(x[x$Sex=='F',]))
  males<- nrow(x[x$Sex=='M',])
  females <- nrow(x[x$Sex=='F',])
  Miss_App_tot <- nrow(x[!is.na(x$Miss_App),])/ncray
  missap <- nrow(x[!is.na(x$Miss_App),])
  nomissap <- nrow(x[is.na(x$Miss_App),])
  CL_mean <- mean(x$CL,na.rm=T)
  CL_sd <- sd(x$CL,na.rm=T)
  CL_max <- max(x$CL, na.rm=T)
  CL_qt90 <- quantile(x$CL, probs=0.90,na.rm=T)
  CL_qt95 <- quantile(x$CL, probs=0.95,na.rm=T)
  CL_SE <- sd(x$CL,na.rm=T)/sqrt(nrow(x[!is.na(x$CL),]))
  Chelae_L_mean <- mean(x$Chelae_L,na.rm=T)
  Chelae_L_SE <- sd(x$Chelae_L,na.rm=T)/sqrt(nrow(x[!is.na(x$Chelae_L),]))
  Chelae_res_mean <-  mean(x$CL_ChelaL_nls_res, na.rm=T)
  Chelae_res_se <- sd(x$CL_ChelaL_nls_res,na.rm=T)/sqrt(nrow(x[!is.na(x$CL_ChelaL_nls_res),]))
  molting <- nrow(x[x$Molting=='Y',])/ncray
  RNADNAratio_mean <- mean(x$RNADNAratio,na.rm=T)
  RNADNAratio_SE <- sd(x$RNADNAratio,na.rm=T)/sqrt(nrow(x[!is.na(x$RNADNAratio),]))
  trophiclevel_mean <- mean(x$trophiclevel,na.rm=T)
  trophiclevel_SE <- sd(x$trophiclevel,na.rm=T)/sqrt(nrow(x[!is.na(x$trophiclevel),]))
  data.frame(ncray, sexratio, males, females, Miss_App_tot,missap, nomissap, CL_mean, CL_sd, CL_max, 
             CL_qt90, CL_qt95, CL_SE, Chelae_L_mean, Chelae_L_SE, Chelae_res_mean, Chelae_res_se, 
             molting, RNADNAratio_mean, RNADNAratio_SE,trophiclevel_mean, trophiclevel_SE)
})

#Process weight statistics differently than the rest as need to only compute it on crayfish that have both chelae
craydat_stat_weight <- ddply(craydatOR[is.na(craydatOR$Miss_App),], .(Site), function(x) {
  Weight_mean <- mean(x$Weight,na.rm=T)
  Weight_SE <- sd(x$Weight,na.rm=T)/sqrt(nrow(x[!is.na(x$Weight),]))
  Weight_res_mean <-  mean(x$CL_weight_lm_res, na.rm=T)
  Weight_res_se <- sd(x$CL_weight_lm_res,na.rm=T)/sqrt(nrow(x[!is.na(x$CL_weight_lm_res),]))
  data.frame(Weight_mean, Weight_SE, Weight_res_mean, Weight_res_se)
}
)

craydat_stat <- merge(craydat_stat, craydat_stat_weight, by="Site", all.x=T)
craydat_stat <- merge(craydat_stat, habdatdistinfo, by="Site", all.x=T)
craydat_stat <- merge(craydat_stat, CPUEOR_statmod, by="Site", all.x=T)
craydat_stat$interval <- craydat_stat$interval/1000
craydat_stat$Spread_dist <- craydat_stat$Spread_dist/1000
craydat_stat$River_Tributary <- as.character(craydat_stat$River_Tributary)
craydat_stat <- craydat_stat[!grepl('Creek',craydat_stat$River_Tributary) & craydat_stat$Site != 106,]
craydat_stat[craydat_stat$Site == 36 | craydat_stat$Site == 107,'River_Tributary'] <- 'Upst Upper mainstem'
craydat_stat[craydat_stat$Site == 35, 'Spread_dist'] <- 0

#To be able to match AFDWmean and RNADNAratio as much as possible. Assign AFDW from site 26 to site 27.
craydat_stat[craydat_stat$Site == 27, 'AFDWmean'] <- craydat_stat[craydat_stat$Site == 26, 'AFDWmean']

#stat.desc(craydat_stat)
#hist.plots(craydat_stat[,sapply(craydat_stat, is.numeric)])
#Transform: sqrt(Miss_App_tot), log(CL_mean), log(CL_mean +- CL_sd), log(Chealae_L_mean), log(Weight_mean), sqrt(weight_mean), 
# sqrt(RNADNAratio_mean), sqrt(greenmean), sqrt(cyanomean)?,sqrt(diatommean),log(AFDWmean),log(Kick_mean),log(Q0001E_MA)

##########################################################################################################################################################
setwd(file.path(resdir,"Dataexplo/"))
#Total number of crayfish collected at sites used in the analysis
sum(CPUEOR[CPUEOR$Site %in% craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'Site'],'CPUE'])

#Number of sampled crayfish in sites used in the analysis
craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'Site']
length(craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'Site'])
sum(craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'ncray'])

#Number of sampled crayfish with trophic position and growth data
nrow(craydat[craydat$Site %in% craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'Site'] &
          !is.na(craydat$trophiclevel),])
nrow(craydat[craydat$Site %in% craydat_stat[craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean), 'Site'] &
               !is.na(craydat$RNADNAratio),])

################################################################################## 2. VISUALIZE DATA #####################################################
############################################################################################
################################################# HABITAT DATA #########################################
############## Primary productivity data########
#Upstream-downstream trend
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = interval, y = greenmean, color = River_Tributary)) + geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,4.5)) +
  geom_errorbar(aes(x=interval, ymin=greenmean-1.96*(greensd/sqrt(10)), ymax=greenmean+1.96*(greensd/sqrt(10)))) 
  #geom_vline(xintercept=JD_towns$town_dist, col='red')
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = interval, y = cyanomean, color = River_Tributary)) + geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,4.5))+
  geom_errorbar(aes(x=interval, ymin= cyanomean-1.96*(cyanosd/sqrt(10)), ymax=cyanomean+1.96*(cyanosd/sqrt(10)))) 
  #geom_vline(xintercept=JD_towns$town_dist, col='red')
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = interval, y = diatommean, color = River_Tributary)) + geom_point(size = 4)+ geom_smooth(span=1, method='gam') +
  theme_bw() + scale_x_sqrt() + coord_cartesian(ylim=c(0,4.5))+
  geom_errorbar(aes(x=interval, ymin= diatommean-1.96*(diatomsd/sqrt(10)), ymax=diatommean+1.96*(diatomsd/sqrt(10)))) 
  #geom_vline(xintercept=JD_towns$town_dist, col='red')
grid.arrange(g,c,d, ncol = 3)
ggsave('habitat_1PP1_spreaddist.png', plot=grid.arrange(g,c,d, ncol = 3), width = 20, height = 12, units='in')

#Invasion distance trend
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Spread_dist, y = greenmean, color = River_Tributary)) + geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Spread_dist, y = cyanomean, color = River_Tributary))+ geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Spread_dist, y = diatommean, color = River_Tributary)) + geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Spread_dist, y = greenmean+cyanomean+diatommean, color = River_Tributary)) + geom_point(size = 4) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))

#Temperature
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = degdays, y = greenmean)) + geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = degdays, y = cyanomean))+ geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = degdays, y = diatommean)) + geom_point(size = 4,aes(color = River_Tributary))+ geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

#Discharge
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Q0001E_8, y = greenmean)) + geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_log10() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Q0001E_8, y = cyanomean))+ geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_log10() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Q0001E_8, y = diatommean)) + geom_point(size = 4,aes(color = River_Tributary))+ geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_log10() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

#Velocity
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = velmean, y = greenmean)) + geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = velmean, y = cyanomean))+ geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = velmean, y = diatommean)) + geom_point(size = 4,aes(color = River_Tributary))+ geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

#Depth
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = depthmean, y = greenmean)) + geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = depthmean, y = cyanomean))+ geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = depthmean, y = diatommean)) + geom_point(size = 4,aes(color = River_Tributary))+ geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

#Crayfish density
g <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Kick_mean, y = greenmean)) + geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
c <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Kick_mean, y = cyanomean))+ geom_point(size = 4,aes(color = River_Tributary)) + geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
d <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$greenmean),], aes(x = Kick_mean, y = diatommean)) + geom_point(size = 4,aes(color = River_Tributary))+ geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") + scale_x_sqrt() + coord_cartesian(ylim=c(0,450))
grid.arrange(g,c,d, ncol = 3)

############## Macroinvertebrate productivity ##################
ggplot(habdatdistinfo, aes(x = Spread_dist, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=Spread_dist, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  geom_text(aes(label=Site), color='black')+
  coord_cartesian(ylim=c(0,1)) + 
  geom_smooth(method='glm') +
  scale_y_continuous(expand=c(0,0))

macrointerval <-ggplot(habdatdistinfo, aes(x = interval, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=interval, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_smooth(method='glm') +
  geom_text(aes(label=Site), color='black') +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw()
macrointerval 
ggsave('habitat_2MACRO_interval.png', plot=macrointerval, width = 20, height = 12, units='in')

#Compare primary and macroinvertebrate prod (coord_cartesian allows to zoom in without losing error bars)
gmac <- ggplot(habdatdistinfo, aes(x = greenmean, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) +  geom_errorbar(aes(x=greenmean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  coord_cartesian(ylim=c(0,1)) +scale_y_continuous(expand=c(0,0))+  geom_smooth(span=1, method='gam') +
  theme_bw() +theme(legend.position = "none") 
cmac <- ggplot(habdatdistinfo, aes(x = cyanomean, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) +  geom_errorbar(aes(x=cyanomean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  coord_cartesian(ylim=c(0,1)) +  scale_y_continuous(expand=c(0,0))+  geom_smooth(span=1, method='gam') +
  theme_bw() +  theme(legend.position = "none") 
dmac <- ggplot(habdatdistinfo, aes(x = diatommean, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) +  geom_errorbar(aes(x=diatommean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  coord_cartesian(ylim=c(0,1)) +scale_y_continuous(expand=c(0,0))+geom_smooth(span=1, method='gam') +
  theme_bw() + theme(legend.position = "none") 
grid.arrange(gmac,cmac,dmac, ncol = 3)

ggplot(habdatdistinfo, aes(x = diatommean+greenmean+cyanomean, y = DWmean, color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=diatommean+greenmean+cyanomean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  coord_cartesian(ylim=c(0,1)) + 
  scale_y_continuous(expand=c(0,0))

#With temperature
ggplot(habdatdistinfo, aes(x = degdays, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=degdays, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  geom_text(aes(label=Site), color='black')+
  coord_cartesian(ylim=c(0,1)) + 
  geom_smooth(method='glm') +
  scale_y_continuous(expand=c(0,0))

#With water velocity
ggplot(habdatdistinfo, aes(x = velmean, y = AFDWmean, color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=velmean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  geom_text(aes(label=Site), color='black')+
  coord_cartesian(ylim=c(0,1)) + 
  scale_y_continuous(expand=c(0,0))

#With crayfish density
ggplot(droplevels(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean),]), aes(x = Kick_mean, y = AFDWmean,color = River_Tributary)) +
  geom_point(size=5) + 
  geom_errorbar(aes(x=Kick_mean, ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn)))) + 
  geom_text(aes(label=Site), color='black')+
  coord_cartesian(ylim=c(0,1)) + 
  geom_smooth(method='glm') +
  scale_y_continuous(expand=c(0,0))

qplot(habdatdistinfo$Kick_mean)
qplot(log(habdatdistinfo$Kick_mean+0.1))
qplot(habdatdistinfo$AFDWmean)

macro_dens <- lmer(AFDWmean~sqrt(Kick_mean), data=habdatdistinfo)
summary(macro_dens)

############## Physical habitat ########
tem <- ggplot(habdatdistinfo, aes(x = interval, y = tempmean, color = River_Tributary)) + geom_point(size = 2) + 
  geom_errorbar(aes(x=interval, ymin=tempmean-1.96*(tempsd/sqrt(10)), ymax=tempmean+1.96*(tempsd/sqrt(10))), alpha=1/2) +
  theme_bw() + theme(legend.position='none') + geom_smooth(method='lm') + scale_x_sqrt()
dep <- ggplot(habdatdistinfo, aes(x = interval, y = depthmean, color = River_Tributary)) + geom_point(size = 2)+ 
  geom_errorbar(aes(x=interval, ymin=depthmean-1.96*(depthsd/sqrt(10)), ymax=depthmean+1.96*(depthsd/sqrt(10))), alpha=1/2) +
  theme_bw() + theme(legend.position='none')+ geom_smooth(method='lm')+ scale_x_sqrt()
vel <- ggplot(habdatdistinfo, aes(x = interval, y = velmean, color = River_Tributary)) + geom_point(size = 2)+ 
  geom_errorbar(aes(x=interval, ymin=velmean-1.96*(velsd/sqrt(10)), ymax=velmean+1.96*(velsd/sqrt(10))), alpha=1/2) +
  theme_bw() + theme(legend.position='none') + geom_smooth(method='lm')+ scale_x_sqrt()
width <- ggplot(habdatdistinfo, aes(x = interval, y = Width, color = River_Tributary)) + geom_point(size = 2)+ 
  geom_smooth(method='lm')+ scale_x_sqrt()+ theme_bw()
grid.arrange(tem,dep, vel,width, ncol = 4)
ggsave('habitat_3physhab_interval.png', plot=grid.arrange(tem,dep, vel,width, ncol = 4), width = 20, height = 12, units='in')


tem <- ggplot(habdatdistinfo, aes(x = Spread_dist, y = tempmean, color = River_Tributary)) + geom_point(size = 2) + theme(legend.position='none') + geom_smooth(method='lm')
dep <- ggplot(habdatdistinfo, aes(x = Spread_dist, y = depthmean, color = River_Tributary)) + geom_point(size = 2)+ theme(legend.position='none')+ geom_smooth(method='lm')
vel <- ggplot(habdatdistinfo, aes(x = Spread_dist, y = velmean, color = River_Tributary)) + geom_point(size = 2)+ geom_smooth(method='lm')
grid.arrange(tem,dep, vel, ncol = 3)

#Compare measured temp to modeled temp
tempmean_Julytempmod<-ggplot(habdatdistinfo, aes(x = X2016.07.01, y = tempmean, color = River_Tributary)) + geom_point(size = 2) + 
  theme_bw() + coord_fixed() + geom_smooth(method='lm') + scale_y_continuous(limits=c(13,27)) +
  scale_y_continuous(limits=c(14, 27), extend=c(0,0))
tempmean_Julytempmod
ggsave('habitat_4julytempmod_tempmean.png', plot=tempmean_Julytempmod, width = 10, height = 12, units='in')

ggplot(habdatdistinfo, aes(x = Spread_dist, y = X2016.07.01, color = River_Tributary)) + geom_point(size = 2)
habitat_5degdays_spreaddist <- ggplot(habdatdistinfo, aes(x = Spread_dist, y = degdays, color = River_Tributary)) + geom_point(size = 3) + geom_smooth(span=1) + theme_classic()
ggsave('habitat_5degdays_spreaddist.png', plot=habitat_5degdays_spreaddist, width = 20, height = 12, units='in')


#Compare velocity and depth to discharge
ggplot(habdatdistinfo, aes(x = Q0001E_MA, y = depthmean, color = River_Tributary)) + geom_point(size = 2)
ggplot(habdatdistinfo, aes(x = Q0001E_8, y = depthmean, color = River_Tributary)) + geom_point(size = 2)

ggplot(habdatdistinfo, aes(x = Q0001E_MA, y = velmean, color = River_Tributary)) + geom_point(size = 2)
ggplot(habdatdistinfo, aes(x = Q0001E_8, y = velmean, color = River_Tributary)) + geom_point(size = 2)

#Check whether there is bias in mesohabitat sampled across invasion gradient
craydat_meso <-ddply(craydat, .(Site, Meso_type, Spread_dist), function(x) nrow(x))
craydat_mesoplot <-ggplot(craydat_meso, aes(x=factor(Site),y=V1)) + geom_bar(aes(fill=Meso_type),stat='identity') + theme_bw() +
  scale_x_discrete(name='Count')
ggsave('site_trophiclevel0_mesohabitatspread.png', plot=craydat_mesoplot, width = 10, height = 12, units='in')

############################################################################################
################################################# CRAYFISH-LEVEL DATA
############## Binned length #########################
#With catch >= 5
site_CL990_spreadist_binnedCPUEo5 <- ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 5,], aes(x=Spread_dist, y=value, fill= variable, color=River_Tributary)) + 
  geom_col(size=1.5,width=2000,position='stack') + 
  geom_text(aes(label=Site)) +
  scale_fill_manual(values = c("black", "grey", "white")) + 
  theme_bw()
ggsave('site_CL990_spreadist_binnedCPUEo5.png', plot=site_CL990_spreadist_binnedCPUEo5, width = 20, height = 12, units='in')

#With catch >= 10
site_CL990_spreadist_binnedCPUEo10 <- ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 10,], aes(x=Spread_dist, y=value, fill= variable, color=River_Tributary)) + 
  geom_col(size=1.5,width=2000,position='stack') + 
  geom_text(aes(label=Site)) +
  scale_fill_manual(values = c("black", "grey", "white")) + 
  theme_bw()
ggsave('site_CL990_spreadist_binnedCPUEo10.png', plot=site_CL990_spreadist_binnedCPUEo10, width = 20, height = 12, units='in')


#Just CL >= 25
site_CL991_spreadist_binnedCPUEo10CLo25 <- ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 10 & CPUEOR_habdat$variable =='o25_relmean',], aes(x=Spread_dist, y=value, fill=River_Tributary)) + 
  geom_col(size=1.5,width=2000) + 
  geom_text(aes(label=totalcatch),nudge_y=0.005) +
  theme_bw()
ggsave('site_CL991_spreadist_binnedCPUEo10CLo25.png', plot=site_CL991_spreadist_binnedCPUEo10CLo25, width = 20, height = 12, units='in')

ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 10 & CPUEOR_habdat$variable =='o25_relmean',], aes(x=meaninvyr, y=value, fill=River_Tributary)) + 
  geom_col(size=1.5,width=0.25,alpha=0.75) + 
  geom_text(aes(label=Site)) +
  theme_bw()
#15<CL<25
site_CL992_spreadist_binnedCPUEo10CLo15u25 <- ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 10 & CPUEOR_habdat$variable =='o15u25_relmean',], aes(x=Spread_dist, y=value, fill=River_Tributary)) + 
  geom_col(size=1.5,width=2000) + 
  geom_text(aes(label=Site)) +
  theme_bw()
ggsave('site_CL992_spreadist_binnedCPUEo10CLo15u25.png', plot=site_CL992_spreadist_binnedCPUEo10CLo15u25, width = 20, height = 12, units='in')

#Just CL < 15
site_CL993_spreadist_binnedCPUEo10CLu15 <- ggplot(CPUEOR_habdat[CPUEOR_habdat$totalcatch >= 10 & CPUEOR_habdat$variable =='u15_relmean',], aes(x=Spread_dist, y=value, fill=River_Tributary)) + 
  geom_col(size=1.5,width=2000) + 
  geom_text(aes(label=Site)) +
  theme_bw()
ggsave('site_CL993_spreadist_binnedCPUEo10CLu15.png', plot=site_CL993_spreadist_binnedCPUEo10CLu15, width = 20, height = 12, units='in')

############## LENGTH ####################################################
#Check size distribution by capture method
CL_method <- ggplot(craydat[craydat$Species=='OR',],aes(x=CL, group=Method))+
  geom_histogram(aes(fill=Method),bins=20)+
  facet_wrap(~Method, scales="free_y", ncol=1)+
  theme_bw()
CL_method
ggsave('cray_CL0_method.png', plot=CL_method, width = 10, height = 12, units='in')

#Check size distribution by tributary
CL_tributary <- ggplot(craydat,aes(x=CL))+geom_histogram(bins=20)+facet_wrap(~River_Tributary, scales="free_y", ncol=)+theme_bw()
CL_tributary 
ggsave('cray_CL0_tributary.png', plot=CL_tributary, width = 20, height = 12, units='in')

#Check size distribution with kicknet
CL_kicktrib <- ggplot(craydatOR,aes(x=CL))+geom_histogram(bins=20)+facet_wrap(~River_Tributary, scales="free_y", ncol=)+theme_bw()
CL_kicktrib
ggsave('cray_CL0_tributarykick.png', plot=CL_kicktrib, width = 20, height = 12, units='in')

CL_kicksite <- ggplot(craydatOR,aes(x=CL))+geom_histogram(bins=20)+facet_wrap(~Site, scales="free_y", ncol=)+theme_bw()
CL_kicksite
ggsave('cray_CL0_sitekick.png', plot=CL_kicksite, width = 20, height = 12, units='in')

#Check size distribution with trap
ggplot(craydat[craydat$Method=='Trap',],aes(x=CL))+geom_histogram(bins=20)+facet_wrap(~Site, scales="free_y", ncol=)+theme_bw()
#Check size distribution with dipnet
ggplot(craydat[craydat$Method=='Dipnet',],aes(x=CL))+geom_histogram(bins=20)+facet_wrap(~Site, scales="free_y", ncol=)+theme_bw()

#Check which sites have data for which methods
ggplot(craydat, aes(x=factor(Site),y=Method, size=CL)) + geom_point()

#Look at number of crayfish caught and measured by method in each site
craydat_meth <-ddply(craydat, .(Site, Method, River_Tributary), function(x) nrow(x))
ggplot(craydat_meth, aes(x=factor(Site),y=V1)) + geom_bar(aes(fill=Method),stat='identity') + theme_bw()


#DISTANCE FROM INVASION AND LENGTH
CL_spreadist_trib <- ggplot(craydat, aes(x = Spread_dist, y = CL, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary), shape=Method), size=2, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
ggsave('cray_CL2_spreadist_trib.png', plot=CL_spreadist_trib)

CL_spreadist2 <- ggplot(craydatOR, aes(x = Spread_dist, y = CL)) + 
  geom_hex(binwidth = c(30,4)) +
  scale_fill_viridis() + 
  theme_bw()
CL_spreadist2
ggsave('cray_CL2_spreadist.png', plot=CL_spreadist2)

CL_spreadist_tribbox <- ggplot(craydatOR[order(craydatOR$Spread_dist),], aes(x = factor(as.integer(Spread_dist)), y = CL, group=factor(Site))) +
  geom_boxplot(aes(color=factor(River_Tributary), shape=Method), size=3, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) +
  theme_bw()
CL_spreadist_tribbox
ggsave('cray_CL2_spreadist_tribbox.png', plot=CL_spreadist_tribbox, width=20, height=12, units='in')


CL_meaninvyr_trib <-ggplot(craydatOR, aes(x = meaninvyr, y = CL)) +
  geom_point(aes(color=River_Tributary, shape=Method), size=2, alpha=1/2) +
  #geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm') + facet_wrap(~River_Tributary)
CL_meaninvyr_trib
ggsave('cray_CL4_meaninvyr_trib.png', plot=CL_meaninvyr_trib)

CL_interval_trib <- ggplot(craydat, aes(x =interval, y = CL, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary), shape=Method), size=2, alpha=1/2) + 
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))
ggsave('cray_CL5_interval_trib.png', plot=CL_interval_trib)

#PRODUCTIVITY AND CARAPACE LENGTH
CL_ADFW <- ggplot(craydatOR, aes(x = AFDWmean, y = CL, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2))
ggsave('cray_CL6_ADFW.png', plot=CL_ADFW)

CL_AFDW_trib <- ggplot(craydat, aes(x = AFDWmean, y = CL)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  #geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary) + geom_smooth(method='lm')
CL_AFDW_trib
ggsave('cray_CL7_AFDW_trib.png', plot=CL_AFDW_trib)

CL_DW_trib <- ggplot(craydatOR, aes(x = DWmean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('cray_CL8_DW_trib.png', plot=CL_DW_trib)

CL_green_trib <- ggplot(craydatOR, aes(x = greenmean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary)
CL_green_trib
ggsave('cray_CL9_green_trib.png', plot=CL_green_trib)

CL_cyano_trib <- ggplot(craydatOR, aes(x = cyanomean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()+ facet_wrap(~River_Tributary)
CL_cyano_trib
ggsave('cray_CL9_cyano_trib.png', plot=CL_cyano_trib)

CL_diatom_trib <- ggplot(craydatOR, aes(x = diatommean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary)
CL_diatom_trib
ggsave('cray_CL9_diatom_trib.png', plot=CL_diatom_trib)

CL_pp_trib <- ggplot(craydatOR, aes(x = greenmean+cyanomean+diatommean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
CL_pp_trib
ggsave('cray_CL9_pp_trib.png', plot=CL_pp_trib)

#COMPETITION AND CARAPACE LENGTH
CL_CPUE_trib <- ggplot(craydat[craydat$Method == 'Kick',], aes(x = Kick_mean, y = CL)) +
  geom_point(aes(color=River_Tributary, shape=Method), size=4, alpha=1/2)+ 
  #geom_text(data=craydat[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(aes(group=River_Tributary), method='glm') +
  facet_wrap(~River_Tributary, scales='free_x')
CL_CPUE_trib
ggsave('cray_CL90_CPUE_trib2.png', plot=CL_CPUE_trib)

#HABITAT AND CARAPACE LENGTH
#Temperature
CL_tempfield_trib <- ggplot(craydatOR, aes(x = tempmean, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) 
ggsave('cray_CL91_tempfield_trib.png', plot=CL_tempfield_trib)

CL_tempmod_trib <- ggplot(craydatOR, aes(x = degdays, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2))
ggsave('cray_CL92_tempmod_trib.png', plot=CL_tempmod_trib)

CL_tempmod_trib2 <- ggplot(craydatOR, aes(x = degdays, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2)) +
  facet_wrap(~River_Tributary) + geom_smooth(method='lm')
CL_tempmod_trib2
ggsave('cray_CL92_tempmod_trib2.png', plot=CL_tempmod_trib)

CL_tempmodAug_trib <- ggplot(craydatOR, aes(x = X2016.07.01, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2))
ggsave('cray_CL93_tempmodAug_trib.png', plot=CL_tempmodAug_trib)

#Discharge
CL_discharge_trib <- ggplot(craydatOR, aes(x = (Q0001E_MA), y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2))
ggsave('cray_CL94_discharge_trib.png', plot=CL_discharge_trib)

#Mesohabitat
CL_Meso_type_trib <- ggplot(craydatOR, aes(x = (Meso_type), y = CL)) +
  geom_boxplot(size=1, alpha=1/2)
CL_Meso_type_trib
ggsave('cray_CL95_Meso_type_trib.png', plot=CL_Meso_type_trib)

#OTHER INTRINSIC PROPERTIES
CL_RNADNA_trib <- ggplot(craydat, aes(x = RNADNAratio, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3)+ 
  geom_smooth(method='gam') + facet_wrap(~Site, scales="free_y")
CL_RNADNA_trib
ggsave('cray_CL96_RNADNA_trib.png', plot=CL_RNADNA_trib, width = 20, height = 12, units='in')

CL_trophiclevel_trib <- ggplot(craydatOR, aes(x = trophiclevel, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('cray_CL97_trophiclevel_trib.png', plot=CL_trophiclevel_trib)

CL_trophiclevel_site <- ggplot(craydat, aes(x = trophiclevel, y = CL)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) + facet_wrap(~Site, scales = "free_y") + geom_smooth(method='lm')
CL_trophiclevel_site 
ggsave('cray_CL97_trophiclevel_site.png', plot=CL_trophiclevel_site, width = 20, height = 12, units='in')

CL_sex_trib <- ggplot(craydatOR, aes(x = Sex, y = CL)) +
  geom_boxplot(size=1, alpha=1/2) + facet_wrap(~Site)
  #geom_text(aes(label = Site, size = 2), alpha=1/3) 
CL_sex_trib
ggsave('cray_CL98_sex_trib.png', plot=CL_sex_trib)



############## WEIGHT ###################################################
#CL VS WEIGHT
ggplot(craydat, aes(x = CL, y = Weight)) +
  geom_point() + 
  geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) 

CL_Weight_loglog <-ggplot(craydat, aes(x = CL, y = Weight)) +
  geom_point() + 
  scale_x_log10(breaks=c(5,10,20,30,40,50)) +
  scale_y_log10(breaks=c(0.1,1,10,25,50))+
  geom_smooth(method='lm') + 
  theme_bw()
  # geom_text(aes(label = Site, size = 3), check_overlap = TRUE) +
  # geom_text(aes(label = Cray_ID, size =3), hjust = 1.5, color = "red",check_overlap = TRUE) +
CL_Weight_loglog
ggsave('cray_Weight0_CL_loglog.png',CL_Weight_loglog)


#DISTANCE FROM INVASION AND WEIGHT 
spreadist_weight <- ggplot(craydat, aes(x = Spread_dist, y = Weight_CL_lm_res, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))+
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='lm') +
  scale_y_continuous(limits=c(-1.5,1.5))
spreadist_weight
ggsave('cray_weight1_spreadist.png', plot=spreadist_weight, width = 20, height = 12, units='in')

meaninvyr_weight <- ggplot(craydat, aes(x = meaninvyr, y = Weight_CL_lm_res)) +
  geom_point(, size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  geom_smooth(method='lm')
meaninvyr_weight
ggsave('cray_weight2_meaninvyr.png', plot=meaninvyr_weight, width = 10, height = 10, units='in')

meaninvyr_weight_trib2 <-ggplot(craydat, aes(x = meaninvyr, y = Weight_CL_lm_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) +
  #geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))+ 
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')+
  facet_wrap(~River_Tributary,scales='free_y') + 
  theme_bw()
meaninvyr_weight_trib2
ggsave('cray_weight3_meaninvyr_trib2.png', plot=meaninvyr_weight_trib2, width = 10, height = 10, units='in')

interval_weight_trib <- ggplot(craydat, aes(x =interval, y = Weight_CL_lm_res, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/4) + 
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')
ggsave('cray_weight4_interval_trib.png', plot=interval_weight_trib, width = 8, height = 6, units='in')

#PRODUCTIVITY AND WEIGHT
ADFW_weight_trib <- ggplot(craydat[!is.na(craydat$AFDWmean),], aes(x = AFDWmean, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary,drop = TRUE, scales='free_y') + geom_smooth(method='lm')
ggsave('cray_weight6_ADFW_trib.png', plot=ADFW_weight_trib, width = 8, height = 6, units='in')

pp_weight_trib <- ggplot(craydat, aes(x = greenmean+cyanomean+diatommean, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary, scales='free') + geom_smooth(method='lm')
ggsave('cray_weight8_pp_trib.png', plot=pp_weight_trib, width = 20, height = 12, units='in')

#COMPETITION AND WEIGHT
CPUE_weight_trib <- ggplot(craydat, aes(x = Kick_mean, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/4)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + geom_smooth(method='lm')
#CPUE_weight_trib
ggsave('cray_weight9_CPUE_trib.png', plot=CPUE_weight_trib, width = 10, height = 10, units='in')

CPUE_weight_trib2 <- ggplot(craydat, aes(x = Kick_mean, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/4)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + geom_smooth(method='lm') + facet_wrap(~River_Tributary, scales='free')
#CPUE_weight_trib
ggsave('cray_weight9_CPUE_trib2.png', plot=CPUE_weight_trib2, width = 10, height = 10, units='in')

#HABITAT AND WEIGHT
#Temperature
tempfield_weight_trib <- ggplot(craydat, aes(x = tempmean, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+
  geom_smooth(method='lm', aes(color=River_Tributary))
ggsave('cray_weight90_tempfield_trib.png', plot=tempfield_weight_trib)

tempmod_weight_trib <- ggplot(craydat, aes(x = degdays, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  geom_smooth()
ggsave('cray_weight91_tempmod_trib.png', plot=tempmod_weight_trib, width = 10, height = 10, units='in')

#Discharge
discharge_weight_trib <- ggplot(craydatOR, aes(x = Q0001E_MA, y = Weight)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydatOR[!duplicated(craydatOR$Site),], aes(label = Site, size = 2))
ggsave('cray_weight93_discharge_trib.png', plot=discharge_weight_trib)

#Mesohabitat
weight_Meso_type_trib <- ggplot(craydat, aes(x = (Meso_type), y = Weight_CL_lm_res)) +
  geom_boxplot(size=1, alpha=1/2)
weight_Meso_type_trib
ggsave('cray_Weight95_Meso_type_trib.png', plot=weight_Meso_type_trib)

#Method
weight_method_trib <- ggplot(craydat, aes(x = Method, y = Weight_CL_lm_res)) +
  geom_boxplot(size=1, alpha=1/2) +
  geom_hline(yintercept=0, color='red')
weight_method_trib
ggsave('cray_Weight97_method_trib.png', plot=weight_method_trib)


#OTHER INTRINSIC PROPERTIES
CL_weight_trib <- ggplot(craydat, aes(x = CL, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=2, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 1.5), alpha=1/3)+
  geom_smooth(method='lm') + facet_wrap(~Site)
CL_weight_trib
ggsave('cray_weight94_CL_trib.png', plot=CL_weight_trib, width = 10, height = 10, units='in')

RNADNA_weight_trib <- ggplot(craydat, aes(x = RNADNAratio, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=2, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 1.5), alpha=1/3)+
  geom_smooth(method='lm')
ggsave('cray_weight94_RNADNA_trib.png', plot=RNADNA_weight_trib, width = 10, height = 10, units='in')

RNADNA_weight_trib2 <- ggplot(craydat, aes(x = RNADNAratio, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=3, alpha=1/2)+ 
  facet_wrap(~Site, scales='free_y') +
  geom_smooth(method='lm')
ggsave('cray_weight94_RNADNA_trib2.png', plot=RNADNA_weight_trib2, width = 20, height = 12, units='in')

trophiclevel_weight_trib <- ggplot(craydat, aes(x = trophiclevel, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=3, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 1.5), alpha=1/3)+
  geom_smooth(method='lm')
ggsave('cray_weight95_trophiclevel_trib.png', plot=trophiclevel_weight_trib, width = 20, height = 12, units='in')

trophiclevel_weight_trib2 <- ggplot(craydat, aes(x = trophiclevel, y = Weight_CL_lm_res)) +
  geom_point(aes(color=River_Tributary), size=3, alpha=1/2)+ 
  facet_wrap(~Site, scales='free_y') +
  geom_smooth(method='lm')
ggsave('cray_weight95_trophiclevel_trib2.png', plot=trophiclevel_weight_trib2, width = 20, height = 12, units='in')

weight_sex_trib <- ggplot(craydat, aes(x = Sex, y = Weight_CL_lm_res)) +
  geom_boxplot(size=1, alpha=1/2)
#geom_text(aes(label = Site, size = 2), alpha=1/3) 
weight_sex_trib
ggsave('cray_weight98_sex_trib.png', plot=weight_sex_trib)


############## CHELA LENGTH ###################################################
#CL vs ChELA LENGTH
chelaratio_CL_trib <- ggplot(craydat, aes(x = CL, y = Chelae_L)) +
  geom_jitter(size=3, alpha=1/2,aes(color=factor(River_Tributary))) +
  geom_smooth(method='lm') + 
  theme_bw()+
  scale_x_log10(breaks=c(5,10,25,50), labels=c(5,10,25,50))+
  scale_y_log10(breaks=c(5,10,25,50))
chelaratio_CL_trib
ggsave('cray_chelaratio0CL_trib.png', plot=chelaratio_CL_trib)

#DISTANCE FROM INVASION AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_spreadist_trib <- ggplot(craydat, aes(x = Spread_dist, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(method='lm')
chelaratio_spreadist_trib
ggsave('cray_chelaratio2_spreadist_trib.png', plot=chelaratio_spreadist_trib, width=20, height=12, units='in')

chelaratio_spreadist_trib2 <- ggplot(craydat, aes(x = Spread_dist, y = CL_ChelaL_nls_res, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm') + facet_wrap(~River_Tributary, scales='free_y')
chelaratio_spreadist_trib2
ggsave('cray_chelaratio2_spreadist_trib2.png', plot=chelaratio_spreadist_trib2, width=20, height=12, units='in')

chelaratio_meaninvyr <- ggplot(craydat, aes(x = meaninvyr, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))+
  geom_smooth(method='lm')
chelaratio_meaninvyr
ggsave('cray_chelaratio3_meaninvyr_trib.png', plot=chelaratio_meaninvyr, width=10, height=6, units='in')

chelaratio_meaninvyr_2 <- ggplot(craydat, aes(x = meaninvyr, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))+
  geom_smooth(method='lm') + 
  facet_wrap(~River_Tributary, scales='free_y')
chelaratio_meaninvyr_2
ggsave('cray_chelaratio3_meaninvyr_trib2.png', plot=chelaratio_meaninvyr_2, width=20, height=12, units='in')


chelaratio_interval_trib <- ggplot(craydat, aes(x =interval, y =  CL_ChelaL_nls_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) + 
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))
chelaratio_interval_trib
ggsave('cray_chelaratio5_interval_trib.png', plot=chelaratio_interval_trib, width=20, height=12, units='in')

#PRODUCTIVITY AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_ADFW <- ggplot(craydat, aes(x = AFDWmean, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  geom_smooth(method='lm')
chelaratio_ADFW 
ggsave('cray_chelaratio5_ADFW.png', plot=chelaratio_ADFW, width=20, height=12, units='in')

chelaratio_ADFW_trib <- ggplot(craydat, aes(x = AFDWmean, y = CL_ChelaL_nls_res, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  facet_wrap(~River_Tributary, scales='free_y')
chelaratio_ADFW_trib
ggsave('cray_chelaratio6_ADFW_trib.png', plot=chelaratio_ADFW_trib, width=20, height=12, units='in')

chelaratio_pp_trib <- ggplot(craydat, aes(x = greenmean+cyanomean+diatommean, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary, scales='free_y') +
  geom_smooth(method='lm')
chelaratio_pp_trib
ggsave('cray_chelaratio8_pp_trib.png', plot=chelaratio_pp_trib, width=20, height=12, units='in')

#COMPETITION AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_CPUE_trib <- ggplot(craydat, aes(x = Kick_mean, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')
chelaratio_CPUE_trib
ggsave('cray_chelaratio9_CPUE_trib.png', plot=chelaratio_CPUE_trib, width=20, height=12, units='in')

chelaratio_CPUE_trib2 <- ggplot(craydat, aes(x = Kick_mean, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')+
  facet_wrap(~River_Tributary, scales='free_y')
chelaratio_CPUE_trib2
ggsave('cray_chelaratio9_CPUE_trib2.png', plot=chelaratio_CPUE_trib2, width=20, height=12, units='in')

#HABITAT AND CHELA LENGTH/CARAPACE LENGTH RATIO
#Temperature
chelaratio_tempfield_trib <- ggplot(craydat, aes(x = tempmean, y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(method='lm')
chelaratio_tempfield_trib 
ggsave('cray_chelaratio90_tempfield_trib.png', plot=chelaratio_tempfield_trib)

#Discharge
chelaratio_discharge_trib <- ggplot(craydat, aes(x = (Q0001E_MA), y = CL_ChelaL_nls_res)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  facet_wrap(~River_Tributary, scale='free') +
  geom_smooth(span=1)
chelaratio_discharge_trib
ggsave('cray_chelaratio93_discharge_trib.png', plot=chelaratio_discharge_trib, width=20, height=12, units='in')

#Mesohabitat
chelaratio_Meso_type_trib <- ggplot(craydat, aes(x = (Meso_type), y = CL_ChelaL_nls_res)) +
  geom_boxplot(size=1, alpha=1/2)
chelaratio_Meso_type_trib
ggsave('cray_chelaratio930_Meso_type_trib.png', plot=chelaratio_Meso_type_trib)

#Method
chelaratio_method_trib <- ggplot(craydat, aes(x = Method, y = CL_ChelaL_nls_res)) +
  geom_boxplot(size=1, alpha=1/2) +
  geom_hline(yintercept=0, color='red')
chelaratio_method_trib
ggsave('cray_chelaratio97_method_trib.png', plot=chelaratio_method_trib)

#OTHER INTRINSIC PROPERTIES
# chelaratio_RNADNA_trib <- ggplot(craydat, aes(x = RNADNAratio, y = CL_ChelaL_nls_res)) +
#   geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
#   geom_text(aes(label = Site, size = 2), alpha=1/3)+ 
#   geom_smooth(span=1)
# chelaratio_RNADNA_trib
# ggsave('cray_chelaratio94_RNADNA_trib.png', plot=chelaratio_RNADNA_trib)

# chelaratio_RNADNA_site <- ggplot(craydat, aes(x = RNADNAratio, y = CL_ChelaL_nls_res)) +
#   geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
#   geom_text(aes(label = Site, size = 2), alpha=1/3)+ 
#   geom_smooth(method='lm') +
#   facet_wrap(~Site, scales='free')
# chelaratio_RNADNA_site
# ggsave('cray_chelaratio94_RNADNA_site.png', plot=chelaratio_RNADNA_site, width=20, height=12, units='in')

chelaratio_sex_trib <- ggplot(craydat, aes(x = Sex, y = CL_ChelaL_nls_res)) +
  geom_boxplot(size=1, alpha=1/2)
chelaratio_sex_trib
ggsave('cray_chelaratio95_sex_trib.png', plot=chelaratio_sex_trib)

chelaratio_weight_trib <- ggplot(craydat, aes(x = Weight, y = CL_ChelaL_nls_res)) +
  geom_point(size=1, alpha=1/2)
chelaratio_weight_trib

############## TROPHIC LEVEL ###################################################
#DISTANCE FROM INVASION AND TROPHIC LEVEL
TL_spreadist <- ggplot(craydat, aes(x = Spread_dist, y = trophiclevel, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL1_spreadist.png', plot=TL_spreadist)

TL_spreadist_trib <- ggplot(craydat, aes(x = Spread_dist, y = trophiclevel, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
ggsave('cray_TL2_spreadist_trib.png', plot=TL_spreadist_trib)

TL_spreadist_trib2 <- ggplot(craydat, aes(x = Spread_dist, y = trophiclevel, group=factor(Site))) +
  geom_violin(aes(color=factor(River_Tributary),width=10), size=1, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
TL_spreadist_trib2
ggsave('cray_TL2_spreadist_trib2.png', plot=TL_spreadist_trib2)

TL_meaninvyr <- ggplot(craydat, aes(x = meaninvyr, y = trophiclevel, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL3_meaninvyr.png', plot=TL_meaninvyr)

TL_meaninvyr_trib <-ggplot(craydat, aes(x = meaninvyr, y = trophiclevel)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
ggsave('cray_TL4_meaninvyr_trib.png', plot=TL_meaninvyr_trib)

TL_meaninvyr_trib2 <- ggplot(craydat, aes(x = meaninvyr, y = trophiclevel, group=factor(Site))) +
  geom_violin(aes(color=factor(River_Tributary),width=2), size=1, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(color=River_Tributary), method='glm')
TL_meaninvyr_trib2
ggsave('cray_TL4_meaninvyr_trib2.png', plot=TL_meaninvyr_trib2)

TL_interval_trib <- ggplot(craydat, aes(x =interval, y = trophiclevel, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) + 
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))
ggsave('cray_TL5_interval_trib.png', plot=TL_interval_trib)

#PRODUCTIVITY AND TROPHIC LEVEL
TL_ADFW_trib <- ggplot(craydat, aes(x = AFDWmean, y = trophiclevel, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('cray_TL7_ADFW_trib.png', plot=TL_ADFW_trib)

TL_ADFW_trib2 <- ggplot(craydat, aes(x = AFDWmean, y = trophiclevel)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + facet_wrap(~River_Tributary) + geom_smooth(method='lm')
TL_ADFW_trib2
ggsave('cray_TL7_ADFW_trib2.png', plot=TL_ADFW_trib2, width=10, height=6, units='in')


TL_pp_trib <- ggplot(craydat, aes(x = greenmean+cyanomean+diatommean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('cray_TL9_pp_trib.png', plot=TL_pp_trib)

TL_green_trib <- ggplot(craydat, aes(x = greenmean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
TL_green_trib
ggsave('cray_TL9_green_trib.png', plot=TL_green_trib)

TL_cyano_trib <- ggplot(craydat, aes(x = cyanomean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
TL_cyano_trib
ggsave('cray_TL9_cyano_trib.png', plot=TL_cyano_trib)

TL_diatom_trib <- ggplot(craydat, aes(x = diatommean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
TL_diatom_trib
ggsave('cray_TL9_diatom_trib.png', plot=TL_diatom_trib)

#COMPETITION AND TROPHIC LEVEL
TL_CPUE_trib <- ggplot(craydat, aes(x = Kick_mean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')
ggsave('cray_TL90_CPUE_trib.png', plot=TL_CPUE_trib)

#HABITAT AND TROPHIC LEVEL
#Temperature
TL_tempfield_trib <- ggplot(craydat, aes(x = tempmean, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL91_tempfield_trib.png', plot=TL_tempfield_trib)

TL_tempmod_trib <- ggplot(craydat, aes(x = degdays, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL92_tempmod_trib.png', plot=TL_tempmod_trib)

TL_tempmodAug_trib <- ggplot(craydat, aes(x = X2016.07.01, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL93_tempmodAug_trib.png', plot=TL_tempmodAug_trib)

#Discharge
TL_discharge_trib <- ggplot(craydat, aes(x = (Q0001E_MA), y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_TL94_discharge_trib.png', plot=TL_discharge_trib)

#OTHER INTRINSIC PROPERTIES
facetorder <- c(35,34,33,32,31,30,29,28,27,26,36,107,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,113,40,41,42,43,20,21,23,25)
craydat_stat$facet <- factor(craydat_stat$Site, levels = facetorder)
craydat_stat <- craydat_stat[order(craydat_stat$facet),]
craydat_stat$Site
craydat$facet <- factor(craydat$Site, levels = facetorder,
                        labels=as.character(round(craydat_stat[(craydat_stat$Site %in% facetorder), 'Spread_dist'],0)))
TL_CL_site <- ggplot(droplevels(craydat[(craydat$Site %in% craydat_stat[!is.na(craydat_stat$trophiclevel_mean),"Site"]) & craydat$Site != 107,]), aes(x = CL, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_smooth(method='lm') + 
  facet_grid(~facet,drop = TRUE) + 
  scale_x_continuous(name='Carapace length (mm)', expand=c(0,0),breaks=c(10,30,50)) +
  scale_y_continuous(name='Trophic position', expand=c(0,0)) +
  theme_bw()
TL_CL_site
ggsave('cray_TL95_CL_site.png', plot=TL_CL_site)


TL_CL_trib <- ggplot(craydat, aes(x = CL, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  geom_smooth(method='lm') 
TL_CL_trib
ggsave('cray_TL95_CL_trib.png', plot=TL_CL_trib)

TL_CL_trib2 <- ggplot(craydat, aes(x = CL, y = trophiclevel)) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2)) +
  geom_smooth(method='lm') + facet_wrap(~River_Tributary)
TL_CL_trib2
ggsave('cray_TL95_CL_trib2.png', plot=TL_CL_trib2, width=20, height=12, units='in')

TL_ChelaL_trib <- ggplot(craydat, aes(x = CL_ChelaL_nls_res, y = trophiclevel)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  #geom_text(aes(label = Site, size = 2)) +
  geom_smooth(method='lm')
TL_ChelaL_trib
ggsave('cray_TL95_ChelaL_trib.png', plot=TL_ChelaL_trib)

TL_ChelaL_trib2 <- ggplot(craydat, aes(x = CL_ChelaL_nls_res, y = trophiclevel)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2)) +
  geom_smooth(method='lm') +
  facet_wrap(~River_Tributary, scales='free')
TL_ChelaL_trib2
ggsave('cray_TL95_ChelaL_trib2.png', plot=TL_ChelaL_trib2)

TL_sex_trib <- ggplot(craydat, aes(x = Sex, y = trophiclevel)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('cray_TL96_sex_trib.png', plot=TL_sex_trib)

TL_sex <- ggplot(craydat, aes(x = Sex, y = trophiclevel)) +
  geom_boxplot(size=2, alpha=1/2)
ggsave('cray_TL97_sex.png', plot=TL_sex)

TL_molt <- ggplot(craydat, aes(x = Molting, y = trophiclevel)) +
  geom_boxplot(size=2, alpha=1/2)
ggsave('cray_TL98_molt.png', plot=TL_molt)

TL_molt_trib <- ggplot(craydat, aes(x = Molting, y = trophiclevel)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~River_Tributary)
ggsave('cray_TL99_molt_trib.png', plot=TL_molt_trib)

TL_molt_site <- ggplot(craydat, aes(x = Molting, y = trophiclevel)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~Site)
ggsave('cray_TL990_molt_site.png', plot=TL_molt_site)

TL_mesohab_trib <- ggplot(craydat, aes(x = Meso_type, y = trophiclevel)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~River_Tributary)
ggsave('cray_TL991_mesohab_trib.png', plot=TL_mesohab_trib)

TL_mesohab <- ggplot(craydat, aes(x = Meso_type, y = trophiclevel)) +
  geom_boxplot(alpha=1/2)
ggsave('cray_TL992_mesohab.png', plot=TL_mesohab)

TL_Miss_App <- ggplot(craydat, aes(x = Miss_App, y = trophiclevel)) +
  geom_boxplot(alpha=1/2)
ggsave('cray_TL993_Miss_App.png', plot=TL_Miss_App)

############## RNA/DNA ###################################################
#DISTANCE FROM INVASION AND RNA/DNA
RNADNAratio_spreadist <- ggplot(craydat, aes(x = Spread_dist, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_RNADNAratio1_spreadist.png', plot=RNADNAratio_spreadist)

RNADNAratio_spreadist_trib <- ggplot(craydat, aes(x = Spread_dist, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
ggsave('cray_RNADNAratio2_spreadist_trib.png', plot=RNADNAratio_spreadist_trib)

RNADNAratio_meaninvyr <- ggplot(craydat, aes(x = meaninvyr, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_RNADNAratio3_meaninvyr.png', plot=RNADNAratio_meaninvyr)

RNADNAratio_meaninvyr_trib <-ggplot(craydat, aes(x = meaninvyr, y = RNADNAratio)) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) +
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(aes(group=River_Tributary), method='glm')
ggsave('cray_RNADNAratio4_meaninvyr_trib.png', plot=RNADNAratio_meaninvyr_trib)

RNADNAratio_interval_trib <- ggplot(craydat, aes(x =interval, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2) + 
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))
ggsave('cray_RNADNAratio5_interval_trib.png', plot=RNADNAratio_interval_trib)

#PRODUCTIVITY AND RNA/DNA
RNADNAratio_ADFW <- ggplot(craydat, aes(x = AFDWmean, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(Site)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_RNADNAratio6_ADFW.png', plot=RNADNAratio_ADFW)

RNADNAratio_ADFW_trib <- ggplot(craydat, aes(x = AFDWmean, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('cray_RNADNAratio7_ADFW_trib.png', plot=RNADNAratio_ADFW_trib)

RNADNAratio_ADFW_trib2 <- ggplot(craydat, aes(x = AFDWmean, y = RNADNAratio, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() + 
  facet_wrap(~River_Tributary)
RNADNAratio_ADFW_trib2
ggsave('cray_RNADNAratio7_ADFW_trib2.png', plot=RNADNAratio_ADFW_trib2, width=10, height=6, units='in')

RNADNAratio_pp_trib <- ggplot(craydat, aes(x = greenmean+cyanomean+diatommean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('cray_RNADNAratio9_pp_trib.png', plot=RNADNAratio_pp_trib)

RNADNAratio_green_trib <- ggplot(craydat, aes(x = greenmean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
RNADNAratio_green_trib
ggsave('cray_RNADNAratio9_green_trib.png', plot=RNADNAratio_green_trib)

RNADNAratio_cyano_trib <- ggplot(craydat, aes(x = cyanomean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
RNADNAratio_cyano_trib
ggsave('cray_RNADNAratio9_cyano_trib.png', plot=RNADNAratio_cyano_trib)

RNADNAratio_diatom_trib <- ggplot(craydat, aes(x = diatommean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
RNADNAratio_diatom_trib
ggsave('cray_RNADNAratio9_diatom_trib.png', plot=RNADNAratio_diatom_trib)

#COMPETITION AND RNA/DNA
RNADNAratio_CPUE_trib <- ggplot(craydat, aes(x = Kick_mean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(aes(group=River_Tributary, color=River_Tributary), method='glm')
ggsave('cray_RNADNAratio90_CPUE_trib.png', plot=RNADNAratio_CPUE_trib)

#HABITAT AND RNA/DNA
#Temperature
RNADNAratio_tempfield_trib <- ggplot(craydat, aes(x = tempmean, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  facet_wrap(~River_Tributary)
RNADNAratio_tempfield_trib
ggsave('cray_RNADNAratio91_tempfield_trib.png', plot=RNADNAratio_tempfield_trib)

RNADNAratio_tempmod_trib <- ggplot(craydat, aes(x = degdays, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_RNADNAratio92_tempmod_trib.png', plot=RNADNAratio_tempmod_trib)

RNADNAratio_tempmodAug_trib <- ggplot(craydat, aes(x = X2016.07.01, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2)) +
  geom_smooth(aes(group=River_Tributary), method='lm')
RNADNAratio_tempmodAug_trib
ggsave('cray_RNADNAratio93_tempmodAug_trib.png', plot=RNADNAratio_tempmodAug_trib)

#Discharge
RNADNAratio_discharge_trib <- ggplot(craydat, aes(x = (Q0001E_MA), y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat[!duplicated(craydat$Site),], aes(label = Site, size = 2))
ggsave('cray_RNADNAratio94_discharge_trib.png', plot=RNADNAratio_discharge_trib)

#OTHER INTRINSIC PROPERTIES
RNADNAratio_trophiclevel_trib <- ggplot(craydat, aes(x = trophiclevel, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('cray_RNADNAratio95_trophiclevel_trib.png', plot=RNADNAratio_trophiclevel_trib)

RNADNAratio_trophiclevel_trib2 <- ggplot(craydat, aes(x = trophiclevel, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) +
  facet_wrap(~Site)
RNADNAratio_trophiclevel_trib2
ggsave('cray_RNADNAratio95_trophiclevel_trib2.png', plot=RNADNAratio_trophiclevel_trib2, width=20, height=12, units='in')

RNADNAratio_sex_trib <- ggplot(craydat, aes(x = Sex, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('cray_RNADNAratio96_sex_trib.png', plot=RNADNAratio_sex_trib)

RNADNAratio_sex <- ggplot(craydat, aes(x = Sex, y = RNADNAratio)) +
  geom_boxplot(size=2, alpha=1/2)
ggsave('cray_RNADNAratio96_sex.png', plot=RNADNAratio_sex)

RNADNAratio_sex_trib <- ggplot(craydat, aes(x = Sex, y = RNADNAratio)) +
  geom_boxplot(aes(color=River_Tributary), size=2, alpha=1/2)
ggsave('cray_RNADNAratio97_sex_trib.png', plot=RNADNAratio_sex_trib)

RNADNAratio_molt <- ggplot(craydat, aes(x = Molting, y = RNADNAratio)) +
  geom_boxplot(size=2, alpha=1/2)
ggsave('cray_RNADNAratio98_molt.png', plot=RNADNAratio_molt)

RNADNAratio_molt_trib <- ggplot(craydat, aes(x = Molting, y = RNADNAratio)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~River_Tributary)
ggsave('cray_RNADNAratio99_molt_trib.png', plot=RNADNAratio_molt_trib)

RNADNAratio_molt_site <- ggplot(craydat, aes(x = Molting, y = RNADNAratio)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~Site)
ggsave('cray_RNADNAratio990_molt_site.png', plot=RNADNAratio_molt_site)


RNADNAratio_mesohab_trib <- ggplot(craydat, aes(x = Meso_type, y = RNADNAratio)) +
  geom_boxplot(alpha=1/2) + facet_wrap(~River_Tributary)
ggsave('cray_RNADNAratio991_mesohab_trib.png', plot=RNADNAratio_mesohab_trib)

RNADNAratio_mesohab <- ggplot(craydat, aes(x = Meso_type, y = RNADNAratio)) +
  geom_boxplot(alpha=1/2)
ggsave('cray_RNADNAratio992_mesohab.png', plot=RNADNAratio_mesohab)

RNADNAratio_Miss_App_site <- ggplot(craydat, aes(x = Miss_App, y = RNADNAratio)) +
  geom_boxplot(alpha=1/2)
ggsave('cray_RNADNAratio993_Miss_App_site.png', plot=RNADNAratio_Miss_App_site)

RNADNAratio_CL_trib <- ggplot(craydat, aes(x = CL, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) +
  geom_smooth(method='lm')
RNADNAratio_CL_trib
ggsave('cray_RNADNAratio994_CL_trib.png', plot=RNADNAratio_CL_trib)

RNADNAratio_CL_trib2 <- ggplot(craydat, aes(x = CL, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) +
  geom_smooth(method='lm') + 
  facet_wrap(~Site, scales='free_y')
RNADNAratio_CL_trib2
ggsave('cray_RNADNAratio994_CL_trib2.png', plot=RNADNAratio_CL_trib2, width=20, height=12, units='in')

RNADNAratio_Weight_trib <- ggplot(craydat, aes(x = Weight_CL_lm_res, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) +
  geom_smooth(method='lm')
RNADNAratio_Weight_trib
ggsave('cray_RNADNAratio994_Weight_trib.png', plot=RNADNAratio_Weight_trib, width=20, height=12, units='in')

RNADNAratio_Weight_trib2 <- ggplot(craydat, aes(x = Weight_CL_lm_res, y = RNADNAratio)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) +
  geom_smooth(method='lm') + 
  facet_wrap(~Site, scales='free_y')
RNADNAratio_Weight_trib2
ggsave('cray_RNADNAratio994_Weight_trib2.png', plot=RNADNAratio_Weight_trib2, width=20, height=12, units='in')

############################################################################################
################################################# SITE-LEVEL DATA
############## GENERAL RELATIONSHIPS#############################
#Kickmean relationsip 
ggplot(craydat_stat, aes(x = Spread_dist, y = Kick_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat, aes(x = meaninvyr, y = Kick_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat, aes(x = tempmean, y = Kick_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat, aes(x = depthmean, y = Kick_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat, aes(x = velmean, y = Kick_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
#Molting
ggplot(craydat_stat[craydat_stat$ncray >5,], aes(x = degdays, y = molting, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat[craydat_stat$ncray >5,], aes(x = X2016.07.01, y = molting, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggplot(craydat_stat[craydat_stat$ncray >5,], aes(x = Spread_dist, y = molting, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()

############## SEX RATIO #######################
#DISTANCE FROM INVASION AND SEX RATIO
sexratio_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = Spread_dist, y = sexratio)) +
  geom_point(aes(color=factor(River_Tributary)), size=3, alpha=1/2)+ 
  geom_text(aes(label = ncray, size = 2), alpha=1/3, nudge_y = 0.01)  +
  geom_hline(yintercept=0.5) +
  geom_smooth(method='loess', aes(color=River_Tributary)) + 
  theme_bw()
sexratio_spreadist_trib
ggsave('site_sexratio1_spreadist_trib.png', plot=sexratio_spreadist_trib)

sexratio_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = Spread_dist, y = sexratio)) +
  geom_col(aes(fill=factor(River_Tributary)), size=1, width=3, alpha=1/2)+ 
  geom_text(aes(label = ncray, size = 2), alpha=1/3, nudge_y = 0.01)  +
  geom_hline(yintercept=0.5) +
  #geom_smooth(method='glm') + 
  theme_bw()
sexratio_spreadist_trib
ggsave('site_sexratio1_spreadist_trib.png', plot=sexratio_spreadist_trib)

sexratio_meaninvyr_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = meaninvyr, y = sexratio)) +
  geom_col(aes(fill=factor(River_Tributary)), size=1, width=0.25, alpha=1/2)+ 
  geom_text(aes(label = ncray, size = 2), alpha=1/3, nudge_y = 0.01)  +
  #geom_smooth(method='glm')+   
  geom_hline(yintercept=0.5) +
  theme_bw()
sexratio_meaninvyr_trib 
ggsave('site_sexratio2_meaninvyr_trib.png', plot=sexratio_meaninvyr_trib)

sexratio_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >5,], aes(x = sqrt(interval), y = sexratio)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm') + 
  theme_bw()
ggsave('site_sexratio3_interval_trib.png', plot=sexratio_interval_trib)

#COMPETITION AND SEX RATIO
sexratio_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >5,], aes(x = log(Kick_mean), y = sexratio)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3) +
  #geom_smooth(method='glm')+ 
  theme_bw()
ggsave('site_sexratio4_CPUE_trib.png', plot=sexratio_CPUE_trib)


############## MISSING APPENDAGES #######################
#DISTANCE FROM INVASION AND MISSING APPENDAGES
Miss_App_tot_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x = Spread_dist, y = Miss_App_tot)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)  +
  geom_smooth(method='glm') + 
  theme_bw()
ggsave('site_Miss_App_tot1_spreadist_trib.png', plot=Miss_App_tot_spreadist_trib)

Miss_App_tot_meaninvyr_trib <- ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x = meaninvyr, y = Miss_App_tot)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)+
  #geom_smooth(aes(group=River_Tributary), method='glm')+
  theme_bw()
ggsave('site_Miss_App_tot2_meaninvyr_trib.png', plot=Miss_App_tot_meaninvyr_trib)

Miss_App_tot_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x = sqrt(interval), y = Miss_App_tot)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3)+
  #geom_smooth(aes(group=River_Tributary), method='glm')+
  theme_bw()
ggsave('site_Miss_App_tot3_interval_trib.png', plot=Miss_App_tot_interval_trib)

#COMPETITION AND MISSING APPENDAGES
Miss_App_tot_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x = sqrt(Kick_mean), y = Miss_App_tot)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2), alpha=1/3) +
  geom_smooth(method='glm')+ 
  theme_bw()
ggsave('site_Miss_App_tot4_CPUE_trib.png', plot=Miss_App_tot_CPUE_trib)


############## CARAPACE LENGTH ################
#DISTANCE FROM INVASION AND LENGTH
  CL_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Spread_dist, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=CL_mean-CL_sd, ymax=CL_mean+CL_sd, color=River_Tributary), alpha=1/3, size=2) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw() + 
  geom_smooth(method='gam') +
  facet_wrap(~River_Tributary)
ggsave('site_CL1_spreadist_trib.png', plot=CL_spreadist_trib, width=10, height=6, units='in')
CL_spreadist_trib

CLsd_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Spread_dist, y = CL_sd)) +
  geom_point(aes(color=River_Tributary), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw() + 
  geom_smooth(method='glm', aes(group=River_Tributary))
ggsave('site_CLsd1_spreadist_trib.png', plot=CLsd_spreadist_trib, width=10, height=6, units='in')
CLsd_spreadist_trib

CL_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = meaninvyr, y = CL_mean)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=CL_mean-CL_sd, ymax=CL_mean+CL_sd), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3)+ 
  theme_bw()
ggsave('site_CL2_meaninvyr_trib.png', plot=CL_meaninvyr_trib)

CLsd_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = meaninvyr, y = CL_sd)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_smooth(alpha=1/3, span=1)+ 
  theme_bw()
ggsave('site_CLsd2_meaninvyr_trib.png', plot=CLsd_meaninvyr_trib)

CLmax_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = meaninvyr, y = CL_max)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_smooth(method='lm')+ 
  theme_bw()
CLmax_meaninvyr_trib

CLmax_spreaddist_trib <-ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = Spread_dist, y = CL_max)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_smooth(method='lm')+ 
  theme_bw()
CLmax_spreaddist_trib

CLqt90_spreaddist_trib <-ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = Spread_dist, y = CL_qt90)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_smooth(method='lm')+ 
  theme_bw()
CLqt90_spreaddist_trib

CLqt95_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = Spread_dist, y = CL_qt95)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_smooth(method='lm')+ 
  theme_bw()
CLqt95_meaninvyr_trib

CL_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x =interval, y = CL_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) + 
  geom_errorbar(aes(ymin=CL_mean-CL_sd, ymax=CL_mean+CL_sd), alpha=1/3) +
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  theme_bw()
ggsave('site_CL3_interval_trib.png', plot=CL_interval_trib)

#PRODUCTIVITY AND CARAPACE LENGTH
CL_AFDW_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = AFDWmean, y = CL_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
CL_AFDW_trib
ggsave('site_CL4_AFDW_trib.png', plot=CL_AFDW_trib)

CL_pp_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = greenmean+cyanomean+diatommean, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_CL5_pp_trib.png', plot=CL_pp_trib)

#COMPETITION AND CARAPACE LENGTH
CL_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Kick_mean, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(method='glm')
ggsave('site_CL6_CPUE_trib.png', plot=CL_CPUE_trib)

#HABITAT AND CARAPACE LENGTH
#Temperature
CL_tempfield_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = tempmean, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_CL7_tempfield_trib.png', plot=CL_tempfield_trib)

CL_tempmod_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = degdays, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_CL8_tempmod_trib.png', plot=CL_tempmod_trib)

#Discharge
CL_discharge_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = (Q0001E_MA), y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_CL9_discharge_trib.png', plot=CL_discharge_trib)

#OTHER INTRINSIC PROPERTIES
CL_RNADNA_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = RNADNAratio_mean, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3)
ggsave('site_CL90_RNADNA_trib.png', plot=CL_RNADNA_trib)

CL_trophiclevel_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = trophiclevel_mean, y = CL_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('site_CL91_trophiclevel_trib.png', plot=CL_trophiclevel_trib)

############## WEIGHT ############
#DISTANCE FROM INVASION AND WEIGHT/CARAPACE LENGTH
Weight_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = Spread_dist, y = Weight_res_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
  #geom_errorbar(aes(ymin=Weight_mean-Weight_SE, ymax=Weight_mean+Weight_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw()
ggsave('site_Weight1_spreadist_trib.png', plot=Weight_spreadist_trib)

Weight_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = meaninvyr, y = Weight_res_mean)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
  #geom_errorbar(aes(ymin=Weight_mean-Weight_SE, ymax=Weight_mean+Weight_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3)+ 
  theme_bw()
ggsave('site_Weight2_meaninvyr_trib.png', plot=Weight_meaninvyr_trib)

Weight_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x =interval, y = Weight_res_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) + 
  #geom_errorbar(aes(ymin=Weight_mean-Weight_SE, ymax=Weight_mean+Weight_SE), alpha=1/3) +
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  theme_bw()
ggsave('site_Weight3_interval_trib.png', plot=Weight_interval_trib)

#PRODUCTIVITY AND WEIGHT/CARAPACE LENGTH
Weight_AFDW_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = AFDWmean, y = Weight_res_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_Weight4_AFDW_trib.png', plot=Weight_AFDW_trib)

Weight_pp_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = greenmean+cyanomean+diatommean, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_Weight5_pp_trib.png', plot=Weight_pp_trib)

#COMPETITION AND WEIGHT/CARAPACE LENGTH
Weight_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = Kick_mean, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  #geom_smooth(method='glm')
ggsave('site_Weight6_CPUE_trib.png', plot=Weight_CPUE_trib)

#HABITAT AND WEIGHT/CARAPACE LENGTH
#Temperature
Weight_tempfield_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = tempmean, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
ggsave('site_Weight7_tempfield_trib.png', plot=Weight_tempfield_trib)

Weight_tempmod_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = degdays, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
ggsave('site_Weight8_tempmod_trib.png', plot=Weight_tempmod_trib)

#Discharge
Weight_discharge_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = (Q0001E_MA), y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
ggsave('site_Weight9_discharge_trib.png', plot=Weight_discharge_trib)

#OTHER INTRINSIC PROPERTIES
Weight_RNADNA_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = RNADNAratio_mean, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
ggsave('site_Weight90_RNADNA_trib.png', plot=Weight_RNADNA_trib)

Weight_trophiclevel_trib <- ggplot(craydat_stat[craydat_stat$ncray >10], aes(x = trophiclevel_mean, y = Weight_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site) & craydat_stat$ncray >10,], aes(label = Site, size = 2)) + 
ggsave('site_Weight91_trophiclevel_trib.png', plot=Weight_trophiclevel_trib)

############## CHELA LENGTH ################
#DISTANCE FROM INVASION AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x=Spread_dist, y =Chelae_res_mean)) +
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se, color=factor(River_Tributary)), alpha=1/3) +
  geom_point(aes(color=factor(River_Tributary)),size=6, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2)) + 
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw()
chelaratio_spreadist_trib 
ggsave('site_chelaratio1_spreadist_trib.png', plot=chelaratio_spreadist_trib)

chelaratio_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = meaninvyr, y =Chelae_res_mean)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3)+ 
  theme_bw()
ggsave('site_chelaratio2_meaninvyr_trib.png', plot=chelaratio_meaninvyr_trib)

chelaratio_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x =interval, y =Chelae_res_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) + 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  theme_bw()
ggsave('site_chelaratio3_interval_trib.png', plot=chelaratio_interval_trib)

#PRODUCTIVITY AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_AFDW_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = AFDWmean, y =Chelae_res_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_chelaratio4_AFDW_trib.png', plot=chelaratio_AFDW_trib)

chelaratio_pp_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10,], aes(x = greenmean+cyanomean+diatommean, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_chelaratio5_pp_trib.png', plot=chelaratio_pp_trib)

#COMPETITION AND CHELA LENGTH/CARAPACE LENGTH RATIO
chelaratio_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >10], aes(x = Kick_mean, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(method='glm')
ggsave('site_chelaratio6_CPUE_trib.png', plot=chelaratio_CPUE_trib)

#HABITAT AND CHELA LENGTH/CARAPACE LENGTH RATIO
#Temperature
chelaratio_tempfield_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10], aes(x = tempmean, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2))
ggsave('site_chelaratio7_tempfield_trib.png', plot=chelaratio_tempfield_trib)

chelaratio_tempmod_trib <- ggplot(craydat_stat[craydat_stat$ncray >=10], aes(x = degdays, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2))
ggsave('site_chelaratio8_tempmod_trib.png', plot=chelaratio_tempmod_trib)

#Discharge
chelaratio_discharge_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = (Q0001E_MA), y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2))
ggsave('site_chelaratio9_discharge_trib.png', plot=chelaratio_discharge_trib)

#OTHER INTRINSIC PROPERTIES
chelaratio_RNADNA_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = RNADNAratio_mean, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2), alpha=1/3)
ggsave('site_chelaratio90_RNADNA_trib.png', plot=chelaratio_RNADNA_trib)

chelaratio_trophiclevel_trib <- ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x = trophiclevel_mean, y =Chelae_res_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_errorbar(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=1/3) +
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('site_chelaratio91_trophiclevel_trib.png', plot=chelaratio_trophiclevel_trib)

############## TROPHIC LEVEL ##########
#DISTANCE FROM INVASION AND LENGTH
trophiclevel_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Spread_dist, y = trophiclevel_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=trophiclevel_mean-1.96*trophiclevel_SE, ymax=trophiclevel_mean+1.96*trophiclevel_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw()
ggsave('site_trophiclevel1_spreadist_trib.png', plot=trophiclevel_spreadist_trib)

trophiclevel_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = meaninvyr, y = trophiclevel_mean)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=trophiclevel_mean-1.96*trophiclevel_SE, ymax=trophiclevel_mean+1.96*trophiclevel_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3)+ 
  theme_bw()
ggsave('site_trophiclevel2_meaninvyr_trib.png', plot=trophiclevel_meaninvyr_trib)

trophiclevel_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x =interval, y = trophiclevel_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) + 
  geom_errorbar(aes(ymin=trophiclevel_mean-1.96*trophiclevel_SE, ymax=trophiclevel_mean+1.96*trophiclevel_SE), alpha=1/3) +
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  theme_bw()
ggsave('site_trophiclevel3_interval_trib.png', plot=trophiclevel_interval_trib)

#PRODUCTIVITY AND CARAPACE LENGTH
trophiclevel_AFDW_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = AFDWmean, y = trophiclevel_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_trophiclevel4_AFDW_trib.png', plot=trophiclevel_AFDW_trib)

trophiclevel_pp_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = greenmean+cyanomean+diatommean, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_trophiclevel5_pp_trib.png', plot=trophiclevel_pp_trib)

#COMPETITION AND CARAPACE LENGTH
trophiclevel_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Kick_mean, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(method='glm')
ggsave('site_trophiclevel6_CPUE_trib.png', plot=trophiclevel_CPUE_trib)

#HABITAT AND CARAPACE LENGTH
#Temperature
trophiclevel_tempfield_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = tempmean, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_trophiclevel7_tempfield_trib.png', plot=trophiclevel_tempfield_trib)

trophiclevel_tempmod_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = degdays, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_trophiclevel8_tempmod_trib.png', plot=trophiclevel_tempmod_trib)

trophiclevel_tempmod_mainstem <- ggplot(craydat_stat[craydat_stat$ncray >=10 & 
                                                   (craydat_stat$River_Tributary =='Upper mainstem' | craydat_stat$River_Tributary =='Lower mainstem'),],
                                    aes(x = degdays, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2))
ggsave('site_trophiclevel8_tempmod_mainstem.png', plot=trophiclevel_tempmod_mainstem)

#Discharge
trophiclevel_discharge_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = (Q0001E_MA), y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_trophiclevel9_discharge_trib.png', plot=trophiclevel_discharge_trib)

#OTHER INTRINSIC PROPERTIES
trophiclevel_RNADNA_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = RNADNAratio_mean, y = trophiclevel_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3)
ggsave('site_trophiclevel90_RNADNA_trib.png', plot=trophiclevel_RNADNA_trib)


############## RNA/DNA############
#DISTANCE FROM INVASION AND RNA/DNA
RNADNAratio_spreadist_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Spread_dist, y = RNADNAratio_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=RNADNAratio_mean-RNADNAratio_SE, ymax=RNADNAratio_mean+RNADNAratio_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3) + 
  theme_bw()
ggsave('site_RNADNAratio1_spreadist_trib.png', plot=RNADNAratio_spreadist_trib)

RNADNAratio_meaninvyr_trib <-ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = meaninvyr, y = RNADNAratio_mean)) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) +
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  geom_errorbar(aes(ymin=RNADNAratio_mean-RNADNAratio_SE, ymax=RNADNAratio_mean+RNADNAratio_SE), alpha=1/3) +
  #geom_smooth(aes(group=River_Tributary), method='glm', alpha=1/3)+ 
  theme_bw()
ggsave('site_RNADNAratio2_meaninvyr_trib.png', plot=RNADNAratio_meaninvyr_trib)

RNADNAratio_interval_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x =interval, y = RNADNAratio_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=6, alpha=1/2) + 
  geom_errorbar(aes(ymin=RNADNAratio_mean-RNADNAratio_SE, ymax=RNADNAratio_mean+RNADNAratio_SE), alpha=1/3) +
  scale_x_sqrt(breaks=c(0,1000,10000,25000,50000,100000,250000,400000))+ 
  theme_bw()
ggsave('site_RNADNAratio3_interval_trib.png', plot=RNADNAratio_interval_trib)

#PRODUCTIVITY AND RNA/DNA
RNADNAratio_AFDW_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = AFDWmean, y = RNADNAratio_mean, group=factor(Site))) +
  geom_point(aes(color=factor(River_Tributary)), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_RNADNAratio4_AFDW_trib.png', plot=RNADNAratio_AFDW_trib)

RNADNAratio_pp_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = greenmean+cyanomean+diatommean, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt()
ggsave('site_RNADNAratio5_pp_trib.png', plot=RNADNAratio_pp_trib)

#COMPETITION AND RNA/DNA
RNADNAratio_CPUE_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = Kick_mean, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2)) + 
  scale_x_sqrt() +
  geom_smooth(method='glm')
ggsave('site_RNADNAratio6_CPUE_trib.png', plot=RNADNAratio_CPUE_trib)

#HABITAT AND RNA/DNA
#Temperature
RNADNAratio_tempfield_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = tempmean, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_RNADNAratio7_tempfield_trib.png', plot=RNADNAratio_tempfield_trib)

RNADNAratio_tempmod_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = degdays, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_RNADNAratio8_tempmod_trib.png', plot=RNADNAratio_tempmod_trib)

#Discharge
RNADNAratio_discharge_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = (Q0001E_MA), y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(data=craydat_stat[!duplicated(craydat_stat$Site),], aes(label = Site, size = 2))
ggsave('site_RNADNAratio9_discharge_trib.png', plot=RNADNAratio_discharge_trib)

#OTHER INTRINSIC PROPERTIES
RNADNAratio_trophiclevel_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = trophiclevel_mean, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('site_RNADNAratio91_trophiclevel_trib.png', plot=RNADNAratio_trophiclevel_trib)

RNADNAratio_molting_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = molting, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('site_RNADNAratio92_molting_trib.png', plot=RNADNAratio_molting_trib)

RNADNAratio_sexratio_trib <- ggplot(craydat_stat[craydat_stat$ncray >=5,], aes(x = sexratio, y = RNADNAratio_mean)) +
  geom_point(aes(color=River_Tributary), size=4, alpha=1/2)+ 
  geom_text(aes(label = Site, size = 2), alpha=1/3) 
ggsave('site_RNADNAratio93_sexratio_trib.png', plot=RNADNAratio_sexratio_trib)


############## Look at kicknet CPUE across sites#############
ggplot(craydat_stat, aes(x = Spread_dist, y = Kick_mean, color = River_Tributary)) + geom_point(size = 2)

#Look at correlation between traps and kick-nets
ggplot(CPUEOR_meandistinfo, aes(x = Kick, y = Trap, color = River_Tributary)) + geom_point(size = 2) + scale_x_sqrt()
ggplot(CPUEOR_meandistinfo, aes(x = Density_estimate_min, y = Trap, color = River_Tributary)) + geom_point(size = 2) + 
  scale_x_continuous(limits=c(0,5))

##########################################################################################################################################################

################################################################################## 3. MODEL DATA AT CRAYFISH LEVEL #####################################################
setwd(file.path(resdir,"Datamod/"))
############################################################################################
################################################## LINEAR MIXED EFFECT MODELS ##############
#Good intro: http://www.ssc.wisc.edu/sscc/pubs/MM/MM_Models.html
#Also read lme4 book
#Might need to do GAMM to take in account sudden drops at last sites or make comparison between core and edge as categorical

######## CARAPACE LENGTH - LME - SIGNIFICANT HETEROSCEDASTICITY #########
#### Using samples from all methods but only from those sites with AFDW - with invasion year
CL_meaninvall_cray_lme <- lmer(log(CL) ~ sqrt(meaninvyr) + degdays + log(AFDWmean) + sqrt(Kick_mean) +
                                 sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                 Method + (1|River_Tributary), 
                        data=craydat[!is.na(craydat$AFDWmean),], REML=F)
summary(CL_meaninvall_cray_lme)
craydat[!is.na(craydat$AFDWmean),'CL_meaninvall_pred'] <- predict(CL_meaninvall_cray_lme, craydat[!is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$AFDWmean),], aes(x=CL_meaninvall_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
#The LRT of mixed models is only approximately ??2 distributed. 
#For tests of fixed effects the p-values will be smaller. 
#Thus if a p-value is greater than the cutoff value, you can be confident that a more accurate test would also retain the null hypothesis. 
#For p-values that are only a little below the cutoff value, a more accurate approach would need to be used.
drop1(CL_meaninvall_cray_lme,test="Chisq")
confint(CL_meaninvall_cray_lme)
plot(CL_meaninvall_cray_lme)
qqnorm(residuals(CL_meaninvall_cray_lme))

##### All samples (no AFDW) - with invasion year without degdays (too collinear and non-sensical relationship with CL)
CL_meaninvall_cray_lme <- lmer(log(CL) ~ meaninvyr + sqrt(Kick_mean) +
                                 sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                 Method + (1|River_Tributary), 
                               data=craydat[!is.na(craydat$Kick_mean),], REML=F)
summary(CL_meaninvall_cray_lme)
craydat[!is.na(craydat$Kick_mean),'CL_meaninvall_pred'] <- predict(CL_meaninvall_cray_lme, craydat[!is.na(craydat$Kick_mean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$Kick_mean),], aes(x=CL_meaninvall_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_meaninvall_cray_lme,test="Chisq")
confint(CL_meaninvall_cray_lme)
plot(CL_meaninvall_cray_lme)
qqnorm(residuals(CL_meaninvall_cray_lme))
qqline(residuals(CL_meaninvall_cray_lme))

##### All samples (no AFDW) - with invasion year without degdays or PP (too collinear and non-sensical relationship with CL)
ggplot(craydat[!is.na(craydat$Kick_mean),], aes(x=meaninvyr, y=Kick_mean)) + geom_point(aes(color=River_Tributary), size=4, alpha=1/2)
CL_meaninvall_cray_lme <- lmer(log(CL) ~ sqrt(meaninvyr) + sqrt(Kick_mean) +
                                 + sqrt(cyanomean) + sqrt(diatommean) +
                                 Method + (1|River_Tributary), 
                               data=craydat[!is.na(craydat$Kick_mean),], REML=F)
summary(CL_meaninvall_cray_lme)
craydat[!is.na(craydat$Kick_mean),'CL_meaninvall_pred'] <- predict(CL_meaninvall_cray_lme, craydat[!is.na(craydat$Kick_mean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$Kick_mean),], aes(x=CL_meaninvall_pred,y=CL, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_meaninvall_cray_lme,test="Chisq")
plot(CL_meaninvall_cray_lme)
qqnorm(residuals(CL_meaninvall_cray_lme))
qqline(residuals(CL_meaninvall_cray_lme))

##### Only samples with kick or snorkel/dipnet (no AFDW) - with invasion year without degdays (too collinear and non-sensical relationship with CL)
CL_meaninv_cray_lme <- lmer(log(CL) ~ sqrt(meaninvyr) + sqrt(Kick_mean) + Method + (1|River_Tributary), 
                                data=craydat[!is.na(craydat$Kick_mean),], REML=F)
summary(CL_meaninv_cray_lme)
craydat[,'CL_meaninv_pred'] <- predict(CL_meaninv_cray_lme,craydat,allow.new.levels=T)
ggplot(craydat, aes(x=CL_meaninv_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_meaninv_cray_lme,test="Chisq")
confint(CL_meaninv_cray_lme)
plot(CL_meaninv_cray_lme)
qqnorm(residuals(CL_meaninv_cray_lme))
qqline(residuals(CL_meaninv_cray_lme))


##### All samples (no AFDW) - without invasion year
CL_meaninvall_cray_lme <- lmer(log(CL) ~ sqrt(Kick_mean) + Method + (1|River_Tributary) + (1|Meso_type), 
                               data=craydat, REML=F)
summary(CL_meaninvall_cray_lme)
craydat[,'CL_meaninvall_pred'] <- predict(CL_meaninvall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=CL_meaninvall_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_meaninvall_cray_lme,test="Chisq")

##### All samples (no AFDW) - without Kick mean
CL_meaninvall_cray_lme <- lmer(log(CL) ~ sqrt(meaninvyr) + Method + (1|River_Tributary), 
                               data=craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',], REML=F)
summary(CL_meaninvall_cray_lme)
craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek','CL_meaninvall_pred'] <- predict(CL_meaninvall_cray_lme, craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',],allow.new.levels=T)
ggplot(craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',], aes(x=CL_meaninvall_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_meaninvall_cray_lme,test="Chisq")
resplot <- plot(CL_meaninvall_cray_lme)
png(filename = 'CL_LMEmeainvonly_residuals.png',
    width = 480, height = 480, units = "px")
resplot
dev.off()


##### All samples (no AFDW) - with spread dist
CL_spreadall_cray_lme <- lmer(log(CL) ~ log(Spread_dist) + degdays + sqrt(Kick_mean) +
                                 sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                 Method + (1|River_Tributary) + (1|Meso_type), 
                               data=craydat, REML=F)
summary(CL_spreadall_cray_lme)
craydat[,'CL_spreadall_pred'] <- predict(CL_spreadall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=CL_spreadall_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_spreadall_cray_lme,test="Chisq")

##### All samples (no AFDW) - without spread dist
CL_nospread_cray_lme <- lmer(log(CL) ~ degdays + sqrt(Kick_mean) +
                                sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                Method + (1|River_Tributary) + (1|Meso_type), 
                              data=craydat, REML=F)
summary(CL_nospread_cray_lme)
craydat[,'CL_nospread_pred'] <- predict(CL_nospread_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=CL_nospread_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_nospread_cray_lme,test="Chisq")

###### All samples (no AFDW) - no random effects
CL_spreadtrimfix_cray_lm <- lm(log(CL) ~ log(Spread_dist) + sqrt(Kick_mean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                    Method, data=craydat)
summary(CL_spreadtrimfix_cray_lm)
craydat[,'CL_spreadtrimfixlm_pred'] <- predict(CL_spreadtrimfix_cray_lm, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=CL_spreadtrimfixlm_pred,y=log(CL), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(CL_spreadtrimfix_cray_lm,test="Chisq")


############### CARAPACE LENGTH - QLM AND QLMM ############
CL_meaninvall_cray_rq <- rq(CL ~ meaninvyr, tau=0.95, data=craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',])
summary(CL_meaninvall_cray_rq, se="nid")
r1 <- resid(CL_meaninvall_cray_rq)
plot(r1)

#############################################With just invasion year and only kick sampling
par(mfrow=c(1,1))
attach(craydatOR[craydatOR$River_Tributary != 'Bridge Creek' & craydatOR$River_Tributary!='Beech Creek',])
plot(meaninvyr,CL,cex=.25,type="n",xlab="Invasion year", ylab="Carapace length")
points(meaninvyr,CL,cex=.5,col="blue")
abline(rq(CL~meaninvyr,tau=.5),col="blue")
abline(lm(CL~meaninvyr),lty=2,col="red") #the dreaded ols line
taus <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.97)
for( i in 1:length(taus)){
  abline(rq(CL~meaninvyr,tau=taus[i]),col="gray")
}

taus <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.97)
CLmeaninyr_lq<- rq(CL~meaninvyr,tau=taus)
Interc <- data.frame(taus)
Slo <- data.frame(taus)
#str(summary(CLmeaninyr_lq)[[1]])
for (x in 1:length(taus)) {
  Interc[x,2] <- summary(CLmeaninyr_lq, se="nid")[[x]]$coefficients[[1,1]]
  Interc[x,3] <- summary(CLmeaninyr_lq)[[x]]$coefficients[[1,2]]
  Interc[x,4] <- summary(CLmeaninyr_lq)[[x]]$coefficients[[1,3]]
  Slo[x,2] <- summary(CLmeaninyr_lq, se="nid")[[x]]$coefficients[[2,1]]
  Slo[x,3] <- summary(CLmeaninyr_lq)[[x]]$coefficients[[2,2]]
  Slo[x,4] <- summary(CLmeaninyr_lq)[[x]]$coefficients[[2,3]]
  Slo[x,5] <-  summary(CLmeaninyr_lq, se="nid")[[x]]$coefficients[[2,4]]
}

int_g <-ggplot(Interc, aes(x=taus, y=V2))+ geom_ribbon(aes(ymin=V3, ymax=V4), fill='orange', alpha=1/4) + geom_line(color='#0868ac', size=1.25) +
  theme_classic() + labs(x='Quantile', y='Intercept')
slo_g <-ggplot(Slo, aes(x=taus, y=V2)) + geom_abline(intercept=0,slope=0,size=1.5,alpha=1/2) + geom_line(color='#0868ac', size=1.25) + 
  geom_ribbon(aes(ymin=V3, ymax=V4), fill='orange', alpha=1/4)  + geom_point(color='red') +
  theme_classic() + labs(x='Quantile', y='Slope') + 
  scale_y_continuous(breaks=c(0.25,0,-0.25,-0.5)) +
  geom_point(data=Slo[Slo$V5 <= 0.10,], size=4, col='red')
grid.arrange(int_g, slo_g)
ggsave('CL_meaninvyr_LQR_2.png',grid.arrange(int_g, slo_g))


###################################### With Spread_dist
plot(CLmeaninyr_lq)
par(mfrow = c(1,1))
attach(craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',])
plot(Spread_dist,CL,cex=.25,type="n",xlab="Distance from introduction point", ylab="Carapace length")
points(Spread_dist,CL,cex=.5,col="blue")
abline(rq(CL~Spread_dist,tau=.5),col="blue")
abline(lm(CL~Spread_dist),lty=2,col="red") #the dreaded ols line
taus <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.97)
for( i in 1:length(taus)){
  abline(rq(CL~Spread_dist,tau=taus[i]),col="gray")
}

CL_spreadmethod <- rq(CL~Spread_dist+Method,tau=taus, data=craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',])
summary(CL_spreadmethod, se='nid')


#Only with kicknet data
par(mfrow = c(1,1))
attach(craydatOR[(craydatOR$River_Tributary == 'Lower mainstem' | craydatOR$River_Tributary == 'Upper mainstem' |
                   craydatOR$River_Tributary == 'North Fork' | craydatOR$River_Tributary == 'South Fork') &
                   craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10, 'Site'] &
                   craydatOR$Site != 107 & craydatOR$Site != 36,])
plot(Spread_dist,CL,cex=.25,type="n",xlab="Distance from introduction point", ylab="Carapace length")
points(Spread_dist,CL,cex=.5,col="blue")
abline(rq(CL~Spread_dist,tau=.5),col="blue")
abline(lm(CL~Spread_dist),lty=2,col="red") #the dreaded ols line
taus <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.97)
for( i in 1:length(taus)){
  abline(rq(CL~Spread_dist,tau=taus[i]),col="gray")
}

taus <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,0.95,0.97)
CLspread_lq<- rq(CL~Spread_dist,tau=taus, 
                 data=craydatOR[(craydatOR$River_Tributary == 'Lower mainstem' | craydatOR$River_Tributary == 'Upper mainstem' |
                                   craydatOR$River_Tributary == 'North Fork' | craydatOR$River_Tributary == 'South Fork') &
                                  (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10, 'Site']) &
                                  craydatOR$Site != 107 & craydatOR$Site != 36,])
Interc <- data.frame(taus)
Slo <- data.frame(taus)
#str(summary(CLmeaninyr_lq)[[1]])
for (x in 1:length(taus)) {
  Interc[x,2] <- summary(CLspread_lq, se="nid")[[x]]$coefficients[[1,1]]
  Interc[x,3] <- summary(CLspread_lq)[[x]]$coefficients[[1,2]]
  Interc[x,4] <- summary(CLspread_lq)[[x]]$coefficients[[1,3]]
  Slo[x,2] <- summary(CLspread_lq, se="nid")[[x]]$coefficients[[2,1]]
  Slo[x,3] <- summary(CLspread_lq)[[x]]$coefficients[[2,2]]
  Slo[x,4] <- summary(CLspread_lq)[[x]]$coefficients[[2,3]]
  Slo[x,5] <-  summary(CLspread_lq, se="nid")[[x]]$coefficients[[2,4]]
}

int_g <-ggplot(Interc, aes(x=taus, y=V2))+ geom_ribbon(aes(ymin=V3, ymax=V4), fill='orange', alpha=1/4) + geom_line(color='#0868ac', size=1.25) +
  theme_classic() + labs(x='Quantile', y='Intercept')
slo_g <-ggplot(Slo, aes(x=taus, y=V2)) + geom_abline(intercept=0,slope=0,size=1.5,alpha=1/2) + geom_line(color='#0868ac', size=1.25) + 
  geom_ribbon(aes(ymin=V3, ymax=V4), fill='orange', alpha=1/4)  + geom_point(color='red') +
  theme_classic() + labs(x='Quantile', y='Slope') + 
  scale_y_continuous(breaks=c(0.05,0.025, 0.01,0,-0.01,-0.025,-0.05)) +
  geom_point(data=Slo[Slo$V5 <= 0.05,], size=4, col='red')
grid.arrange(int_g, slo_g)
ggsave('CL_spreaddist_LQR_2.png',grid.arrange(int_g, slo_g))


##### TRY with random effect
craydat$Spread_dist_sc<- log(craydat$Spread_dist)
craydat$CL_sc <- log(craydat$CL)

str(craydat)
CL_spread_lqmm <- lqmm(fixed = CL~ Spread_dist, random = ~ 1, group = River_Tributary, tau = taus, nK = 9, type = "normal",
     data = craydat[craydat$River_Tributary != 'Bridge Creek' & craydat$River_Tributary!='Beech Creek',], 
     control = lqmmControl(method = "df", UP_max_iter = 50))
summary(CL_spread_lqmm)

CL_invyr_lqmm <- lqmm(fixed = CL~ meaninvyr, random = ~ 1, group =Method, tau = c(0.25,0.5,0.75,0.9,0.95,0.97), nK = 7, type = "normal",
                       data = craydat[craydat$River_Tributary =='Upper mainstem' | craydat$River_Tributary =='Lower mainstem' |
                                      craydat$River_Tributary =='North Fork'| craydat$River_Tributary == 'South Fork',], 
                       control = lqmmControl(method = "df", UP_max_iter = 200))
summary(CL_invyr_lqmm, R=50)
uhat.lqmm <- ranef(CL_invyr_lqmm)[["0.50"]]
uhat.lqmm
craydatsub <- craydat[craydat$River_Tributary =='Upper mainstem' | craydat$River_Tributary =='Lower mainstem' |
                        craydat$River_Tributary =='North Fork'| craydat$River_Tributary == 'South Fork',]
craydatsub$pred.lqmm <- predict(CL_invyr_lqmm, level = 1)[,4]

ggplot(craydatsub,
       aes(x=Spread_dist, y=CL, color=Method, shape=River_Tributary)) + 
  geom_point(size=4, alpha=1/2)+
  geom_point(aes(y=pred.lqmm), color='black')
  
  

############################################################################################
##### WEIGHT ###########
#Using samples from all methods but only from those sites with AFDW - with invasion year
Weight_meaninvall_cray_lme <- lmer(Weight_CL_lm_res ~sqrt(meaninvyr) + degdays + log(AFDWmean) + sqrt(Kick_mean) +
                                 sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) + 
                                 Method + (1|River_Tributary) + (1|Meso_type), 
                               data=craydat[!is.na(craydat$AFDWmean),], REML=F)
summary(Weight_meaninvall_cray_lme)
craydat[!is.na(craydat$AFDWmean),'Weight_meaninvall_pred'] <- predict(Weight_meaninvall_cray_lme, craydat[!is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$AFDWmean),], aes(x=Weight_meaninvall_pred,y=Weight_CL_lm_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(Weight_meaninvall_cray_lme,test="Chisq")
confint(Weight_meaninvall_cray_lme)
plot(Weight_meaninvall_cray_lme)
qqnorm(residuals(Weight_meaninvall_cray_lme))
qqline(residuals(Weight_meaninvall_cray_lme))


#Using samples from all methods for all sites with Kick CPUE
Weight_meaninvall_cray_lme <- lmer(Weight_CL_lm_res ~sqrt(meaninvyr) + (1|River_Tributary), 
                                   data=craydat, REML=F)
summary(Weight_meaninvall_cray_lme)
craydat[,'Weight_meaninvall_pred'] <- predict(Weight_meaninvall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=Weight_meaninvall_pred,y=Weight_CL_lm_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(Weight_meaninvall_cray_lme,test="Chisq")
confint(Weight_meaninvall_cray_lme)
plot(Weight_meaninvall_cray_lme)
qqnorm(residuals(Weight_meaninvall_cray_lme))
qqline(residuals(Weight_meaninvall_cray_lme))


#Using samples from all methods for all sites with Kick CPUE
qplot(sqrt(craydat$meaninvyr))
qplot(craydat$Weight_CL_lm_res)
Weight_meaninvall_cray_lm <- lm(Weight_CL_lm_res ~ sqrt(meaninvyr), data=craydat)
summary(Weight_meaninvall_cray_lm)
craydat[,'Weight_meaninvall_pred'] <- predict(Weight_meaninvall_cray_lm, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=Weight_meaninvall_pred,y=Weight_CL_lm_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(Weight_meaninvall_cray_lm,test="Chisq")

Weight_meaninvall_cray_lm <- lm(Weight_CL_lm_res ~ sqrt(meaninvyr)+sqrt(Kick_mean), data=craydat)
summary(Weight_meaninvall_cray_lm)
craydat[,'Weight_meaninvall_pred'] <- predict(Weight_meaninvall_cray_lm, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=Weight_meaninvall_pred,y=Weight_CL_lm_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(Weight_meaninvall_cray_lm,test="Chisq")

#Using spreading distance rather than invasion year
Weight_spreadall_cray_lm <- lm(Weight_CL_lm_res ~ log(Spread_dist) + sqrt(Kick_mean) +
                                     sqrt(diatommean),data=craydat[!is.na(craydat$Kick_mean),])
summary(Weight_spreadall_cray_lm)
craydat[!is.na(craydat$Kick_mean),'Weight_spreadall_pred'] <- predict(Weight_spreadall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=Weight_spreadall_pred,y=Weight_CL_lm_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(Weight_spreadall_cray_lme,test="Chisq")


##### CHELA LENGTH ###########
#Using samples from all methods but only from those sites with AFDW - with invasion year
chelaratio_meaninvall_cray_lme <- lmer(CL_ChelaL_nls_res ~ sqrt(meaninvyr) + sqrt(Kick_mean) + Sex +
                                     (1|River_Tributary), 
                                   data=craydat[!is.na(craydat$Chelae_L) & !is.na(craydat$Kick_mean)& !is.na(craydat$Sex),], REML=F)
summary(chelaratio_meaninvall_cray_lme)
craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean)& !is.na(craydat$Sex),'chelaratio_meaninvall_pred'] <- 
  predict(chelaratio_meaninvall_cray_lme, craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean)& !is.na(craydat$Sex),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean)& !is.na(craydat$Sex),], aes(x=chelaratio_meaninvall_pred,y=CL_ChelaL_nls_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)

drop1(chelaratio_meaninvall_cray_lme,test="Chisq")
confint(chelaratio_meaninvall_cray_lme )
plot(chelaratio_meaninvall_cray_lme )
qqnorm(residuals(chelaratio_meaninvall_cray_lme ))
qqline(residuals(chelaratio_meaninvall_cray_lme ))

chelaratio_meaninvall_cray_lme <- lmer(CL_ChelaL_nls_res ~ sqrt(meaninvyr) + sqrt(Kick_mean) + 
                                         (1|River_Tributary), 
                                       data=craydat[!is.na(craydat$Chelae_L) & !is.na(craydat$Kick_mean),], REML=F)
summary(chelaratio_meaninvall_cray_lme)
craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean),'chelaratio_meaninvall_pred'] <- predict(chelaratio_meaninvall_cray_lme, craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$Chelae_L)& !is.na(craydat$Kick_mean),], aes(x=chelaratio_meaninvall_pred,y=CL_ChelaL_nls_res, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)

drop1(chelaratio_meaninvall_cray_lme,test="Chisq")
confint(chelaratio_meaninvall_cray_lme )
plot(chelaratio_meaninvall_cray_lme )
qqnorm(residuals(chelaratio_meaninvall_cray_lme ))
qqline(residuals(chelaratio_meaninvall_cray_lme ))


chelaratio_meaninvall_cray_lm <- lm(CL_ChelaL_nls_res ~ sqrt(meaninvyr),
                                       data=craydat[!is.na(craydat$Chelae_L) & !is.na(craydat$Kick_mean),], REML=F)
summary(chelaratio_meaninvall_cray_lm)


##### TROPHIC LEVEL ############
qplot(craydat$X2016.07.01)
TL_meaninvall_cray_lme <- lmer(trophiclevel ~ sqrt(meaninvyr) + log(CL) + log(AFDWmean) + sqrt(Kick_mean) + 
                                 (1|River_Tributary) + (1|Meso_type), 
                               data=craydat[!is.na(craydat$AFDWmean),], REML=F)
summary(TL_meaninvall_cray_lme)
craydat[!is.na(craydat$AFDWmean),'TL_meaninvall_pred'] <- predict(TL_meaninvall_cray_lme, craydat[!is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$AFDWmean),], aes(x=TL_meaninvall_pred,y=trophiclevel, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(TL_meaninvall_cray_lme,test="Chisq")
confint(TL_meaninvall_cray_lme,oldNames = FALSE)
plot(TL_meaninvall_cray_lme)
qqnorm(residuals(TL_meaninvall_cray_lme))
qqline(residuals(TL_meaninvall_cray_lme))


TL_meaninvall_cray_lme <- lmer(trophiclevel ~ sqrt(meaninvyr) + log(CL) + log(AFDWmean) +
                                 (1+meaninvyr|River_Tributary) + (1|Meso_type), 
                               data=craydat[!is.na(craydat$AFDWmean),], REML=F)
summary(TL_meaninvall_cray_lme)
craydat[!is.na(craydat$AFDWmean),'TL_meaninvall_pred'] <- predict(TL_meaninvall_cray_lme, craydat[!is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$AFDWmean),], aes(x=TL_meaninvall_pred,y=trophiclevel, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(TL_meaninvall_cray_lme,test="Chisq")
confint(TL_meaninvall_cray_lme)

TL_meaninvall_cray_lme <- lmer(trophiclevel ~ log(CL) + sqrt(greenmean) +
                                 (1+meaninvyr|River_Tributary) + (1|Meso_type), 
                               data=craydat[!is.na(craydat$greenmean),], REML=F)
summary(TL_meaninvall_cray_lme)
craydat[,'TL_meaninvall_pred'] <- predict(TL_meaninvall_cray_lme, craydat[,],allow.new.levels=T)
ggplot(craydat, aes(x=TL_meaninvall_pred,y=trophiclevel, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(TL_meaninvall_cray_lme,test="Chisq")
confint(TL_meaninvall_cray_lme)


#Temperature not needed
TL_meaninvall_cray_lme <- lmer(trophiclevel ~ sqrt(meaninvyr) + log(CL) +
                                 (1+meaninvyr|River_Tributary) + (1|Meso_type), 
                               data=craydat, REML=F)
summary(TL_meaninvall_cray_lme)
craydat[,'TL_meaninvall_pred'] <- predict(TL_meaninvall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=TL_meaninvall_pred,y=trophiclevel, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(TL_meaninvall_cray_lme,test="Chisq")
confint(TL_meaninvall_cray_lme,oldNames = FALSE)


TL_meaninvall_cray_lme <- lmer(trophiclevel ~ sqrt(meaninvyr)
                                 (1+meaninvyr|River_Tributary), 
                               data=craydat, REML=F)
summary(TL_meaninvall_cray_lme)

#Plot random slopes and intercepts for meaninvyr (from https://www.r-bloggers.com/getting-the-most-of-mix-models-with-random-slopes/)
#extract fixed effects
a=fixef(TL_meaninvall_cray_lme)
#extract random effects
b=ranef(TL_meaninvall_cray_lme, condVar=TRUE)
# Extract the variances of the random effects
qq <- attr(b[[1]], "postVar")
e=(sqrt(qq)) 
e=e[2,2,] #here we want to access the Petal.Weigth, which is stored in column 2 in b[[1]], that's why I use the [,2,2]
#calculate CI's
liminf=(b[[1]][2]+a[2])-(e*2)
mean_=(b[[1]][2]+a[2])
limsup=(b[[1]][2]+a[2])+(e*2)
#Plot betas and its errors
dotchart(mean_$meaninvyr, labels = rownames(mean_), cex = 0.5,         xlim = c(-2,2), xlab = "betas")
#add CI's...
for (i in 1:nrow(mean_)){
  lines(x = c(liminf[i,1], limsup[i,1]), y = c(i,i)) 
}



TL_spreadall_cray_lme <- lmer(trophiclevel ~ log(Spread_dist) + log(CL) + 
                                (1|River_Tributary) + (1|Meso_type), 
                              data=craydat, REML=F)
summary(TL_spreadall_cray_lme)
craydat[,'TL_spreadall_pred'] <- predict(TL_spreadall_cray_lme, craydat,allow.new.levels=T)
ggplot(craydat, aes(x=TL_spreadall_pred,y=trophiclevel, color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(TL_spreadall_cray_lme,test="Chisq")
confint(TL_meaninvall_cray_lme,oldNames = FALSE)

##### RNA/DNA - POSITIVE BUT TEMP and ADFW COLLINEAR BUT NOT MAINLY DRIVEN BY TEMP###########
#Using samples from all methods but only from those sites with AFDW - with invasion year -- too much collinearity (e.g. invasion year negatively correlated)
RNADNAratio_meaninvall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(meaninvyr) + X2016.07.01 + log(CL) + sqrt(Kick_mean) +
                                          sqrt(greenmean) + sqrt(cyanomean) +
                                         log(AFDWmean)+ (1|River_Tributary) + (1|Meso_type), 
                                       data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$Kick_mean) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_meaninvall_cray_lme)
craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),'RNADNAratio_meaninvall_pred'] <- 
  predict(RNADNAratio_meaninvall_cray_lme, craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_meaninvall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_meaninvall_cray_lme,test="Chisq")
confint(RNADNAratio_meaninvall_cray_lme)
plot(RNADNAratio_meaninvall_cray_lme)
qqnorm(residuals(RNADNAratio_meaninvall_cray_lme))
qqline(residuals(RNADNAratio_meaninvall_cray_lme))


#Remove some collinearity but mean invasion year not significant anymore
RNADNAratio_meaninvall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(meaninvyr) + log(CL) + sqrt(greenmean) + 
                                          log(AFDWmean)+ (1|River_Tributary) + (1|Meso_type), 
                                        data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_meaninvall_cray_lme)
craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),'RNADNAratio_meaninvall_pred'] <- 
  predict(RNADNAratio_meaninvall_cray_lme, craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_meaninvall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_meaninvall_cray_lme,test="Chisq")

#Without green algae
RNADNAratio_meaninvall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(meaninvyr) + log(CL) + 
                                          log(AFDWmean)+ (1|River_Tributary) + (1|Meso_type), 
                                        data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_meaninvall_cray_lme)
craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),'RNADNAratio_meaninvall_pred'] <- 
  predict(RNADNAratio_meaninvall_cray_lme, craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_meaninvall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_meaninvall_cray_lme,test="Chisq")



#Just invasion year
RNADNAratio_meaninvall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(meaninvyr) + (1|River_Tributary) + (1|Meso_type), 
                                        data=craydat[!is.na(craydat$RNADNAratio),], REML=F)
summary(RNADNAratio_meaninvall_cray_lme)
craydat[!is.na(craydat$RNADNAratio),'RNADNAratio_meaninvall_pred'] <- 
  predict(RNADNAratio_meaninvall_cray_lme, craydat[!is.na(craydat$RNADNAratio),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio),], 
       aes(x=RNADNAratio_meaninvall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_meaninvall_cray_lme,test="Chisq")



RNADNAratio_spreadall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(Spread_dist) + degdays+ log(CL) + sqrt(Kick_mean) +
                                          sqrt(greenmean) + sqrt(cyanomean) + sqrt(diatommean) +
                                          log(AFDWmean)+ Method + (1|River_Tributary) + (1|Meso_type), 
                                        data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$Kick_mean) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_spreadall_cray_lme)
craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),'RNADNAratio_spreadall_pred'] <- 
  predict(RNADNAratio_spreadall_cray_lme, craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_spreadall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_spreadall_cray_lme,test="Chisq")


RNADNAratio_spreadall_cray_lme <- lmer(log(RNADNAratio) ~  sqrt(Spread_dist) + 
                                          log(AFDWmean)+ (1|River_Tributary), 
                                        data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_spreadall_cray_lme)
craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),'RNADNAratio_spreadall_pred'] <- 
  predict(RNADNAratio_spreadall_cray_lme, craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_spreadall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_spreadall_cray_lme,test="Chisq")

#Without Spread_dist but with deg_days
RNADNAratio_spreadall_cray_lme <- lmer(log(RNADNAratio) ~  degdays + 
                                         log(AFDWmean)+ (1|River_Tributary), 
                                       data=craydat[!is.na(craydat$RNADNAratio) & !is.na(craydat$AFDWmean),], REML=F)
summary(RNADNAratio_spreadall_cray_lme)
craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),'RNADNAratio_spreadall_pred'] <- 
  predict(RNADNAratio_spreadall_cray_lme, craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),],allow.new.levels=T)
ggplot(craydat[!is.na(craydat$RNADNAratio)& !is.na(craydat$Kick_mean)& !is.na(craydat$AFDWmean),], 
       aes(x=RNADNAratio_spreadall_pred,y=log(RNADNAratio), color=River_Tributary,shape=Method))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(RNADNAratio_spreadall_cray_lme,test="Chisq")



##### SEX RATIO #####
qplot(craydat_stat$sexratio)

sexratio_meaninvall_cray_lme <- lmer(sexratio ~  sqrt(meaninvyr) + sqrt(Kick_mean) + 
                                       (1|River_Tributary), 
                                        data=craydat_stat[craydat_stat$ncray >10,], REML=F)
summary(sexratio_meaninvall_cray_lme)
craydat_stat[craydat_stat$ncray >10,'sexratio_meaninvall_pred'] <- 
  predict(sexratio_meaninvall_cray_lme, craydat_stat[craydat_stat$ncray >10,],allow.new.levels=T)
ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x=sexratio_meaninvall_pred,y=sexratio, color=River_Tributary))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(sexratio_meaninvall_cray_lme,test="Chisq")

sexratio_spreadall_cray_lme <- lmer(sexratio ~  sqrt(Spread_dist) + sqrt(Kick_mean) + 
                                       (1|River_Tributary), 
                                     data=craydat_stat[craydat_stat$ncray >10,], REML=F)
summary(sexratio_spreadall_cray_lme)
craydat_stat[craydat_stat$ncray >10,'sexratio_spreadall_pred'] <- 
  predict(sexratio_spreadall_cray_lme, craydat_stat[craydat_stat$ncray >10,],allow.new.levels=T)
ggplot(craydat_stat[craydat_stat$ncray >10,], aes(x=sexratio_spreadall_pred,y=sexratio, color=River_Tributary))+
  geom_point()+geom_abline(slope=1, intercept=0)
drop1(sexratio_spreadall_cray_lme,test="Chisq")



##########################################################################################################################################################


################################################################################## 4. MODEL DATA AT SITE LEVEL #####################################################
#Discussion on testing for normality: https://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r/7788452#7788452
#For discussion of GAM see Chapter 3 and 4 in Wood (2006), for implementation of mgcv, see chapter 5 (p 217) 
#- see p222 and p224 for advantages and disadvantages of somoothing bases

#p239 for predictions

#Why use mgcv? By default, estimation of the degree of smoothness of model terms is part of model fitting. 
#Better confidence interval calculation, parametric part of the model can be penalized, 

#Why use thin plate regression spline (rather the B-spline, P-spline, or cubic splines)? (Chapter 4, Wood 2006)
#No need to define the position of the knots, cheap to compute, and defines the smoothing parameter)

craydat_statsub <- craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                  craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),]
# qplot(sqrt(craydat_statsub$meaninvyr))
# qplot(craydat_statsub$Chelae_res_mean)
# qplot(log10(craydat_statsub$Spread_dist))
# qplot(log10(craydat_statsub$CL_mean))
# qplot(craydat_statsub$Chelae_res_mean)
# qplot(craydat_statsub$Weight_res_mean)
#qplot(craydat_statsub$RNADNAratio_mean)
#qplot(log(craydat_statsub$RNADNAratio_mean))
# qplot(sqrt(craydat_statsub$greenmean))
# qplot(sqrt(craydat_statsub$cyanomean))
# qplot(sqrt(craydat_statsub$diatommean))
# qplot(log(craydat_statsub$AFDWmean))
# qplot(sqrt(craydat_statsub$Kick_mean))

################# Differences in environmental conditions for North Fork and South Fork ########################
#Sex ratio
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=Chelae_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=sexratio, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=sexratio, y=CL_sd)) + geom_col()

#CPUE
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=Chelae_res_mean)) + geom_col() 
ggplot(craydat_stat[craydat_stat$Site == 40 | craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=Chelae_res_mean)) + geom_col() + geom_text(aes(label=Site))

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=Kick_mean, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=Kick_mean, y=CL_sd)) + geom_col()

#Temperature
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=Chelae_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=degdays, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=degdays, y=CL_sd)) + geom_col()

#Green algae
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=Chelae_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=greenmean, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=greenmean, y=CL_sd)) + geom_col()


#Diatoms
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=Chelae_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=diatommean, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=diatommean, y=CL_sd)) + geom_col()


#AFDW macroinv. 
ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=Chelae_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=Chelae_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=trophiclevel_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=trophiclevel_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=RNADNAratio_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=RNADNAratio_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=CL_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=CL_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=Weight_res_mean)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=Weight_res_mean)) + geom_col()

ggplot(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], aes(x=AFDWmean, y=CL_sd)) + geom_col()
ggplot(craydat_stat[craydat_stat$Site == 41 | craydat_stat$Site == 43,], aes(x=AFDWmean, y=CL_sd)) + geom_col()


###################################### MODELS FOR CHELA LENGTH RESIDUALS ####################
##### All samples####
chela_sex_cray_lme <- lmer(CL_ChelaL_nls_res ~ Sex  + (1|Site), data=craydatOR, REML=F)
summary(CL_sex_cray_lme)
anova(CL_sex_cray_lme)

ggplot(craydatOR[craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site'],], aes(x=Sex, y=CL_ChelaL_nls_res)) + geom_boxplot() + facet_wrap(~Site)

chela_sex_ttest<- ddply(craydatOR[(craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site']),], 'Site', function(x) {
  print(c(unique(x$Site), unique(x$Spread_dist)))
  test <- t.test(x[x$Sex=='F','CL_ChelaL_nls_res'], x[x$Sex=='M','CL_ChelaL_nls_res'])
  data.frame(Site=x$Site[1], Spread_dist=x$Spread_dist[1], pvalue=test$p.value,estimate=test$estimate[2]-test$estimate[1])
})

x<-craydatOR[craydatOR$Site==2,]
test$estimate[2]

test <- t.test(craydatOR[craydatOR$Sex=='F','CL_ChelaL_nls_res'], craydatOR[craydatOR$Sex=='M','CL_ChelaL_nls_res'])
test

#####MAINSTEM######
#NULL MODEL
mgcv_chela_inter <- mgcv::gam(Chelae_res_mean~1, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_inter)
AIC(mgcv_chela_inter)

#With n=21
mgcv_chela_inter <- mgcv::gam(Chelae_res_mean~1, family=gaussian(link='identity'),
                              data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                  craydat_stat$ncray >10,])
summary(mgcv_chela_inter)
AIC(mgcv_chela_inter)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_chela_sp <- mgcv::gam(Chelae_res_mean~s(Spread_dist), family=gaussian(link='identity'),
                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                         craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_sp)
plot(mgcv_chela_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_sp)
par(mfrow=c(2,2))
gam.check(mgcv_chela_sp)
par(mfrow=c(1,1))


# gam_chela_sp <-  gam::gam(Chelae_res_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'),
#                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
#                                          craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_chela_sp)
# plot(gam_chela_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_chela_sp'] <- predict(mgcv_chela_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_chela_sp_se'] <- predict(mgcv_chela_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Chelae_res_mean))+
  geom_ribbon(aes(ymin = (gam_chela_sp-1.96*gam_chela_sp_se), ymax = (gam_chela_sp+1.96*gam_chela_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_chela_sp)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for invasion year (thin plate regression spline, penalized regression spline)
mgcv_chela_yr <- mgcv::gam(Chelae_res_mean~s(meaninvyr), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_yr)
plot(mgcv_chela_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_yr)
par(mfrow=c(2,2))
gam.check(mgcv_chela_yr)
par(mfrow=c(1,1))

#LOESS smoother for meaninvyr
# gam_chela_yr <-  gam(Chelae_res_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
#                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                          craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean),])
# summary(gam_chela_yr)
# plot(gam_chela_yr)

# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_chela_yr'] <- predict(gam_chela_yr)
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_chela_yr_se'] <- predict(gam_chela_yr, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 5 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Chelae_res_mean))+
  geom_ribbon(aes(ymin = (gam_chela_yr-1.96*gam_chela_yr_se), ymax = (gam_chela_yr+1.96*gam_chela_yr_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_chela_yr)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

##TS basis for environment (thin plate regression spline, penalized regression spline)
mgcv_chela_degdays <- mgcv::gam(Chelae_res_mean~s(degdays), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_degdays)
plot(mgcv_chela_degdays,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_degdays)
par(mfrow=c(2,2))
gam.check(mgcv_chela_degdays)
par(mfrow=c(1,1))

mgcv_chela_AFDW <- mgcv::gam(Chelae_res_mean~s(AFDWmean, k=3), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_AFDW)
plot(mgcv_chela_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_chela_AFDW)
par(mfrow=c(1,1))


# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
#                craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean), 'mgcv_chela_AFDW'] <- predict(mgcv_chela_AFDW)
# ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
#                       +                craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),],
#        aes(x=AFDWmean, y=Chelae_res_mean)) + geom_point() + geom_line(aes(y=mgcv_chela_AFDW)) + geom_text(aes(label=Site))
# 

mgcv_chela_AFDW <- mgcv::gam(Chelae_res_mean~s(AFDWmean, k=3), family=gaussian(link='identity'),
                                data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                    craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_AFDW)
plot(mgcv_chela_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_chela_AFDW)
par(mfrow=c(1,1))

mgcv_chela_green <- mgcv::gam(Chelae_res_mean~s(greenmean, k=3), family=gaussian(link='identity'),
                                data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                    craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_green)
plot(mgcv_chela_green,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_green)
par(mfrow=c(2,2))
gam.check(mgcv_chela_green)
par(mfrow=c(1,1))

mgcv_chela_diatom <- mgcv::gam(Chelae_res_mean~s(diatommean), family=gaussian(link='identity'),
                                data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                    craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_diatom )
plot(mgcv_chela_diatom,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_diatom) 
par(mfrow=c(2,2))
gam.check(mgcv_chela_diatom)
par(mfrow=c(1,1))

mgcv_chela_degAFDW <- mgcv::gam(Chelae_res_mean~s(degdays)+s(AFDWmean,k=3), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_degAFDW)
plot(mgcv_chela_degAFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_degAFDW)

#LOESS smooth for Environment
# gam_chela_env <-  gam(Chelae_res_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
#                       data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                           craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_chela_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_chela_env, resid=T, pch=16)

#Use anova to approximate significance of model terms and compute AIC


###WITH CPUE
mgcv_chela_CPUE <- mgcv::gam(Chelae_res_mean~s(Kick_mean,k=3), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10& !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_CPUE)
plot(mgcv_chela_CPUE,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_chela_CPUE)
par(mfrow=c(1,1))

#CPUE With all sites
mgcv_chela_CPUE <- mgcv::gam(Chelae_res_mean~s(Kick_mean,k=3), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10,])
summary(mgcv_chela_CPUE)
plot(mgcv_chela_CPUE,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_chela_CPUE)
par(mfrow=c(1,1))

#CPUE and Spread dist with all sites
mgcv_chela_CPUEsp <- mgcv::gam(Chelae_res_mean~s(Kick_mean,k=4) + s(Spread_dist), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_CPUEsp)
par(mfrow=c(1,2))
plot(mgcv_chela_CPUEsp,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_chela_CPUEsp)
par(mfrow=c(2,2))
gam.check(mgcv_chela_CPUEsp)
par(mfrow=c(1,1))
concurvity(mgcv_chela_CPUEsp)

###WITH SEX
mgcv_chela_sex <- mgcv::gam(Chelae_res_mean~s(sexratio), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_sex)
plot(mgcv_chela_sex,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_sex)
par(mfrow=c(2,2))
gam.check(mgcv_chela_sex)
par(mfrow=c(1,1))

#SEX WITH ALL SITES
mgcv_chela_sex <- mgcv::gam(Chelae_res_mean~s(sexratio), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10,])
summary(mgcv_chela_sex)
plot(mgcv_chela_sex,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_chela_sex)
par(mfrow=c(2,2))
gam.check(mgcv_chela_sex)
par(mfrow=c(1,1))

###WITH SEX AND SPREAD DIST
mgcv_chela_sexsp <- mgcv::gam(Chelae_res_mean~s(sexratio, k=4) + s(Spread_dist, k=3), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_chela_sexsp)
par(mfrow=c(1,2))
plot(mgcv_chela_sexsp,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_chela_sexsp)
par(mfrow=c(2,2))
gam.check(mgcv_chela_sexsp)
par(mfrow=c(1,1))
concurvity(mgcv_chela_sexsp)



#####NF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'North Fork' & craydat_stat$ncray >10,], aes(x=Spread_dist, y=Chelae_res_mean))+geom_col()+geom_errorbar(aes(ymin=Chelae_res_mean-1.96*(Chelae_res_mean/sqrt(ncray)), ymax=Chelae_res_mean+1.96*(Chelae_res_mean/sqrt(ncray))))
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' &  
                   craydatOR$Site %in%  craydat_stat[craydat_stat$ncray > 10,'Site'],], aes(x=CL_ChelaL_nls_res)) + geom_histogram()
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20,'CL_ChelaL_nls_res'])
shapiro.test(craydatOR[craydatOR$Site == 25,'CL_ChelaL_nls_res'])
#Test for equal variances
leveneTest(CL_ChelaL_nls_res ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
wilcox_chela_sp <- wilcox.test(craydatOR[craydatOR$Site==20,'CL'], craydatOR[craydatOR$Site==25,'CL_ChelaL_nls_res'],  var.equal = T)
wilcox_chela_sp


#####SF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'South Fork' & 
                      craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Chelae_res_mean))+geom_col()+geom_errorbar(aes(ymin=Chelae_res_mean-1.96*(CL_sd/sqrt(ncray)), ymax=Chelae_res_mean+1.96*(CL_sd/sqrt(ncray))))
ggplot(craydatOR[craydat_stat$River_Tributary == 'South Fork' & 
                   craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),], aes(x=CL_ChelaL_nls_res)) + geom_histogram()

ggplot(craydatOR[craydat_stat$Site ==41,], aes(x=CL_ChelaL_nls_res)) + geom_histogram()


#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'CL_ChelaL_nls_res'])
shapiro.test(craydatOR[craydatOR$Site == 43,'CL_ChelaL_nls_res'])
#Test for equal variances
leveneTest(CL_ChelaL_nls_res ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#T-test with welch's correction
wilcox_chela_sp <- wilcox.test(CL_ChelaL_nls_res ~ Site, 
                         data=craydatOR[craydatOR$River_Tributary == 'South Fork' & !is.na(craydatOR$AFDWmean) &  (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),],
                         var.equal=F)
wilcox_chela_sp

#####ALL SITES#####
gam_chela_sp <-  gam(Chelae_res_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
                     data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_chela_sp)
plot(gam_chela_sp)


gam_chela_yr <-  gam(Chelae_res_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                     data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_chela_yr)
plot(gam_chela_yr)


gam_chela_env <-  gam(Chelae_res_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                      data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_chela_env)

ggplot(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 5,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=2)
ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=1)

##############################################################################################

###################################### MODELS FOR TP #########################################
#Test the influence of Sex
TP_sex_cray_lme <- lmer(trophiclevel ~ Sex  + (1|Site), data=craydatOR, REML=F)
summary(TP_sex_cray_lme)
anova(TP_sex_cray_lme)

ggplot(craydatOR[craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site'],], aes(x=Sex, y=trophiclevel)) + geom_boxplot() + facet_wrap(~Site)

TP_sex_ttest<- ddply(craydatOR[!is.na(craydatOR$trophiclevel),], 'Site', function(x) {
  #print(c(unique(x$Site), unique(x$Spread_dist)))
  test <- tryCatch(t.test(x[x$Sex=='F','trophiclevel'], x[x$Sex=='M','trophiclevel']), error=function(e) NA)
  if (!is.na(test)) {
    data.frame(Site=x$Site[1], Spread_dist=x$Spread_dist[1], pvalue=test$p.value,estimate=test$estimate[2]-test$estimate[1])
  }
})

ggplot(craydatOR[craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site'],], aes(x=Sex, y=trophiclevel)) + geom_boxplot()
test <- t.test(craydatOR[craydatOR$Sex=='F','trophiclevel'], craydatOR[craydatOR$Sex=='M','trophiclevel'], var.equal=T)
test

#Test the influence of carapace length
TP_inter_cray_lme <- lmer(trophiclevel ~ 1  + (1|Site), data=craydatOR, REML=F)
summary(TP_inter_cray_lme)

TP_CL_cray_lme <- lmer(trophiclevel ~ CL  + (1|Site), data=craydatOR, REML=F)
summary(TP_CL_cray_lme)
anova(TP_CL_cray_lme)

#####MAINSTEM#####
#NULL MODEL
mgcv_TP_inter <- mgcv::gam(trophiclevel_mean~1, family=gaussian(link='identity'),
                              data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                  craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_inter)
AIC(mgcv_TP_inter)

qplot((craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 10,'trophiclevel_mean']))

#Thin plate regression smooter for spread_dist (basis dimension was determined iteratively)
mgcv_TP_sp <- gam(trophiclevel_mean~s(Spread_dist, bs='tp',k=4), family=gaussian(link='identity'),
                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                         craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_sp)
plot(mgcv_TP_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_sp)
par(mfrow=c(2,2))
gam.check(mgcv_TP_sp)
par(mfrow=c(1,1))

#Thin plate regression smooter for CL (basis dimension was determined iteratively)
mgcv_TP_CL <- gam(trophiclevel_mean~s(CL_mean, k=3), family=gaussian(link='identity'),
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_CL)
plot(mgcv_TP_CL,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_CL)
par(mfrow=c(2,2))
gam.check(mgcv_TP_CL)
par(mfrow=c(1,1))

mgcv_TP_spCL <- gam(trophiclevel_mean~s(Spread_dist, bs='tp',k=4) + s(CL_mean, k=3), family=gaussian(link='identity'),
                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_spCL)
par(mfrow=c(1,2))
plot(mgcv_TP_spCL,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_sp)
par(mfrow=c(2,2))
gam.check(mgcv_TP_sp)
par(mfrow=c(1,1))


#LOESS smoother for Spread_dist
# gam_TP_sp <-  gam(trophiclevel_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
#                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                         craydat_stat$ncray >10 & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_TP_sp)
#plot(gam_TP_sp)


craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_TP_sp'] <- predict(mgcv_TP_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_TP_sp_se'] <- predict(mgcv_TP_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=trophiclevel_mean))+
  geom_ribbon(aes(ymin = (gam_TP_sp-1.96*gam_TP_sp_se), ymax = (gam_TP_sp+1.96*gam_TP_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=gam_TP_sp), size=1.1) +
  scale_y_continuous(name='Carapace length standard deviation (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 



#Thin plate regression smooter for mean invasion year (basis dimension was determined iteratively)
mgcv_TP_yr <- gam(trophiclevel_mean~s(meaninvyr, bs='tp',k=4), family=gaussian(link='identity'),
                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_yr)
plot(mgcv_TP_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_yr)


#Thin plate regression smooter for invasion year (basis dimension was determined iteratively)
mgcv_TP_yr <- gam(trophiclevel_mean~s(meaninvyr, bs='tp',k=4), family=gaussian(link='identity'),
                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_yr)
plot(mgcv_TP_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_yr)
par(mfrow=c(2,2))
gam.check(mgcv_TP_yr)
par(mfrow=c(1,1))

#LOESS smoother for meaninvyr
# gam_TP_yr <-  gam(trophiclevel_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
#                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                         craydat_stat$ncray >10   & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_TP_yr) #Not significant

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean), 'mgcv_TP_yr'] <- predict(mgcv_TP_yr)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean), 'mgcv_TP_yr_se'] <- predict(mgcv_TP_yr, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10  & !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),], aes(x=meaninvyr, y=trophiclevel_mean))+
  geom_ribbon(aes(ymin = (gam_TP_yr-1.96*gam_TP_yr_se), ymax = (gam_TP_yr+1.96*gam_TP_yr_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=gam_TP_yr), size=1.1) +
  scale_y_continuous(name='Carapace length standard deviation (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 


#Thin plate regression smoother for environment
mgcv_TP_degdays <- gam(trophiclevel_mean~s(degdays, k=4), family=gaussian(link='identity'),
                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_degdays)
plot(mgcv_TP_degdays,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_degdays)
par(mfrow=c(2,2))
gam.check(mgcv_TP_degdays)
par(mfrow=c(1,1))


mgcv_TP_AFDW <- gam(trophiclevel_mean~s(AFDWmean), family=gaussian(link='identity'),
                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                       craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_AFDW)
plot(mgcv_TP_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_TP_AFDW)
par(mfrow=c(1,1))

mgcv_TP_green <- gam(trophiclevel_mean~s(greenmean), family=gaussian(link='identity'),
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_green)
plot(mgcv_TP_green,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_green)
par(mfrow=c(2,2))
gam.check(mgcv_TP_green)
par(mfrow=c(1,1))

mgcv_TP_diatom<- gam(trophiclevel_mean~s(diatommean,k=3), family=gaussian(link='identity'),
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_diatom)
plot(mgcv_TP_diatom,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_diatom)
par(mfrow=c(2,2))
gam.check(mgcv_TP_diatom)
par(mfrow=c(1,1))


#Add null space penalty such that the smoothing parameter selection can remove terms from the model altogether
qplot(craydat_statsub$degdays, craydat_statsub$diatommean)
qplot(craydat_statsub$degdays, craydat_statsub$AFDWmean)
qplot(craydat_statsub$AFDWmean, craydat_statsub$diatommean)

mgcv_TP_degAFDW <- gam(trophiclevel_mean~s(degdays, k=4)+s(AFDWmean, k=3), family=gaussian(link='identity'), select=T,
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_degAFDW)
mgcv_TP_degAFDW
par(mfrow=c(1,2))
plot(mgcv_TP_degAFDW,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_degAFDW)
par(mfrow=c(2,2))
gam.check(mgcv_TP_degAFDW)
par(mfrow=c(1,1))
concurvity(mgcv_TP_degAFDW)

mgcv_TP_degdia <- gam(trophiclevel_mean~s(degdays, k=4)+s(diatommean, k=3), family=gaussian(link='identity'), select=T,
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_degdia)
par(mfrow=c(1,2))
plot(mgcv_TP_degdia,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_degdia)
par(mfrow=c(2,2))
gam.check(mgcv_TP_degdia)
par(mfrow=c(1,1))
#0 indicates no problem, and 1 indicates total lack of identifiability
#The three indices are all based on the idea that a smooth term, f, in the model can be decomposed into a part, 
#g, that lies entirely in the space of one or more other terms in the model,
#and a remainder part that is completely within the term's own space. 
#If g makes up a large part of f then there is a concurvity problem.
concurvity(mgcv_TP_degdia)
#Diatom and degree days highly collinear

mgcv_TP_AFDWdia <- gam(trophiclevel_mean~s(AFDWmean, k=3)+s(diatommean, k=3), family=gaussian(link='identity'), select=T,
                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                          craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_AFDWdia)
par(mfrow=c(1,2))
plot(mgcv_TP_AFDWdia,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_AFDWdia)
par(mfrow=c(2,2))
gam.check(mgcv_TP_AFDWdia)
par(mfrow=c(1,1))
concurvity(mgcv_TP_AFDWdia)

#LOESS smoother for Environment
# gam_TP_env <-  gam(trophiclevel_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
#                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                          craydat_stat$ncray >10 &  !is.na(craydat_stat$trophiclevel_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_TP_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_TP_env, resid=T, pch=16)

#Check relationship with sex ratio
mgcv_TP_sex <- gam(trophiclevel_mean~s(sexratio,k=4), family=gaussian(link='identity'),
                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                       craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_sex)
plot(mgcv_TP_sex,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_sex)
par(mfrow=c(2,2))
gam.check(mgcv_TP_sex)
par(mfrow=c(1,1))


mgcv_TP_spsex <- gam(trophiclevel_mean~s(Spread_dist,k=5) + s(sexratio, k=5), family=gaussian(link='identity'),
                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                         craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_spsex)
par(mfrow=c(1,2))
plot(mgcv_TP_spsex,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_spsex)
par(mfrow=c(2,2))
gam.check(mgcv_TP_spsex)
par(mfrow=c(1,1))
concurvity(mgcv_TP_spsex)


#Check relationship with CPUE
mgcv_TP_CPUE <- gam(trophiclevel_mean~s(Kick_mean), family=gaussian(link='identity'),
                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                       craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_CPUE)
plot(mgcv_TP_CPUE,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_TP_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_TP_CPUE)
par(mfrow=c(1,1))

mgcv_TP_CPUEsp <- gam(trophiclevel_mean~s(Spread_dist,k=3)+s(Kick_mean,k=3), family=gaussian(link='identity'),
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_CPUEsp)
par(mfrow=c(1,2))
plot(mgcv_TP_CPUEsp,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_CPUEsp)
par(mfrow=c(2,2))
gam.check(mgcv_TP_CPUEsp)
par(mfrow=c(1,1))
concurvity(mgcv_TP_CPUEsp)

mgcv_TP_CPUEdeg <- gam(trophiclevel_mean~s(degdays,k=3)+s(Kick_mean,k=4), family=gaussian(link='identity'),
                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                          craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_TP_CPUEdeg)
par(mfrow=c(1,2))
plot(mgcv_TP_CPUEdeg,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_TP_CPUEdeg)
par(mfrow=c(2,2))
gam.check(mgcv_TP_CPUEdeg)
par(mfrow=c(1,1))
concurvity(mgcv_TP_CPUEdeg)

#####NORTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' & !is.na(craydatOR$AFDWmean) & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=trophiclevel)) + geom_bar() + facet_wrap(~Site)
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20,'trophiclevel'])
shapiro.test(craydatOR[craydatOR$Site == 25,'trophiclevel'])
#Test for equal variances
leveneTest(trophiclevel ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
ttest_TP_sp <- t.test(craydatOR[craydatOR$Site==20,'trophiclevel'], craydatOR[craydatOR$Site==25,'trophiclevel'],var.equal = T)
ttest_TP_sp

#####SOUTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=trophiclevel)) + geom_bar() + facet_wrap(~Site)
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'trophiclevel'])
shapiro.test(craydatOR[craydatOR$Site == 43,'trophiclevel'])
#Test for equal variances
leveneTest(trophiclevel ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
wilcox_TP_sp <- wilcox.test(craydatOR[craydatOR$Site==41,'trophiclevel'], craydatOR[craydatOR$Site==43,'trophiclevel'])
wilcox_TP_sp


#####ALL SITES#####
gam_TP_sp <-  gam(trophiclevel_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                    data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$trophiclevel_mean) & !is.na(craydat_stat$AFDWmean),])
summary(gam_TP_sp)


gam_TP_sp <-  gam(trophiclevel_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                     data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$trophiclevel_mean) & !is.na(craydat_stat$AFDWmean),])
summary(gam_TP_env)

ggplot(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 5,], aes(x=meaninvyr, y=trophiclevel_mean))+geom_point()+geom_smooth(span=2)
ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x=meaninvyr, y=trophiclevel_mean))+geom_point()+geom_smooth(span=1)


#############################################################################################

###################################### MODELS FOR RNA/DNA ###################################
#####MAINSTEM#####
qplot((craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 10,'RNADNAratio_mean']))
qplot(log(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 10,'RNADNAratio_mean']))

#NULL MODEL
mgcv_RD_inter <- mgcv::gam(RNADNAratio_mean~1, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_inter)
AIC(mgcv_RD_inter)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_RD_sp <- mgcv::gam(RNADNAratio_mean~Spread_dist, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_sp)
confint(mgcv_RD_sp)
plot(mgcv_RD_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_sp)
par(mfrow=c(2,2))
gam.check(mgcv_RD_sp)
par(mfrow=c(1,1))
0.006335-(1.96*0.004612)
0.006335+(1.96*0.004612)

mgcv_RD_sp <- glm(RNADNAratio_mean~(Spread_dist), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_sp)
confint(mgcv_RD_sp)

qplot(x=c(1,1,1), y=c(-0.002703, 0.006, 0.015)) + geom_hline(yintercept=0)

#LOESS smoother
# gam_RD_sp <-  gam(RNADNAratio_mean~lo(Spread_dist, span=1), family=gaussian(link='log'),
#                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
#                                       craydat_stat$ncray >10 & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_RD_sp)
# plot(gam_RD_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_RD_sp'] <- predict(mgcv_RD_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_RD_sp_se'] <- predict(mgcv_RD_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=RNADNAratio_mean))+
  geom_ribbon(aes(ymin = (gam_RD_sp-1.96*gam_RD_sp_se), ymax = (gam_RD_sp+1.96*gam_RD_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_RD_sp)), size=1.1) +
  scale_y_continuous(name='Carapace length standard deviation (mm)') +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for invasion year (thin plate regression spline, penalized regression spline)
mgcv_RD_yr <- mgcv::gam(RNADNAratio_mean~s(meaninvyr), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_yr)
plot(mgcv_RD_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_yr)
par(mfrow=c(2,2))
gam.check(mgcv_RD_yr)
par(mfrow=c(1,1))

#LOESS smoother for invasion year
# gam_RD_yr <-  gam(RNADNAratio_mean~lo(meaninvyr, span=1), family=gaussian(link='log'), 
#                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                       craydat_stat$ncray >10 & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_RD_yr)
# #plot(gam_RD_sp)
# 
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_RD_yr'] <- predict(gam_RD_yr)
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean), 'gam_RD_yr_se'] <- predict(gam_RD_yr, se.fit=T)$se.fit
# 
# ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                       craydat_stat$ncray >10  & !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean),], aes(x=meaninvyr, y=RNADNAratio_mean))+
#   geom_ribbon(aes(ymin = exp(gam_RD_yr-1.96*gam_RD_yr_se), ymax = exp(gam_RD_yr+1.96*gam_RD_yr_se), fill=River_Tributary), alpha=0.3)+
#   geom_point(aes(color=River_Tributary), size=4)+ 
#   geom_line(aes(y=exp(gam_RD_yr)), size=1.1) +
#   scale_y_continuous(name='Carapace length standard deviation (mm)', breaks=c(10,15,20,25,30)) +
#   scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
#   geom_text(aes(label=Site)) +
#   scale_fill_discrete(name='Tributary (95% CI)') + 
#   scale_color_discrete(name='Tributary (95% CI)') +
#   theme_classic() 

#TS basis for environment
mgcv_RD_degdays <- mgcv::gam(RNADNAratio_mean~s(degdays), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_degdays)
plot(mgcv_RD_degdays,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_degdays)
par(mfrow=c(2,2))
gam.check(mgcv_RD_degdays)
par(mfrow=c(1,1))

mgcv_RD_AFDW <- mgcv::gam(RNADNAratio_mean~s(AFDWmean), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_AFDW)
plot(mgcv_RD_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_RD_AFDW)
par(mfrow=c(1,1))

mgcv_RD_green <- mgcv::gam(RNADNAratio_mean~s(greenmean, k=2), family=gaussian(link='identity'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_green)
plot(mgcv_RD_green,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_green)
par(mfrow=c(2,2))
gam.check(mgcv_RD_green)
par(mfrow=c(1,1))

mgcv_RD_diatom <- mgcv::gam(RNADNAratio_mean~s(diatommean), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_diatom)
plot(mgcv_RD_diatom,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_diatom)
par(mfrow=c(2,2))
gam.check(mgcv_RD_diatom)
par(mfrow=c(1,1))

mgcv_RD_deggreen <- mgcv::gam(RNADNAratio_mean~s(greenmean, k=3)+s(degdays), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_deggreen)
par(mfrow=c(1,2))
plot(mgcv_RD_deggreen,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_RD_deggreen)
par(mfrow=c(2,2))
gam.check(mgcv_RD_deggreen)
par(mfrow=c(1,1))
concurvity(mgcv_RD_deggreen)

mgcv_RD_green1 <- mgcv::gam(RNADNAratio_mean~s(AFDWmean)+s(greenmean, k=3), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_green1)
par(mfrow=c(1,2))
plot(mgcv_RD_green1,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_RD_green1)
par(mfrow=c(2,2))
gam.check(mgcv_RD_green1)
par(mfrow=c(1,1))
concurvity(mgcv_RD_green1)

mgcv_RD_green2 <- mgcv::gam(RNADNAratio_mean~s(greenmean, k=3)+AFDWmean, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_green2)
par(mfrow=c(1,2))
plot(mgcv_RD_green2,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_RD_green2)
par(mfrow=c(2,2))
gam.check(mgcv_RD_green2)
par(mfrow=c(1,1))
concurvity(mgcv_RD_green2)

#Remove outlier site 15
mgcv_RD_green3 <- mgcv::gam(RNADNAratio_mean~s(AFDWmean)+s(greenmean, k=3), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean) & craydat_stat$Site != 15,])
summary(mgcv_RD_green3)
par(mfrow=c(1,2))
plot(mgcv_RD_green3,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_RD_green3)


#Environment
# gam_RD_env <-  gam(RNADNAratio_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='log'), 
#                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                        craydat_stat$ncray >10 &  !is.na(craydat_stat$RNADNAratio_mean) &  !is.na(craydat_stat$AFDWmean),])
# summary(gam_RD_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_RD_env, resid=T, pch=16)

#Sex ratio
mgcv_RD_sex <- mgcv::gam(RNADNAratio_mean~s(sexratio), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_sex)
plot(mgcv_RD_sex,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_sex)
par(mfrow=c(2,2))
gam.check(mgcv_RD_sex)
par(mfrow=c(1,1))

#CPUE
mgcv_RD_CPUE <- mgcv::gam(RNADNAratio_mean~s(Kick_mean), family=gaussian(link='identity'),
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_CPUE)
plot(mgcv_RD_CPUE,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_RD_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_RD_CPUE)
par(mfrow=c(1,1))

mgcv_RD_CPUE <- mgcv::gam(RNADNAratio_mean~s(Spread_dist) + s(Kick_mean,k=5), family=gaussian(link='identity'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_RD_CPUE)
par(mfrow=c(1,2))
plot(mgcv_RD_CPUE,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_RD_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_RD_CPUE)
par(mfrow=c(1,1))
concurvity(mgcv_RD_CPUE)

#####NORTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' & !is.na(craydatOR$AFDWmean) & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=RNADNAratio)) + geom_bar() + facet_wrap(~Site)
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20,'RNADNAratio'])
shapiro.test(craydatOR[craydatOR$Site == 25,'RNADNAratio'])
#Test for equal variances
leveneTest(RNADNAratio ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
ttest_RD_sp <- t.test(craydatOR[craydatOR$Site==20,'RNADNAratio'], craydatOR[craydatOR$Site==25,'RNADNAratio'])
ttest_RD_sp


#####SOUTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=RNADNAratio)) + geom_bar() + facet_wrap(~Site)
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'RNADNAratio'])
shapiro.test(craydatOR[craydatOR$Site == 43,'RNADNAratio'])
#Test for equal variances
leveneTest(RNADNAratio ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
wilcox_RD_sp <- wilcox.test(craydatOR[craydatOR$Site==41,'RNADNAratio'], craydatOR[craydatOR$Site==43,'RNADNAratio'])
wilcox_RD_sp

#####ALL SITES#####
gam_RD_sp <-  gam(RNADNAratio_mean~lo(Spread_dist, span=1), family=gaussian(link='log'), 
                  data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$RNADNAratio_mean) & !is.na(craydat_stat$AFDWmean),])
summary(gam_RD_sp)


gam_RD_sp <-  gam(RNADNAratio_mean~lo(meaninvyr, span=1), family=gaussian(link='log'), 
                  data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$RNADNAratio_mean) & !is.na(craydat_stat$AFDWmean),])
summary(gam_RD_sp)


gam_RD_sp <-  gam(RNADNAratio_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                  data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$RNADNAratio_mean) & !is.na(craydat_stat$AFDWmean),])
summary(gam_TP_env)

ggplot(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 5,], aes(x=meaninvyr, y=RNADNAratio_mean))+geom_point()+geom_smooth(span=2)
ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x=meaninvyr, y=trophiclevel_mean))+geom_point()+geom_smooth(span=1)


##############################################################################################

###################################### MODELS FOR CL MEAN ####################################
##### All samples
CL_sex_cray_lme <- lmer(log(CL) ~ Sex  + (1|Site), data=craydatOR, REML=F)
summary(CL_sex_cray_lme)
anova(CL_sex_cray_lme)

ggplot(craydatOR[craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site'],], aes(x=Sex, y=CL)) + geom_boxplot() + facet_wrap(~Site)

ddply(craydatOR[(craydatOR$Site %in% craydat_stat[craydat_stat$ncray>=10,'Site']),], 'Site', function(x) {
  print(c(unique(x$Site), unique(x$Spread_dist)))
  test <- t.test(x[x$Sex=='F','CL'], x[x$Sex=='M','CL'])
  print(test)
})

#####MAINSTEM######
#NULL MODEL
mgcv_CL_inter <- mgcv::gam(CL_mean~1, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_inter)
AIC(mgcv_CL_inter)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_CL_sp <- mgcv::gam(CL_mean~s(Spread_dist, k=4), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_sp)
plot(mgcv_CL_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_sp)
par(mfrow=c(2,2))
gam.check(mgcv_CL_sp)
par(mfrow=c(1,1))


#Check whether the pattern (even if weak) holds with all sites -- it does even though ti is weaker
# mgcv_CL_sp <- mgcv::gam(CL_mean~s(Spread_dist, k=4), family=gaussian(link='identity'),
#                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                             craydat_stat$ncray >10,])
# summary(mgcv_CL_sp)
# plot(mgcv_CL_sp,residuals=TRUE,shade=T, cex=6)

#LOESS smoother for Spread_dist
# gam_CL_sp <-  gam(CL_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
#                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                       craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_CL_sp)
# plot(gam_CL_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_CL_sp'] <- predict(mgcv_CL_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_CL_sp_se'] <- predict(mgcv_CL_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=CL_mean))+
  geom_ribbon(aes(ymin = (gam_CL_sp-1.96*gam_CL_sp_se), ymax = (gam_CL_sp+1.96*gam_CL_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_CL_sp)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for mean invasion year (thin plate regression spline, penalized regression spline)
mgcv_CL_yr <- mgcv::gam(CL_mean~s(meaninvyr, k=4), family=gaussian(link='identity'),
                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_yr)
plot(mgcv_CL_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_yr)
par(mfrow=c(2,2))
gam.check(mgcv_CL_yr)
par(mfrow=c(1,1))

##LOESS smoother for meaninvyr
# gam_CL_yr <-  gam(CL_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
#                   data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                       craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean),])
# summary(gam_CL_yr)
# plot(gam_CL_yr)
# 
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_CL_yr'] <- predict(gam_CL_yr)
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_CL_yr_se'] <- predict(gam_CL_yr, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 5 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=CL_mean))+
  geom_ribbon(aes(ymin = (gam_CL_yr-1.96*gam_CL_yr_se), ymax = (gam_CL_yr+1.96*gam_CL_yr_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_CL_yr)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for environment (thin plate regression spline, penalized regression spline)
mgcv_CL_deg <- mgcv::gam(CL_mean~s(degdays, k=8), family=gaussian(link='identity'),
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_deg)
plot(mgcv_CL_deg,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_deg)
par(mfrow=c(2,2))
gam.check(mgcv_CL_deg)
par(mfrow=c(1,1))

mgcv_CL_AFDW <- mgcv::gam(CL_mean~s(AFDWmean), family=gaussian(link='identity'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_AFDW)
plot(mgcv_CL_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_CL_AFDW)
par(mfrow=c(1,1))

mgcv_CL_green <- mgcv::gam(CL_mean~s(greenmean, k=3), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_green)
plot(mgcv_CL_green,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_green)
par(mfrow=c(2,2))
gam.check(mgcv_CL_green)
par(mfrow=c(1,1))

mgcv_CL_diatom <- mgcv::gam(CL_mean~s(diatommean), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_diatom)
plot(mgcv_CL_diatom,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CL_diatom)
par(mfrow=c(2,2))
gam.check(mgcv_CL_diatom)
par(mfrow=c(1,1))

mgcv_CL_AFDWgreen <- mgcv::gam(CL_mean~s(AFDWmean,k=5) + s(greenmean, k=3), family=gaussian(link='identity'),
                               data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                   craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_AFDWgreen)
par(mfrow=c(1,2))
plot(mgcv_CL_AFDWgreen,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_CL_AFDWgreen)
par(mfrow=c(2,2))
gam.check(mgcv_CL_AFDWgreen)
par(mfrow=c(1,1))
concurvity(mgcv_CL_AFDWgreen)

mgcv_CL_AFDWgreen <- mgcv::gam(CL_mean~s(AFDWmean, k=7) + s(diatommean, k=7), family=gaussian(link='identity'),
                               data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                   craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_AFDWgreen)
par(mfrow=c(1,2))
plot(mgcv_CL_AFDWgreen,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_CL_AFDWgreen)
par(mfrow=c(2,2))
gam.check(mgcv_CL_AFDWgreen)
par(mfrow=c(1,1))
concurvity(mgcv_CL_AFDWgreen)


#LOESS smoother for environment
# #Environment
# gam_CL_env <-  gam(CL_mean~degdays+AFDWmean+lo(greenmean, span=0.8)+diatommean, family=gaussian(link='identity'), 
#                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_CL_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_CL_env, resid=T, pch=16)

mgcv_CL_sex <- lme(CL~sex + , family=gaussian(link='identity'),
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_sex)

#Sex ratio
mgcv_CL_sex <- mgcv::gam(CL_mean~s(sexratio,k=3)+s(Spread_dist), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_sex)
par(mfrow=c(1,2))
plot(mgcv_CL_sex,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_CL_sex)
par(mfrow=c(2,2))
gam.check(mgcv_CL_sex)
par(mfrow=c(1,1))


#CPUE
mgcv_CL_CPUE <- mgcv::gam(CL_mean~s(Kick_mean), family=gaussian(link='identity'),gamma=2,
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_CPUE)
par(mfrow=c(1,2))
plot(mgcv_CL_CPUE,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_CL_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_CL_CPUE)
par(mfrow=c(1,1))

mgcv_CL_CPUE <- mgcv::gam(CL_mean~s(Kick_mean,k=4)+s(Spread_dist,k=4), family=gaussian(link='identity'),gamma=2,
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CL_CPUE)
par(mfrow=c(1,2))
plot(mgcv_CL_CPUE,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_CL_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_CL_CPUE)
par(mfrow=c(1,1))
concurvity(mgcv_CL_CPUE)


#####NF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'North Fork' & craydat_stat$ncray >10,], aes(x=Spread_dist, y=CL_mean))+geom_col()+geom_errorbar(aes(ymin=CL_mean-1.96*(CL_sd/sqrt(ncray)), ymax=CL_mean+1.96*(CL_sd/sqrt(ncray))))
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' &  
                   craydatOR$Site %in%  craydat_stat[craydat_stat$ncray > 10,'Site'],], aes(x=CL)) + geom_histogram()
qplot(craydatOR[craydatOR$River_Tributary == 'North Fork' & !is.na(craydatOR$AFDWmean) & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 5,'Site']),'CL'])

#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20,'CL'])
shapiro.test(craydatOR[craydatOR$Site == 25,'CL'])
#Test for equal variances
leveneTest(CL ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Wilcox t-test
wilcox_CL_sp <- wilcox.test(craydatOR[craydatOR$Site==20,'CL'], craydatOR[craydatOR$Site==25,'CL'])
wilcox_CL_sp


#####SF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'South Fork' & 
                      craydat_stat$ncray > 10,], aes(x=Spread_dist, y=CL_mean))+geom_col()+geom_errorbar(aes(ymin=CL_mean-1.96*(CL_sd/sqrt(ncray)), ymax=CL_mean+1.96*(CL_sd/sqrt(ncray))))
ggplot(craydat_stat[craydat_stat$River_Tributary == 'South Fork' & 
                      craydat_stat$ncray > 10,], aes(x=Spread_dist, y=CL_mean))+geom_point(aes(color=River_Tributary))+geom_smooth(method='lm')

#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'CL'])
shapiro.test(craydatOR[craydatOR$Site == 43,'CL'])
#Test for equal variances
leveneTest(CL ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])


#wilcox
wilcox_CL_sp <- wilcox.test(CL ~ Site, data=craydatOR[craydatOR$River_Tributary == 'South Fork' & !is.na(craydatOR$AFDWmean) &
                                                  (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),],var.equal=F)
wilcox_CL_sp

#####ALL SITES#####
gam_CL_sp <-  gam(CL_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
                  data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CL_sp)
plot(gam_CL_sp)


gam_CL_yr <-  gam(CL_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                  data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CL_yr)
plot(gam_CL_yr)


gam_CL_env <-  gam(CL_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                   data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CL_env)

ggplot(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 5,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=2)
ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=1)







#####DIRICHLET REGRESSION#####
library(DirichletReg)
resp <- DR_data(craydat_stat[craydat_stat$ncray>10 & !is.na(craydat_stat$u15_rel),c("u15_rel","o15u25_rel","o25_rel")])
multinom_glm <- DirichReg(resp~Spread_dist, data=craydat_stat[craydat_stat$ncray>10 & !is.na(craydat_stat$u15_rel),])
coef(multinom_glm)
multinom_glm2 <- DirichReg(resp~1, data=craydat_stat[craydat_stat$ncray>10 & !is.na(craydat_stat$u15_rel),])
anova(multinom_glm, multinom_glm2)
summary(multinom_glm)

ggplot(craydat_stat[craydat_stat$ncray>10 & !is.na(craydat_stat$u15_rel),], aes(Spread_dist)) +
     geom_point(aes(y=u15_rel), color='orange')+
     geom_point(aes(y=o15u25_rel), color='red') +
     geom_point(aes(y=o25_rel), color= 'black') +
     geom_line(aes(y=predict(multinom_glm)[,1]), color='orange') + 
     geom_line(aes(y=predict(multinom_glm)[,2]), color='red') + 
     geom_line(aes(y=predict(multinom_glm)[,3]), color='black') + 
     theme_classic()

##############################################################################################

###################################### MODELS FOR CL SD, MAX, and QT90 #######################
#####MAINSTEM#####
qplot(log(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 10,'CL_sd']))

#NULL MODEL
mgcv_CLsd_inter <- mgcv::gam(CL_sd~1, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLsd_inter)
AIC(mgcv_CLsd_inter)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_CLsd_sp <- mgcv::gam(CL_sd~s(Spread_dist), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLsd_sp)
plot(mgcv_CLsd_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CLsd_sp)
par(mfrow=c(2,2))
gam.check(mgcv_CLsd_sp)
par(mfrow=c(1,1))

#Pattern holds with all sites
# mgcv_CLsd_sp <- mgcv::gam(CL_sd~s(Spread_dist), family=gaussian(link='identity'),
#                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                               craydat_stat$ncray >10,])
# summary(mgcv_CLsd_sp)
# plot(mgcv_CLsd_sp,residuals=TRUE,shade=T, cex=6)

# #LOESS smoother for spread_dist
# gam_CLsd_sp <-  gam(CL_sd~Spread_dist, family=gaussian(link='identity'), 
#                     data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                         craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_CLsd_sp)
#plot(gam_CLsd_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean), 'gam_CLsd_sp'] <- predict(mgcv_CLsd_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean), 'gam_CLsd_sp_se'] <- predict(mgcv_CLsd_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean),], aes(x=Spread_dist, y=CL_sd))+
  geom_ribbon(aes(ymin = (gam_CLsd_sp-1.96*gam_CLsd_sp_se), ymax = (gam_CLsd_sp+1.96*gam_CLsd_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=gam_CLsd_sp), size=1.1) +
  scale_y_continuous(name='Carapace length standard deviation (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 


#TS basis for invasion year (thin plate regression spline, penalized regression spline)
mgcv_CLsd_yr <- mgcv::gam(CL_sd~s(meaninvyr), family=gaussian(link='identity'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLsd_yr)
plot(mgcv_CLsd_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CLsd_yr)
par(mfrow=c(2,2))
gam.check(mgcv_CLsd_yr)
par(mfrow=c(1,1))

#LOESS smoother for invasion year
gam_CLsd_yr <-  gam(CL_sd~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                    data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                        craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean),])
summary(gam_CLsd_yr)
#plot(gam_CLsd_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean), 'gam_CLsd_yr'] <- predict(gam_CLsd_yr)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean), 'gam_CLsd_yr_se'] <- predict(gam_CLsd_yr, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean),], aes(x=meaninvyr, y=CL_sd))+
  geom_ribbon(aes(ymin = (gam_CLsd_yr-1.96*gam_CLsd_yr_se), ymax = (gam_CLsd_yr+1.96*gam_CLsd_yr_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=gam_CLsd_yr), size=1.1) +
  scale_y_continuous(name='Carapace length standard deviation (mm)', breaks=c(10,15,20,25,30)) +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 


# #Environment
# gam_CLsd_env <-  gam(CL_sd~lo(degdays,span=1)+AFDWmean+lo(greenmean, span=1)+diatommean, family=gaussian(link='identity'), 
#                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                          craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean)& !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_CLsd_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_CLsd_env, resid=T, pch=16)

mgcv_CLsd_CPUE <- mgcv::gam(CL_sd~s(Kick_mean,k=5), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLsd_CPUE )
plot(mgcv_CLsd_CPUE ,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CLsd_CPUE )
par(mfrow=c(2,2))
gam.check(mgcv_CLsd_CPUE )
par(mfrow=c(1,1))


###### TEST FOR CL MAX
#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_CLmax_sp <- mgcv::gam(CL_max~s(Spread_dist), family=gaussian(link='identity'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLmax_sp)
plot(mgcv_CLmax_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CLmax_sp)
par(mfrow=c(2,2))
gam.check(mgcv_CLmax_sp)
par(mfrow=c(1,1))

mgcv_CLqt90_sp <- mgcv::gam(CL_qt90~s(Spread_dist), family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_CLqt90_sp)
plot(mgcv_CLqt90_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CLqt90_sp)
par(mfrow=c(2,2))
gam.check(mgcv_CL90_sp)
par(mfrow=c(1,1))

#####NORTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' & !is.na(craydatOR$AFDWmean) & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=CL)) + geom_bar() + facet_wrap(~Site)
#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20,'CL'])
shapiro.test(craydatOR[craydatOR$Site == 25,'CL'])
#Test for equal variances
leveneTest(CL ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#####SOUTH FORK#####
ggplot(craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),], aes(x=CL)) + geom_bar() + facet_wrap(~Site)
ggplot(craydatOR[craydatOR$River_Tributary == 'South Fork' &  
                   craydatOR$Site %in%  craydat_stat[craydat_stat$ncray > 10,'Site'],], aes(x=CL_weight_lm_res)) + geom_histogram()

#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'CL'])
shapiro.test(craydatOR[craydatOR$Site == 43,'CL'])

#Levene's test to test homogeneity of variance in CL among sites
leveneTest(CL ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'South Fork' & 
                                               !is.na(craydatOR$AFDWmean) &
                                               (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#####ALL SITES#####
gam_CLsd_sp <-  gam(CL_sd~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
                    data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CLsd_sp)
plot(gam_CLsd_sp)


gam_CLsd_sp <-  gam(CL_sd~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                    data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CLsd_sp)
plot(gam_CLsd_sp)


gam_CLsd_env <-  gam(CL_sd~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                     data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_CLsd_env)

ggplot(craydat_stat[craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem' & craydat_stat$ncray > 5,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=2)
ggplot(craydat_stat[craydat_stat$ncray > 10,], aes(x=meaninvyr, y=CL_sd))+geom_point()+geom_smooth(span=1)


##############################################################################################

###################################### MODELS FOR weight LENGTH RESIDUALS ####################
#####MAINSTEM######
#NULL MODEL
mgcv_weight_inter <- mgcv::gam(Weight_res_mean~1, family=gaussian(link='identity'),
                           data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_inter)
AIC(mgcv_weight_inter)


#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_weight_sp <- mgcv::gam(Weight_res_mean~s(Spread_dist), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_sp)
plot(mgcv_weight_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_sp)
par(mfrow=c(2,2))
gam.check(mgcv_weight_sp)
par(mfrow=c(1,1))

#Pattern is even stronger with all sites
# mgcv_weight_sp <- mgcv::gam(Weight_res_mean~s(Spread_dist), family=gaussian(link='identity'),
#                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                                 craydat_stat$ncray >10,])
# summary(mgcv_weight_sp)
# plot(mgcv_weight_sp,residuals=TRUE,shade=T, cex=6)


#LOESS smoother for Spread_dist
# gam_weight_sp <-  gam(Weight_res_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
#                       data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                           craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_weight_sp)
# plot(gam_weight_sp)

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_weight_sp'] <- predict(mgcv_weight_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_weight_sp_se'] <- predict(mgcv_weight_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Weight_res_mean))+
  geom_ribbon(aes(ymin = (gam_weight_sp-1.96*gam_weight_sp_se), ymax = (gam_weight_sp+1.96*gam_weight_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_weight_sp)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)') +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for mean invasion year (thin plate regression spline, penalized regression spline)
mgcv_weight_yr <- mgcv::gam(Weight_res_mean~s(meaninvyr), family=gaussian(link='identity'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_yr)
plot(mgcv_weight_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_yr)
par(mfrow=c(2,2))
gam.check(mgcv_weight_yr)
par(mfrow=c(1,1))

#LOESS smoother for meaninvyr
# gam_weight_yr <-  gam(Weight_res_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
#                       data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                           craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean),])
# summary(gam_weight_yr)
# plot(gam_weight_yr)
# 
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_weight_yr'] <- predict(gam_weight_yr)
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_weight_yr_se'] <- predict(gam_weight_yr, se.fit=T)$se.fit
# 
# ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
#                        craydat_stat$River_Tributary == 'Lower mainstem') & 
#                       craydat_stat$ncray > 5 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Weight_res_mean))+
#   geom_ribbon(aes(ymin = (gam_weight_yr-1.96*gam_weight_yr_se), ymax = (gam_weight_yr+1.96*gam_weight_yr_se), fill=River_Tributary), alpha=0.3)+
#   geom_point(aes(color=River_Tributary), size=4)+ 
#   geom_line(aes(y=(gam_weight_yr)), size=1.1) +
#   scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
#   scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
#   #geom_text(aes(label=Site)) +
#   scale_fill_discrete(name='Tributary (95% CI)') + 
#   scale_color_discrete(name='Tributary (95% CI)') +
#   theme_classic() 


#TS basis for environment (thin plate regression spline, penalized regression spline)
mgcv_weight_degdays <- mgcv::gam(Weight_res_mean~s(degdays), family=gaussian(link='identity'),
                                 data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                     craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_degdays)
plot(mgcv_weight_degdays,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_degdays)
par(mfrow=c(2,2))
gam.check(mgcv_weight_degdays)
par(mfrow=c(1,1))

mgcv_weight_AFDW <- mgcv::gam(Weight_res_mean~s(AFDWmean), family=gaussian(link='identity'),
                              data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                  craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_AFDW)
plot(mgcv_weight_AFDW,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_AFDW)
par(mfrow=c(2,2))
gam.check(mgcv_weight_AFDW)
par(mfrow=c(1,1))

mgcv_weight_green <- mgcv::gam(Weight_res_mean~s(greenmean,k=4), family=gaussian(link='identity'),
                               data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                   craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_green)
plot(mgcv_weight_green,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_green)
par(mfrow=c(2,2))
gam.check(mgcv_weight_green)
par(mfrow=c(1,1))

mgcv_weight_diatom <- mgcv::gam(Weight_res_mean~s(diatommean,k=9), family=gaussian(link='identity'),
                                data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                    craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_diatom)
plot(mgcv_weight_diatom,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_diatom)
par(mfrow=c(2,2))
gam.check(mgcv_weight_diatom)
par(mfrow=c(1,1))


mgcv_weight_degAFDW <- mgcv::gam(Weight_res_mean~s(degdays)+s(AFDWmean, k=5), family=gaussian(link='identity'),
                                 data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                     craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_degAFDW)
par(mfrow=c(1,2))
plot(mgcv_weight_degAFDW,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_degAFDW)
par(mfrow=c(2,2))
gam.check(mgcv_weight_degAFDW)
par(mfrow=c(1,1))
concurvity(mgcv_weight_degAFDW)

mgcv_weight_degAFDWgreen <- mgcv::gam(Weight_res_mean~s(degdays)+s(AFDWmean, k=3)+s(greenmean, k=3), family=gaussian(link='identity'),
                                      data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                          craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_degAFDWgreen)
par(mfrow=c(1,3))
plot(mgcv_weight_degAFDWgreen,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_degAFDWgreen)
par(mfrow=c(2,2))
gam.check(mgcv_weight_degAFDWgreen)
par(mfrow=c(1,1))
concurvity(mgcv_weight_degAFDWgreen)



mgcv_weight_deggreen <- mgcv::gam(Weight_res_mean~s(degdays)+s(greenmean, k=3), family=gaussian(link='identity'),
                                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_deggreen)
par(mfrow=c(1,3))
plot(mgcv_weight_deggreen,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_deggreen)
par(mfrow=c(2,2))
gam.check(mgcv_weight_deggreen)
par(mfrow=c(1,1))
concurvity(mgcv_weight_deggreen)

mgcv_weight_degdia <- mgcv::gam(Weight_res_mean~s(degdays)+s(diatommean, k=3), family=gaussian(link='identity'),
                                data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                    craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_degdia)
par(mfrow=c(1,3))
plot(mgcv_weight_degdia,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_degdia)
par(mfrow=c(2,2))
gam.check(mgcv_weight_degdia)
par(mfrow=c(1,1))
concurvity(mgcv_weight_degdia)

mgcv_weight_greendia <- mgcv::gam(Weight_res_mean~s(greenmean, k=3)+s(diatommean, k=3), family=gaussian(link='identity'),
                                  data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_greendia)
par(mfrow=c(1,3))
plot(mgcv_weight_greendia,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_greendia)
par(mfrow=c(2,2))
gam.check(mgcv_weight_greendia)
par(mfrow=c(1,1))
concurvity(mgcv_weight_greendia)


#LOESS smoother for environment
#Environment
# gam_weight_env <-  gam(Weight_res_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
#                        data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                            craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$greenmean) & !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_weight_env) #Not significant
# par(mfrow=c(2,2))
# plot(gam_weight_env, resid=T, pch=16)

#With Sex ratio
mgcv_weight_sex <- mgcv::gam(Weight_res_mean~s(sexratio), family=gaussian(link='identity'),
                                 data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                     craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_sex)
plot(mgcv_weight_sex,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_sex)
par(mfrow=c(2,2))
gam.check(mgcv_weight_sex)
par(mfrow=c(1,1))

#With CPUE
mgcv_weight_CPUE <- mgcv::gam(Weight_res_mean~s(Kick_mean), family=gaussian(link='identity'),
                             data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                 craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_CPUE)
plot(mgcv_weight_CPUE,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_weight_CPUE)
par(mfrow=c(2,2))
gam.check(mgcv_weight_CPUE)
par(mfrow=c(1,1))

mgcv_weight_CPUEsp <- mgcv::gam(Weight_res_mean~s(Kick_mean,k=7)+s(Spread_dist,k=7), family=gaussian(link='identity'),
                              data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                  craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_weight_CPUEsp)
par(mfrow=c(1,2))
plot(mgcv_weight_CPUEsp,residuals=TRUE,shade=T, cex=6)
par(mfrow=c(1,1))
AIC(mgcv_weight_CPUEsp)
par(mfrow=c(2,2))
gam.check(mgcv_weight_CPUEsp)
par(mfrow=c(1,1))
concurvity(mgcv_weight_CPUEsp)



#####NF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'North Fork' & craydat_stat$ncray >10,], aes(x=Spread_dist, y=Weight_res_mean))+geom_col()+
  geom_errorbar(aes(ymin=Weight_res_mean-1.96*Weight_res_se, ymax=Weight_res_mean+1.96*Weight_res_se))
ggplot(craydatOR[craydatOR$River_Tributary == 'North Fork' &  
                   craydatOR$Site %in%  craydat_stat[craydat_stat$ncray > 10,'Site'],], aes(x=CL_weight_lm_res)) + geom_histogram()

#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 20 & is.na(craydatOR$Miss_App),'CL_weight_lm_res'])
shapiro.test(craydatOR[craydatOR$Site == 25,'CL_weight_lm_res'])
#Test for equal variances
leveneTest(CL_weight_lm_res ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'North Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#Welch's t-test
ttest_weight_sp <- t.test(craydatOR[craydatOR$Site==20,'CL'], craydatOR[craydatOR$Site==25,'CL_weight_lm_res'], var.equal = T)
ttest_weight_sp


#####SF######
ggplot(craydat_stat[craydat_stat$River_Tributary == 'South Fork' & 
                      craydat_stat$ncray > 10,], aes(x=Spread_dist, y=Weight_res_mean))+geom_col()+
  geom_errorbar(aes(ymin=Weight_res_mean-1.96*Weight_res_se, ymax=Weight_res_mean+1.96*Weight_res_se))
ggplot(craydatOR[craydatOR$River_Tributary == 'South Fork' &  
                   craydatOR$Site %in%  craydat_stat[craydat_stat$ncray > 10,'Site'],], aes(x=CL_weight_lm_res)) + geom_histogram()

#Test normality of data with Shapiro-Wilk test
shapiro.test(craydatOR[craydatOR$Site == 41,'CL_weight_lm_res'])
shapiro.test(craydatOR[craydatOR$Site == 43,'CL_weight_lm_res'])
#Test for equal variances
leveneTest(CL_weight_lm_res ~ factor(Site), data=craydatOR[craydatOR$River_Tributary == 'South Fork' & (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),])

#t-test with welch's correction
wilcox_weight_sp <- wilcox.test(CL_weight_lm_res ~ Site, data=craydatOR[craydatOR$River_Tributary == 'South Fork' & !is.na(craydatOR$AFDWmean) &                                                                    (craydatOR$Site %in% craydat_stat[craydat_stat$ncray > 10,'Site']),],var.equal=T)
wilcox_weight_sp

#####ALL SITES#####
gam_weight_sp <-  gam(Weight_res_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
                      data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_weight_sp)
plot(gam_weight_sp)


gam_weight_yr <-  gam(Weight_res_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
                      data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_weight_yr)
plot(gam_weight_yr)


gam_weight_env <-  gam(Weight_res_mean~degdays+AFDWmean+greenmean+diatommean, family=gaussian(link='identity'), 
                       data=craydat_stat[craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean),])
summary(gam_weight_env)


##############################################################################################

###################################### MODELS FOR SEX RATIO ##################################
#####MAINSTEM######
cbind(craydat_stat$males, craydat_stat$females)
#NULL MODEL WITH 21 SITES
mgcv_sex_inter<- mgcv::gam(cbind(males, females)~1, family=binomial(link='logit'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10,])
summary(mgcv_sex_inter)
AIC(mgcv_sex_inter)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_sex_sp <- mgcv::gam(cbind(males, females)~s(Spread_dist, k=3), family=binomial(link='logit'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
summary(mgcv_sex_sp)
plot(mgcv_sex_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_sex_sp)
par(mfrow=c(2,2))
gam.check(mgcv_sex_sp)
par(mfrow=c(1,1))


craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_sex_sp'] <- predict(mgcv_sex_sp)
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean), 'gam_sex_sp_se'] <- predict(mgcv_sex_sp, se.fit=T)$se.fit

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),], aes(x=Spread_dist, y=sexratio))+
  geom_ribbon(aes(ymin = (gam_sex_sp-1.96*gam_sex_sp_se), ymax = (gam_sex_sp+1.96*gam_sex_sp_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_sex_sp)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)') +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#With all sites with more than 10 crayfish 
mgcv_sex_sp2 <- mgcv::gam(cbind(males, females)~s(Spread_dist), family=binomial(link='logit'),
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10,])
summary(mgcv_sex_sp2)
plot(mgcv_sex_sp2,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_sex_sp2)
par(mfrow=c(2,2))
gam.check(mgcv_sex_sp2)
par(mfrow=c(1,1))

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10, 'gam_sex_sp2'] <- inv.logit(predict(mgcv_sex_sp2))
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem')& 
               craydat_stat$ncray >10, 'gam_sex_sp2_se'] <- predict(mgcv_sex_sp2, se.fit=T)$se.fit

#LOESS smoother for Spread_dist
# gam_weight_sp <-  gam(Weight_res_mean~lo(Spread_dist, span=1), family=gaussian(link='identity'), 
#                       data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                           craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),])
# summary(gam_weight_sp)
# plot(gam_weight_sp)

ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
                       craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray > 10,], aes(x=Spread_dist, y=sexratio))+
  geom_ribbon(aes(ymin = (gam_sex_sp2-1.96*gam_sex_sp2_se), ymax = (gam_sex_sp2+1.96*gam_sex_sp2_se), fill=River_Tributary), alpha=0.3)+
  geom_point(aes(color=River_Tributary), size=4)+ 
  geom_line(aes(y=(gam_sex_sp2)), size=1.1) +
  scale_y_continuous(name='Carapace length (mm)') +
  scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
  #geom_text(aes(label=Site)) +
  scale_fill_discrete(name='Tributary (95% CI)') + 
  scale_color_discrete(name='Tributary (95% CI)') +
  theme_classic() 

#TS basis for invasion year
mgcv_sex_yr <- mgcv::gam(cbind(males, females)~s(meaninvyr), family=binomial(link='logit'),
                          data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10,])
summary(mgcv_sex_yr)
plot(mgcv_sex_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_sex_yr)
par(mfrow=c(2,2))
gam.check(mgcv_sex_yr)
par(mfrow=c(1,1))

#LOESS smoother for meaninvyr
# gam_weight_yr <-  gam(Weight_res_mean~lo(meaninvyr, span=1), family=gaussian(link='identity'), 
#                       data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                                           craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean),])
# summary(gam_weight_yr)
# plot(gam_weight_yr)
# 
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_weight_yr'] <- predict(gam_weight_yr)
# craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
#                craydat_stat$ncray >5 & !is.na(craydat_stat$AFDWmean), 'gam_weight_yr_se'] <- predict(gam_weight_yr, se.fit=T)$se.fit
# 
# ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | 
#                        craydat_stat$River_Tributary == 'Lower mainstem') & 
#                       craydat_stat$ncray > 5 & !is.na(craydat_stat$AFDWmean),], aes(x=Spread_dist, y=Weight_res_mean))+
#   geom_ribbon(aes(ymin = (gam_weight_yr-1.96*gam_weight_yr_se), ymax = (gam_weight_yr+1.96*gam_weight_yr_se), fill=River_Tributary), alpha=0.3)+
#   geom_point(aes(color=River_Tributary), size=4)+ 
#   geom_line(aes(y=(gam_weight_yr)), size=1.1) +
#   scale_y_continuous(name='Carapace length (mm)', breaks=c(10,15,20,25,30)) +
#   scale_x_continuous(name='Distance from invasion source', expand=c(0,0)) +
#   #geom_text(aes(label=Site)) +
#   scale_fill_discrete(name='Tributary (95% CI)') + 
#   scale_color_discrete(name='Tributary (95% CI)') +
#   theme_classic() 
#####NF#####
sextable_NF <- with(craydatOR[craydatOR$Site == 20 | craydatOR$Site == 25,], table(Sex, Site))
sextable_NF
prop.test(sextable_NF, alternative='two.sided',correct=T) 
#####SF####
sextable_SF <- with(craydatOR[craydatOR$Site == 40 | craydatOR$Site == 43,], table(Sex, Site))
sextable_SF
prop.test(sextable_SF, correct=T) 
##############################################################################################

###################################### MODELS FOR MISSING CHELA ##############################
#####MAINSTEM#######
mgcv_missap_int <- mgcv::gam(cbind(missap, nomissap)~1, family=binomial(link='logit'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10,])
summary(mgcv_missap_int)
AIC(mgcv_missap_int)

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_missap_sp <- mgcv::gam(cbind(missap,nomissap)~s(Spread_dist), family=binomial(link='logit'),
                         data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10,])
summary(mgcv_missap_sp)
plot(mgcv_missap_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_missap_sp)
par(mfrow=c(2,2))
gam.check(mgcv_missap_sp)
par(mfrow=c(1,1))

#TS basis for Spread_dist (thin plate regression spline, penalized regression spline)
mgcv_missap_yr <- mgcv::gam(cbind(missap,nomissap)~s(meaninvyr), family=binomial(link='logit'),
                            data=craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                                craydat_stat$ncray >10,])
summary(mgcv_missap_yr)
plot(mgcv_missap_yr,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_missap_yr)
par(mfrow=c(2,2))
gam.check(mgcv_missap_yr)
par(mfrow=c(1,1))

#####NF#####
with(craydat_stat[craydat_stat$Site == 20 | craydat_stat$Site == 25,], qplot(x=Site, y=Miss_App_tot)) + geom_point()
sextable_NF <- with(craydatOR[craydatOR$Site == 20 | craydatOR$Site == 25,], table(missap_binary, Site))
sextable_NF
prop.test(sextable_NF, alternative='two.sided',correct=T) 
#####SF####
with(craydat_stat[craydat_stat$Site == 40 | craydat_stat$Site == 41 | craydat_stat$Site == 43,], qplot(x=Site, y=Miss_App_tot)) + geom_point()
sextable_SF <- with(craydatOR[craydatOR$Site == 41 | craydatOR$Site == 43,], table(missap_binary, Site))
sextable_SF
prop.test(sextable_SF, correct=T) 


###################################### FIGURES ################################################################################
##########PREPROCESS ################
# Make figures of tributaries using edges_18_updist.shp in HexSim/Crayfish_model/Network_module/Network_18
#Find in python after join

#Get river kilomter of tributary confluences (modify JDNF to match with site 2 even though Site 2 interval is slightly off - to help reader see confluence)
JDNF2 <- 292349
confluences <- c(0, 379.179-c(JDSF/1000, JDNF2/1000))

#Get images of edges
img <- readPNG(file.path(figdir,"Minimaps/mainstem3.png"))
g <- rasterGrob(img, interpolate=TRUE)
img_nf <- readPNG(file.path(figdir,"Minimaps/northfork2.png"))
g_nf <- rasterGrob(img_nf, interpolate=TRUE)
img_sf <- readPNG(file.path(figdir,"Minimaps/southfork2.png"))
g_sf <- rasterGrob(img_sf, interpolate=TRUE)
# img_um <- readPNG(file.path(figdir,"Minimaps/upstmain.png"))
# g_um <- rasterGrob(img_um, interpolate=TRUE)


##########CHELA LENGTH FIGURE###########
Chela_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                              craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) &!is.na(craydat_stat$trophiclevel_mean),], 
                               aes(x = Spread_dist, y = Chelae_res_mean)) +
  geom_ribbon(aes(ymin = (gam_chela_sp-1.96*gam_chela_sp_se), ymax = (gam_chela_sp+1.96*gam_chela_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_linerange(aes(ymin= Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=0.75, size=1, color='#67a9cf') +
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_chela_sp)), size=1, alpha=0.75, color='#016c59') +
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=0, label='Downstream', angle=90, size=2.5, color='#016c59')+
  annotate('text', x=160, y=-2.4, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=0.5, ymax=3) +
  #geom_text(aes(label=Site)) +
  #geom_smooth(span=1, color='black') + 
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(-3,3), expand=c(0,0), breaks=c(-3,0,3),name='Relative chela length (mm)') + 
  theme(axis.title.y = element_text(hjust=0.05, vjust=-1.25, size=9),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0,0.5,0,0), "cm"))

#Chela_spreadist_main

Chelasegdat <- data.frame(x1=90, x2=112, y1=2.3, y2=2.3)
Chela_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                               (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),], 
                                aes(x = Spread_dist, y = Chelae_res_mean)) +
  geom_linerange(aes(ymin=Chelae_res_mean-1.96*Chelae_res_se, ymax=Chelae_res_mean+1.96*Chelae_res_se), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+ 
  theme_classic() + 
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1, alpha=0.75) + 
  annotate('text', x=260,y=0, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  annotation_custom(g_nf, xmin=105, xmax=155, ymin=0.5, ymax=3) +
  annotation_custom(g_sf, xmin=40, xmax=90, ymin=0.5, ymax=3) +
  annotate('text', x=83, y=-1, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=-1, label='South Fork', angle=90, size=2.5) +
  geom_segment(data=Chelasegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=101, y=2.5, label='***', size=4)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(-3,3), expand=c(0,0), breaks=c(0,3),name='Carapace length (mm)') + 
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        plot.margin=unit(c(0.25,0.5,0,0), "cm"),
        axis.title.y = element_blank())
#Chela_spreadist_upst
# 
# grid.newpage()
# grid.draw(rbind(ggplotGrob(Chela_spreadist_upst), ggplotGrob(Chela_spreadist_main), size='last'))
##########TP Figure###########
#No obvious bias in mesohabitat sampling that seems to correlate with trophic position (check out cray_TL991_mesohab_trib.png and site_trophiclevel0_mesohabitatspread.png)

#Get river kilometer of tributary confluences
TP_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                           craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$trophiclevel_mean),], 
                            aes(x = Spread_dist, y = trophiclevel_mean)) +

  geom_linerange(aes(ymin=trophiclevel_mean-trophiclevel_SE, ymax=trophiclevel_mean+trophiclevel_SE), alpha=0.75, size=1, color='#67a9cf') +
  geom_ribbon(aes(ymin = (gam_TP_sp-1.96*gam_TP_sp_se), ymax = (gam_TP_sp+1.96*gam_TP_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_TP_sp)), size=1, color='#016c59') +
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=3.5, label='Downstream', angle=90, size=2.5, color='#016c59') +
  #geom_text(aes(label=Site), nudge_y=-0.2)+
  annotate('text', x=160, y=2.3, label='Mainstem', size=2.5) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(1.5,5), expand=c(0,0), name='Trophic position') + 
  theme(plot.margin=unit(c(0,0.1,0.05,0.05), "cm"), 
        axis.title.y = element_text(hjust=-1.4, vjust=-1.25, size=9),
        axis.title.x = element_blank())
#TP_spreadist_main

TPsegdat <- data.frame(x1=90, x2=112, y1=3.5, y2=3.5)

TP_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                              (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),], 
                               aes(x = Spread_dist, y = trophiclevel_mean)) +
  geom_linerange(aes(ymin=trophiclevel_mean-trophiclevel_SE, ymax=trophiclevel_mean+trophiclevel_SE), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+ 
  #geom_text(aes(label=Site), nudge_y=-0.2)+
  theme_classic() + 
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1) + 
  annotate('text', x=83, y=3.1, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=3.1, label='South Fork', angle=90, size=2.5) +
  annotate('text', x=260,y=3.5, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  geom_segment(data=TPsegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=101, y=3.6, label=c('**'), size=4)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(1.5,5.2), expand=c(0,0), breaks=c(3,4,5),name='Trophic position') + 
  geom_vline(xintercept = confluences) +
  theme(plot.margin=unit(c(0.25,0,-0.1,0), "cm"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
#TP_spreadist_upst

#grid.newpage()
#grid.draw(rbind(ggplotGrob(TP_spreadist_upst), ggplotGrob(TP_spreadist_main)))

##########RNADNA Figure###########
#Get river kilometer of tributary confluences
RNADNA_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) & !is.na(craydat_stat$trophiclevel_mean),], 
                                aes(x = Spread_dist, y = RNADNAratio_mean)) +
  geom_linerange(aes(ymin=RNADNAratio_mean-RNADNAratio_SE, ymax=RNADNAratio_mean+RNADNAratio_SE), alpha=0.75, size=1, color='#67a9cf') +
  geom_ribbon(aes(ymin =(gam_RD_sp-1.96*gam_RD_sp_se), ymax = (gam_RD_sp+1.96*gam_RD_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_RD_sp)), size=1, color='#016c59') +
  geom_vline(xintercept = confluences) +
  #geom_text(aes(label=Site), nudge_y=-0.2)+
  annotate('text', x=260,y=3, label='Downstream', angle=90, size=2.5, color='#016c59') +
  annotate('text', x=160, y=0.7, label='Mainstem', size=2.5) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(0,7), expand=c(0,0), breaks=c(0,2,4,6), name='Physiological fitness (RNA/DNA)') + 
  theme(plot.margin=unit(c(0,0,0.1,0), "cm"), 
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(vjust=0.5, hjust=-0.1, size=9))
RNADNA_spreadist_main

RDsegdat <- data.frame(x1=c(57, 90), x2=c(71, 112), y1=c(5.8, 3.4), y2=c(5.8, 3.4))

RNADNA_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                              (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),], 
                               aes(x = Spread_dist, y = RNADNAratio_mean)) +
  geom_linerange(aes(ymin=RNADNAratio_mean-RNADNAratio_SE, ymax=RNADNAratio_mean+RNADNAratio_SE), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+ 
  #geom_text(aes(label=Site), nudge_y=-0.2)+
  theme_classic() + 
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1) +   
  annotate('text', x=260,y=3.5, label='Upstream', angle=90, size=2.5, color='#045a8d') +
  annotate('text', x=83, y=2.2, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=2.2, label='South Fork', angle=90, size=2.5) +
  geom_segment(data=RDsegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=c(64,101), y=c(6, 3.6), label=c('**', '**'), size=4)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(0,7), expand=c(0,0), breaks=c(0,2,4,6),name='Physiological fitness (RNA/DNA)') + 
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0.25,0,-0.1,0), "cm"))
#RNADNA_spreadist_upst

#grid.newpage()
#grid.draw(rbind(ggplotGrob(RNADNA_spreadist_upst), ggplotGrob(RNADNA_spreadist_main)))

##########EXPORT FIGURE ###############
pdf(file.path(figdir,'Trends/Chela_TP_RNADNA_7.pdf'), width=3.14961, height=5.90551)
grid.arrange(rbind(ggplotGrob(Chela_spreadist_upst), ggplotGrob(Chela_spreadist_main),
                ggplotGrob(TP_spreadist_upst), ggplotGrob(TP_spreadist_main),
                ggplotGrob(RNADNA_spreadist_upst), ggplotGrob(RNADNA_spreadist_main), size='last'))
dev.off()

#############################################################################################

##########CL Figure###########
CL_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                           craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) &  !is.na(craydat_stat$trophiclevel_mean),], 
                            aes(x = Spread_dist, y = CL_mean)) +
  geom_ribbon(aes(ymin = (gam_CL_sp-1.96*gam_CL_sp_se), ymax = (gam_CL_sp+1.96*gam_CL_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_linerange(aes(ymin=CL_mean-1.96*CL_SE, ymax=CL_mean+1.96*CL_SE), alpha=0.75, size=1, color='#67a9cf') +
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_CL_sp)), size=1, alpha=0.75, color='#016c59')+
  geom_vline(xintercept = confluences) +
  annotate('text', x=160, y=13, label='Mainstem', size=2.5) +
  annotate('text', x=260,y=25, label='Downstream', angle=90, size=2.5, color='#016c59') +
  annotation_custom(g, xmin=200, xmax=250, ymin=27.5, ymax=42) +
  theme_classic()+ 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(10,40), expand=c(0,0), , breaks=c(10,20,40), name='Carapace Length (mm)') + 
  theme(plot.margin=unit(c(0,0.05,0,0), "cm"), 
        axis.title.y = element_text(hjust=1.25),
        axis.title.x = element_blank())
#CL_spreadist_main


CL_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                            (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),], 
                             aes(x = Spread_dist, y = CL_mean)) +
  geom_linerange(aes(ymin=CL_mean-1.96*CL_SE, ymax=CL_mean+1.96*CL_SE), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+   
  theme_classic() + 
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1) + 
  annotation_custom(g_nf, xmin=105, xmax=155, ymin=27.5, ymax=42) +
  annotation_custom(g_sf, xmin=40, xmax=90, ymin=27.5, ymax=42) +
  annotate('text', x=83, y=15, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=15, label='South Fork', angle=90, size=2.5) +
  annotate('text', x=260,y=25, label='Upstream', angle=90, size=2.5, color='#045a8d') +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(10,42), expand=c(0,0), , breaks=c(20,40), name='Carapace length (mm)') + 
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        plot.margin=unit(c(0.25,0,-0.1,0), "cm"),
        axis.title.y = element_blank())
#CL_spreadist_upst

#grid.newpage()
#grid.draw(rbind(ggplotGrob(CL_spreadist_upst), ggplotGrob(CL_spreadist_main)))

# pdf(file.path(figdir,'Trends/CL_sp.pdf'), width=5, height=3)
# grid.draw(rbind(ggplotGrob(CL_spreadist_upst), ggplotGrob(CL_spreadist_main)))
# dev.off()
##########CLsd Figure ######
CLsdsegdat <- data.frame(x1=c(57, 90), x2=c(71, 112), y1=c(9, 10), y2=c(9, 10))

CLsd_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                                             craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) &  !is.na(craydat_stat$trophiclevel_mean),], 
                              aes(x = Spread_dist, y = CL_sd)) +
  geom_ribbon(aes(ymin =(gam_CLsd_sp-1.96*gam_CLsd_sp_se), ymax = (gam_CLsd_sp+1.96*gam_CLsd_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_linerange(aes(ymin=CL_sd, ymax=CL_sd), alpha=0.75, size=1, color='#67a9cf') +
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_CLsd_sp)), size=1, alpha=0.75, color='#016c59')+
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=8, label='Downstream', angle=90, size=2.5, color='#016c59') +
  annotate('text', x=160, y=4.8, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=28, ymax=48) +
  #geom_text(aes(label=Site)) +
  #geom_smooth(method='gam',alpha=1/3, se=F) + 
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0),breaks=c(0,100,200), name='Distance from invasion source') +
  scale_y_continuous(limits=c(3,12), expand=c(0,0),breaks=c(0,5,10), name='CL standard deviation (mm)') + 
  coord_cartesian(ylim=c(4,12)) +
  theme(plot.margin=unit(c(0,0.05,0,0), "cm"), 
        axis.title.y=element_text(hjust=1.25),
        axis.title.x = element_blank())
#CLsd_spreadist_main

CLsd_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                              (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),], 
                               aes(x = Spread_dist, y = CL_sd)) +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+   
  theme_classic() + 
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1) + 
  annotate('text', x=83, y=6, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=6, label='South Fork', angle=90, size=2.5) +
  annotate('text', x=260,y=8, label='Upstream', angle=90, size=2.5, color='#045a8d') +
  geom_segment(data=CLsdsegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=101, y=10.5, label='**', size=4)+
  annotate('text', x=64, y=11, label='.', size=6)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0),breaks=c(0,100,200) , name='Distance from invasion source') +
  scale_y_continuous(limits=c(3,12), expand=c(0,0),breaks=c(0,5,10), name='CL standard deviation (mm)') + 
  coord_cartesian(ylim=c(4,12)) +
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        plot.margin=unit(c(0.25,0,-0.1,0), "cm"),
        axis.title.y = element_blank())
#CLsd_spreadist_upst
#g1<-ggplot_gtable(ggplot_build(CLsd_spreadist_main))
#g2<-ggplot_gtable(ggplot_build(CLsd_spreadist_upst))

##########WEIGHT FIGURE###########
weight_spreadist_main <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
                                               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean) &  !is.na(craydat_stat$trophiclevel_mean),],
                                aes(x = Spread_dist, y = Weight_res_mean)) +
  geom_ribbon(aes(ymin = (gam_weight_sp-1.96*gam_weight_sp_se), ymax = (gam_weight_sp+1.96*gam_weight_sp_se)), alpha=0.75, fill='#bdc9e1')+
  geom_linerange(aes(ymin=Weight_res_mean-1.96*Weight_res_se, ymax=Weight_res_mean+1.96*Weight_res_se), alpha=0.75, size=1, color='#67a9cf') +
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_weight_sp)), size=1, alpha=0.75, color='#016c59')+
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=0, label='Downstream', angle=90, size=2.5, color='#016c59') +
  annotate('text', x=160, y=-0.27, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=28, ymax=48) +
  theme_classic() +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(-0.3,0.3), expand=c(0,0), breaks=c(-0.3,0,0.3),name='Carapace length (mm)') +
  theme(plot.margin=unit(c(0,0.05,0,0), "cm"),
        axis.title.y = element_text(hjust=1.25),
        axis.title.x = element_blank())
#weight_spreadist_main

# 
weight_spreadist_upst  <- ggplot(craydat_stat[craydat_stat$ncray >=10  & !is.na(craydat_stat$AFDWmean) &
                                                (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),],
                                 aes(x = Spread_dist, y = Weight_res_mean)) +
  geom_linerange(aes(ymin=Weight_res_mean-1.96*Weight_res_se, ymax=Weight_res_mean+1.96*Weight_res_se), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+   
  theme_classic() +
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1) +
  annotate('text', x=260,y=0, label='Upstream', angle=90, size=2.5,  color='#045a8d') +
  # annotation_custom(g_nf, xmin=85, xmax=135, ymin=32, ymax=52) +
  # annotation_custom(g_sf, xmin=45, xmax=95, ymin=32, ymax=52) +
  # annotation_custom(g_um, xmin=3, xmax=53, ymin=32, ymax=52) +
  annotate('text', x=83, y=-0.2, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=-0.2, label='South Fork', angle=90, size=2.5) +
  # annotate('text', x=4, y=49, label='MS', angle=90, size=3) +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(-0.3,0.3), expand=c(0,0), breaks=c(-0.3,0,0.3),name='Carapace length (mm)') +
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        plot.margin=unit(c(0.25,0,-0.1,0), "cm"),
        axis.title.y = element_blank())
#weight_spreadist_upst
# 
# grid.newpage()
# grid.draw(rbind(ggplotGrob(weight_spreadist_upst), ggplotGrob(weight_spreadist_main)))








##########EXPORT FIGURE #######
pdf(file.path(figdir,'Trends/CL_CLsd_weight.pdf'), width=3.14961, height=5.90551)
grid.arrange(rbind(ggplotGrob(CL_spreadist_upst), ggplotGrob(CL_spreadist_main),
                   ggplotGrob(CLsd_spreadist_upst), ggplotGrob(CLsd_spreadist_main),
                   ggplotGrob(weight_spreadist_upst), ggplotGrob(weight_spreadist_main)))
dev.off()
#############################################################################################
#########CPUE FIGURE###########
crayabs <- habdatdistinfo[habdatdistinfo$OR_PresAbs=='Absence' & !is.na(habdatdistinfo$Spread_dist), c('Site', 'River_Tributary', 'Spread_dist')]
crayabs$Spread_dist <- crayabs$Spread_dist/1000
crayabs$Kick_mean <- 0
crayabs$Kick_std <- 0
crayabs$ncray <- 0

craydat_CPUEplot <- rbind(craydat_stat[,c(colnames(crayabs))], crayabs)

craydat_CPUEplot$plotrib <- craydat_CPUEplot$River_Tributary
craydat_CPUEplot[craydat_CPUEplot$River_Tributary == 'Lower mainstem' | craydat_CPUEplot$River_Tributary == 'Upper mainstem' |
                       craydat_CPUEplot$River_Tributary == 'Upst Upper mainstem','plotrib'] <- 'Mainstem'

craydat_CPUEplot <- craydat_CPUEplot[(craydat_CPUEplot$plotrib=='Mainstem' | craydat_CPUEplot$River_Tributary=='North Fork' | craydat_CPUEplot$River_Tributary=='South Fork') &
                                     !(craydat_CPUEplot$Site %in% c(37,104,105)),]
craydat_CPUEplot$SE_H <- with(craydat_CPUEplot, Kick_mean+1.96*Kick_std/sqrt(6))
craydat_CPUEplot$SE_L <- with(craydat_CPUEplot, Kick_mean-1.96*Kick_std/sqrt(6))


#Get map for plot

img <- readPNG(file.path(figdir,"Minimaps/CPUE_minimap2.png"))
gcpue <- rasterGrob(img, interpolate=TRUE)


#All sites 
mgcv_CPUE_sp <- mgcv::gam(Kick_mean~s(Spread_dist), family=gaussian(link='identity'),
                         data=craydat_CPUEplot)
summary(mgcv_CPUE_sp)
plot(mgcv_CPUE_sp,residuals=TRUE,shade=T, cex=6)
AIC(mgcv_CPUE_sp)
par(mfrow=c(2,2))
gam.check(mgcv_CPUE_sp)
par(mfrow=c(1,1))
craydat_CPUEplot[, 'gam_CPUE_sp'] <- predict(mgcv_CPUE_sp)
craydat_CPUEplot[, 'gam_CPUE_sp_seL'] <- predict(mgcv_CPUE_sp)-1.96*predict(mgcv_CPUE_sp, se.fit=T)$se.fit
craydat_CPUEplot[, 'gam_CPUE_sp_seU'] <- predict(mgcv_CPUE_sp)+1.96*predict(mgcv_CPUE_sp, se.fit=T)$se.fit
craydat_CPUEplot[craydat_CPUEplot$gam_CPUE_sp_seL<0, 'gam_CPUE_sp_seL'] <- 0
craydat_CPUEplot[craydat_CPUEplot$gam_CPUE_sp <0, 'gam_CPUE_sp'] <- 0
craydat_CPUEplot$plotrib <- factor(craydat_CPUEplot$plotrib, levels=c("Mainstem", "South Fork", "North Fork"), labels=c("Mainstem", "South Fork", "North Fork"))


CPUE_spreadist <- ggplot(craydat_CPUEplot,aes(x = Spread_dist, y = Kick_mean)) +
  geom_ribbon(aes(ymin = gam_CPUE_sp_seL, ymax = gam_CPUE_sp_seU), alpha=0.5, fill='#bdc9e1')+
  geom_line(aes(y=(gam_CPUE_sp)), size=1, alpha=0.75)+
  geom_linerange(data=craydat_CPUEplot[craydat_CPUEplot$ncray>0,],aes(ymin=Kick_mean-1.96*(Kick_std/sqrt(ncray)), ymax=Kick_mean+1.96*(Kick_std/sqrt(ncray)),color=plotrib),size=1.5, alpha=1/3) +
  geom_point(data=craydat_CPUEplot[craydat_CPUEplot$ncray>0,],size=3, alpha=0.75, aes(color=plotrib))+
  scale_color_discrete(name='Tributary') +
  geom_point(data=craydat_CPUEplot[craydat_CPUEplot$ncray==0,],aes(color=plotrib) ,size=3, alpha=0.75,shape=13)+
  annotate('text', x=81, y=5.5,label='NF conflu.', angle=90, size=2) +
  annotate('text', x=38, y=5.5, label='SF conflu.', angle=90, size=2) +
  geom_vline(xintercept = confluences) +
  #geom_smooth(span=0.5) +
  theme_classic() +
  scale_x_continuous(limits=c(0,320), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-20,45), expand=c(0.02,0),breaks=c(0,10,20,30,40),name=expression(paste('Kick-seining CPUE (crayfish/',~m^2,')',sep=""))) +
  geom_hline(yintercept=0) +
  coord_cartesian(y=c(0,40)) +
  theme(plot.margin=unit(c(0.5,0.2,0,0), "cm"),
        legend.position = c(0.6,0.6),
        axis.line.x= element_blank(),
        text= element_text(size=8))
CPUE_spreadist

pdf(file.path(figdir,'Trends/CPUE_small_expand1.pdf'), width=3.14961, height=2)
CPUE_spreadist
dev.off()


CPUE_spreadist <- ggplot(craydat_CPUEplot,
                         aes(x = Spread_dist, y = Kick_mean)) +
  geom_ribbon(aes(ymin = gam_CPUE_sp_seL, ymax = gam_CPUE_sp_seU), alpha=0.3, fill='#bdc9e1')+
  geom_line(aes(y=(gam_CPUE_sp)), size=1, alpha=0.75, color='black')+
  geom_linerange(data=craydat_CPUEplot[craydat_CPUEplot$ncray>0,],aes(ymin=Kick_mean-1.96*(Kick_std/sqrt(ncray)), ymax=Kick_mean+1.96*(Kick_std/sqrt(ncray)),color=plotrib),size=1.5, alpha=1/3) +
  geom_point(data=craydat_CPUEplot[craydat_CPUEplot$ncray>0,],size=3, alpha=0.75, aes(color=plotrib))+
  scale_color_manual(name='', values=c('#1b9e77','#7570b3', '#d95f02')) +
  geom_point(data=craydat_CPUEplot[craydat_CPUEplot$ncray==0,],aes(color=plotrib) ,size=3, alpha=0.75,shape=13)+
  annotate('text', x=81, y=36,label='North Fork', angle=90, size=4) +
  annotate('text', x=38, y=36, label='South Fork', angle=90, size=4) +
  annotation_custom(gcpue, xmin=240, xmax=290, ymin=28, ymax=43) +
  geom_vline(xintercept = confluences, linetype = "dashed") +
  #geom_smooth(span=0.5) +
  theme_classic() +
  scale_x_continuous(limits=c(0,320), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-20,45), expand=c(0.02,0),breaks=c(0,10,20,30,40),name=expression(paste('Kick-seining CPUE (crayfish/',~m^2,')',sep=""))) +
  geom_hline(yintercept=0) +
  coord_cartesian(y=c(0,40)) +
  theme(plot.margin=unit(c(0.5,0.2,0,0), "cm"),
        legend.position = c(0.85,0.6),
        axis.line.x= element_blank(),
        text= element_text(size=12),
        legend.text = element_text(size=13),
        legend.background = element_blank())
CPUE_spreadist

pdf(file.path(figdir,'Trends/CPUE_small_expand1large_20181008.pdf'), width=6, height=4)
CPUE_spreadist
dev.off()

#########SEX RATIO FIGURE###########
sex_spreadist_main2 <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') &
                                            craydat_stat$ncray >10,],
                             aes(x = Spread_dist, y = sexratio)) +
  geom_ribbon(aes(ymin = (gam_sex_sp2-1.96*gam_sex_sp2_se), ymax = (gam_sex_sp2+1.96*gam_sex_sp2_se)), alpha=0.75, fill='#bdc9e1')+
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_line(aes(y=(gam_sex_sp2)), size=1, alpha=0.75, color='#016c59')+
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=0.55, label='Downstream', angle=90, size=2.5, color='#016c59') +
  annotate('text', x=160, y=0.3, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=28, ymax=48) +
  theme_classic() +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-2,2), expand=c(0,0), breaks=c(0.2, 0.4,0.6,0.8), labels=c(20,40,60,80),name='% males') +
  coord_cartesian(ylim=c(0.2,0.8)) +
  theme(plot.margin=unit(c(0,0.05,0,0), "cm"),
        axis.title.y = element_text(hjust=1.2))

sex_spreadist_upst2  <- ggplot(craydat_stat[craydat_stat$ncray >=10 &
                                             (craydat_stat$River_Tributary == 'North Fork' | craydat_stat$River_Tributary == 'South Fork'),],
                              aes(x = Spread_dist, y = sexratio)) +
  #geom_linerange(aes(ymin=Weight_res_mean-1.96*Weight_res_se, ymax=Weight_res_mean+1.96*Weight_res_se), alpha=0.75, size=1, color='#74a9cf') +
  geom_point(size=3, alpha=0.75,color='#2b8cbe')+   
  theme_classic() +
  geom_smooth(method='lm', aes(group=River_Tributary), color='#045a8d', size=1, se=F) +
  annotate('text', x=260,y=0.55, label='Upstream', angle=90, size=2.5,  color='#045a8d') +
  # annotation_custom(g_nf, xmin=85, xmax=135, ymin=32, ymax=52) +
  # annotation_custom(g_sf, xmin=45, xmax=95, ymin=32, ymax=52) +
  # annotation_custom(g_um, xmin=3, xmax=53, ymin=32, ymax=52) +
  annotate('text', x=83, y=0.4, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=0.4, label='South Fork', angle=90, size=2.5) +
  # annotate('text', x=4, y=49, label='MS', angle=90, size=3) +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from invasion source') +
  scale_y_continuous(limits=c(-2,2), expand=c(0,0), breaks=c(0.4,0.6,0.8),labels=c(40,60,80),name='% males') +
  coord_cartesian(ylim=c(0.2,0.8)) +
  geom_vline(xintercept = confluences) +
  theme(axis.title.x = element_blank(),
        axis.text.x= element_blank(),
        plot.margin=unit(c(0.25,0,-0.1,0), "cm"),
        axis.title.y = element_blank())

# grid.newpage()
# grid.draw(rbind(ggplotGrob(sex_spreadist_upst), ggplotGrob(sex_spreadist_main), ggplotGrob(sex_spreadist_main2)))

pdf(file.path(figdir,'Trends/sexratio.pdf'), width=3.14961, height=2.952755)
grid.arrange(rbind(ggplotGrob(sex_spreadist_upst2),
                   ggplotGrob(sex_spreadist_main2), size='last'))
dev.off()

#########CHELA LENGTH + SPREAD MODEL DIST FIGURE #############
mgcv_chela_sexsp_pred <- predict(mgcv_chela_sexsp, type="terms", se.fit=TRUE)[[1]]
mgcv_chela_sexsp_se <- predict(mgcv_chela_sexsp, type="terms", se.fit=TRUE)[[2]]

craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),'gam_chela_sex_pred'] <- mgcv_chela_sexsp_pred[,1]
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),'gam_chela_sp_pred'] <- mgcv_chela_sexsp_pred[,2]
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),'gam_chela_sex_se'] <- mgcv_chela_sexsp_se[,1]
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),'gam_chela_sp_se'] <- mgcv_chela_sexsp_se[,2]
craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
               craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),'gam_chela_sexsp_res'] <- residuals(mgcv_chela_sexsp)


chelasexsp_gam1  <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),],
       aes(x=sexratio, y=gam_chela_sex_pred)) +
  geom_ribbon(aes(ymin=gam_chela_sex_pred-1.96*gam_chela_sex_se, ymax=gam_chela_sex_pred+1.96*gam_chela_sex_se), alpha=0.75, fill='#bdc9e1') +
  geom_line(size=1, alpha=0.75, color='#016c59') + 
  geom_point(aes(y=gam_chela_sex_pred+gam_chela_sexsp_res),size=3, alpha=0.75, color='#1c9099') +
  scale_y_continuous(name='Chela length partial residuals (mm)', limits=c(-3.1,3), breaks=c(-3,-2,-1,0,1,2,3), expand=c(0,0)) +
  scale_x_continuous(name='% males', limits=c(0.33,0.8), breaks=c(seq(0.4,0.8,0.1)), expand=c(0,0)) +
  coord_cartesian(ylim=c(-2,1.5)) +
  theme_classic()

chelasexsp_gam2 <- ggplot(craydat_stat[(craydat_stat$River_Tributary == 'Upper mainstem' | craydat_stat$River_Tributary == 'Lower mainstem') & 
                      craydat_stat$ncray >10 & !is.na(craydat_stat$AFDWmean)& !is.na(craydat_stat$trophiclevel_mean),],
       aes(x=Spread_dist, y=gam_chela_sp_pred)) +
  geom_ribbon(aes(ymin=gam_chela_sp_pred-1.96*gam_chela_sp_se, ymax=gam_chela_sp_pred+1.96*gam_chela_sp_se), alpha=0.75, fill='#bdc9e1') +
  geom_line(size=1, alpha=0.75, color='#016c59') + 
  geom_point(aes(y=gam_chela_sp_pred+gam_chela_sexsp_res),size=3, alpha=0.75, color='#1c9099') + 
  scale_y_continuous(name='Chela length partial residuals (mm)', limits=c(-3.1,3), breaks=c(-3,-2,-1,0,1,2,3), expand=c(0,0)) +
  scale_x_continuous(name='Distance from initial introduction (km)', limits=c(0,250)) +
  coord_cartesian(ylim=c(-2,1.5)) +
  theme_classic()

pdf(file.path(figdir,'Trends/chela_sexsp.pdf'), width=3.14961, height=5.8)
grid.arrange(chelasexsp_gam1,chelasexsp_gam2)
dev.off()

#############################################################################################
#########ENVIRONMENTAL VARIABLES FIGURE########################
habdatdistinfo$plotrib <- as.character(habdatdistinfo$River_Tributary)
habdatdistinfo[(habdatdistinfo$River_Tributary == 'Lower mainstem' | habdatdistinfo$River_Tributary == 'Upper mainstem')
               & !is.na(habdatdistinfo$River_Tributary),'plotrib'] <- 'Mainstem'
habdatdistinfo$plotrib <- factor(habdatdistinfo$plotrib, levels=c("Mainstem", "South Fork", "North Fork"), labels=c("Mainstem", "South Fork", "North Fork"))

############ DEGREE DAYS ##################
degdays_spreadist_main <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'Lower mainstem' | 
                                                                                     habdatdistinfo$River_Tributary == 'Upper mainstem'),], 
                                 aes(x = Spread_dist/1000, y = degdays)) +
  geom_smooth(aes(group=River_Tributary), color='darkgrey', alpha=0.3, span=1) + 
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=4400, label='Downstream', angle=90, size=2.5, color='#016c59')+
  annotate('text', x=160, y=3650, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=0.5, ymax=3) +
  #geom_text(aes(label=Site)) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(3500,5300), expand=c(0,0), breaks=c(3500,4000,4500, 5000), labels=c(3.5, 4.0, 4.5,5.0),name='Degree days (10???C)') + 
  theme(axis.title.x = element_blank(),
        plot.margin=unit(c(-0.1,0.5,0,0), "cm"))
#degdays_spreadist_main 

degdays_spreadist_upst <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'North Fork' | 
                                                                                     habdatdistinfo$River_Tributary == 'South Fork'),], 
                                 aes(x = Spread_dist/1000, y = degdays)) +
  geom_smooth(method='lm',aes(group=River_Tributary), color='darkgrey', alpha=0.3) + 
  geom_linerange(aes(ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_point(size=3, alpha=0.75, color='#2b8cbe')+ 
  geom_vline(xintercept = confluences) +
  theme_classic() + 
  annotate('text', x=260,y=4400, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  # annotation_custom(g_nf, xmin=105, xmax=155, ymin=0.5, ymax=3) +
  # annotation_custom(g_sf, xmin=40, xmax=90, ymin=0.5, ymax=3) +
  annotate('text', x=83, y=3800, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=3800, label='South Fork', angle=90, size=2.5) +
  geom_segment(data=Chelasegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=101, y=2.5, label='***', size=4)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(3500,5300), expand=c(0,0), breaks=c(3500, 4000,4500,5000), labels=c(3.5, 4.0,4.5, 5.0),name='Degree days (???C)') + 
  geom_vline(xintercept = confluences) +
  theme(plot.margin=unit(c(0.25,0.5,0,0), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

#grid.draw(rbind(ggplotGrob(degdays_spreadist_upst),ggplotGrob(degdays_spreadist_main), size='last'))

############# AFDW ############
AFDW_spreadist_main <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'Lower mainstem' | 
                                                                                   habdatdistinfo$River_Tributary == 'Upper mainstem'),], 
                               aes(x = Spread_dist/1000, y = AFDWmean)) +
  geom_smooth(method='lm',aes(group=River_Tributary), color='darkgrey', alpha=0.3) + 
  geom_linerange(aes(ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn))), alpha=0.6, size=1, color='#67a9cf') +
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=0.25, label='Downstream', angle=90, size=2.5, color='#016c59')+
  annotate('text', x=160, y=0.04, label='Mainstem', size=2.5) +
  annotation_custom(g, xmin=200, xmax=250, ymin=0.5, ymax=3) +
  #geom_text(aes(label=Site)) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-1,1), expand=c(0,0), breaks=c(0,0.20,0.40),name='Macroinvertebrate biomass (AFDW, g)') + 
  coord_cartesian(ylim=c(0,0.6)) + 
  theme(axis.title.y = element_text(hjust=0.05, vjust=-1.25, size=9),
        axis.title.x = element_text(size=9),
        plot.margin=unit(c(-0.1,0.5,0,0), "cm"))


AFDW_spreadist_upst <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'North Fork' | 
                                                                                  habdatdistinfo$River_Tributary == 'South Fork'),], 
                              aes(x = Spread_dist/1000, y = AFDWmean)) +
  geom_smooth(method='lm',aes(group=River_Tributary), color='darkgrey', alpha=0.3) + 
  geom_linerange(aes(ymin=DWmean-1.96*(DWsd/sqrt(macn)), ymax=DWmean+1.96*(DWsd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_point(size=3, alpha=0.75, color='#2b8cbe')+ 
  geom_vline(xintercept = confluences) +
  theme_classic() + 
  annotate('text', x=260,y=0.25, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  # annotation_custom(g_nf, xmin=105, xmax=155, ymin=0.5, ymax=3) +
  # annotation_custom(g_sf, xmin=40, xmax=90, ymin=0.5, ymax=3) +
  annotate('text', x=83, y=0.08, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=0.08, label='South Fork', angle=90, size=2.5) +
  geom_segment(data=Chelasegdat, aes(x=x1, xend=x2, y=y1, yend=y2))  +
  annotate('text', x=101, y=2.5, label='***', size=4)+
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-1,1), expand=c(0,0), breaks=c(0,0.20,0.40),name='Macroinvertebrate biomass (AFDW, g)') + 
  coord_cartesian(ylim=c(0,0.6)) + 
  geom_vline(xintercept = confluences) +
  theme(plot.margin=unit(c(0.25,0.5,0,0), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
  
#grid.draw(rbind(ggplotGrob(AFDW_spreadist_upst),ggplotGrob(AFDW_spreadist_main), size='last'))

  
######################## GREEN ALGAE AND DIATOMS ########################
green_spreadist_main <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'Lower mainstem' | 
                                                                                     habdatdistinfo$River_Tributary == 'Upper mainstem'),], 
                                 aes(x = Spread_dist/1000, y = greenmean)) +
  geom_linerange(aes(ymin=greenmean-1.96*(greensd/sqrt(macn)), ymax=greenmean+1.96*(greensd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_smooth(aes(group=River_Tributary), color='darkgrey', alpha=0.3, method='lm') + 
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=50, label='Downstream', angle=90, size=2.5, color='#016c59')+
  annotate('text', x=160, y=7, label='Mainstem', size=2.5) +
  #annotation_custom(g, xmin=200, xmax=250, ymin=0.5, ymax=3) +
  #geom_text(aes(label=Site)) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-3,5), breaks=c(0,0.25,0.50,0.75,1.00),name='Green algae concentration (Ch-a ug/cm-2)') + 
  coord_cartesian(ylim=c(0,1.00))+
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(),
        axis.line.x=element_blank(),
        plot.margin=unit(c(-0.1,0.5,0,0), "cm"))

green_spreadist_upst <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'North Fork' | 
                                                                                     habdatdistinfo$River_Tributary == 'South Fork'),], 
                                 aes(x = Spread_dist/1000, y = greenmean)) +
  geom_smooth(method='lm',aes(group=River_Tributary), color='darkgrey', alpha=0.3) + 
  geom_linerange(aes(ymin=greenmean-1.96*(greensd/sqrt(macn)), ymax=greenmean+1.96*(greensd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_point(size=3, alpha=0.75, color='#2b8cbe')+ 
  geom_vline(xintercept = confluences) +
  theme_classic() + 
  annotate('text', x=260,y=50, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  annotate('text', x=83, y=15, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=15, label='South Fork', angle=90, size=2.5) +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-0.50,2.00), breaks=c(0,.25,.50,.75,1.00),name='Green algae concentration (Ch-a ug/cm-2)') + 
  coord_cartesian(ylim=c(0,1.00))+
  geom_hline(yintercept=0) +
  theme(plot.margin=unit(c(0.25,0.5,0,0), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x=element_blank())


diatom_spreadist_main <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'Lower mainstem' | 
                                                                                   habdatdistinfo$River_Tributary == 'Upper mainstem'),], 
                               aes(x = Spread_dist/1000, y = diatommean)) +
  geom_linerange(aes(ymin=diatommean-1.96*(diatomsd/sqrt(macn)), ymax=diatommean+1.96*(diatomsd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_smooth(aes(group=River_Tributary), color='darkgrey', alpha=0.3, method='lm') + 
  geom_point(size=3, alpha=0.75, color='#1c9099')+ 
  geom_vline(xintercept = confluences) +
  annotate('text', x=260,y=200, label='Downstream', angle=90, size=2.5, color='#016c59')+
  annotate('text', x=160, y=300, label='Mainstem', size=2.5) +
  #annotation_custom(g, xmin=200, xmax=250, ymin=0.5, ymax=3) +
  #geom_text(aes(label=Site)) +
  theme_classic() + 
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-0.50,6.00), breaks=c(0,1.00,2.00,3.00,4.00),name='Diatom concentration (Ch-a ug/cm-2)') + 
  coord_cartesian(ylim=c(0,4.00))+
  geom_hline(yintercept = 0) +
  theme(axis.line.x=element_blank(),
        plot.margin=unit(c(-0.1,0.5,0,0), "cm"))

diatom_spreadist_upst <- ggplot(habdatdistinfo[!is.na(habdatdistinfo$AFDWmean) & (habdatdistinfo$River_Tributary == 'North Fork' | 
                                                                                   habdatdistinfo$River_Tributary == 'South Fork'),], 
                               aes(x = Spread_dist/1000, y = diatommean)) +
  geom_smooth(method='lm',aes(group=River_Tributary), color='darkgrey', alpha=0.3) + 
  geom_linerange(aes(ymin=diatommean-1.96*(diatomsd/sqrt(macn)), ymax=diatommean+1.96*(diatomsd/sqrt(macn))), alpha=0.75, size=1, color='#2b8cbe') +
  geom_point(size=3, alpha=0.75, color='#2b8cbe')+ 
  geom_vline(xintercept = confluences) +
  theme_classic() + 
  annotate('text', x=260,y=200, label='Upstream', angle=90, size=2.5,color='#045a8d') +
  # annotation_custom(g_nf, xmin=105, xmax=155, ymin=0.5, ymax=3) +
  # annotation_custom(g_sf, xmin=40, xmax=90, ymin=0.5, ymax=3) +
  annotate('text', x=83, y=50, label='North Fork', angle=90, size=2.5) +
  annotate('text', x=40, y=50, label='South Fork', angle=90, size=2.5) +
  scale_x_continuous(limits=c(0,270), expand=c(0,0), name='Distance from initial introduction (km)') +
  scale_y_continuous(limits=c(-0.50,6.00), breaks=c(0,1.00,2.00,3.00,4.00),name='Diatom concentration (Ch-a ug/cm-2)') + 
  coord_cartesian(ylim=c(0,4.00))+
  geom_vline(xintercept = confluences) +
  geom_hline(yintercept=0) +
  theme(plot.margin=unit(c(0.25,0.5,0,0), "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x=element_blank())


pdf(file.path(figdir,'Trends/ENV2.pdf'), width=7, height=5)
grid.arrange(rbind(ggplotGrob(degdays_spreadist_upst), ggplotGrob(degdays_spreadist_main),
                   ggplotGrob(AFDW_spreadist_upst), ggplotGrob(AFDW_spreadist_main), size='last'),
             rbind(ggplotGrob(green_spreadist_upst),ggplotGrob(green_spreadist_main),
                   ggplotGrob(diatom_spreadist_upst),ggplotGrob(diatom_spreadist_main),
                   size='last'),
             ncol=2)
dev.off()

#################################################################################################################################
####################### COMPARE DATA WITH SORENSON 2012 ###################
sorenson_082010<- read.csv(file.path(datadir,"Orusticus_netdata_Sorenson_082010.csv"), colClasses = c('factor', 'numeric', 'numeric','character','numeric','character'))
sorenson_082011<- read.csv(file.path(datadir,"Orusticus_netdata_Sorenson_082011.csv"), colClasses = c('factor', 'numeric', 'numeric','character','numeric','character'))
clydeholliday_2016 <- craydatOR[craydatOR$Site == 35,]
str(sorenson_082010)
colnames(clydeholliday_2016[,1:30])

sorenson_082010$Year <- 2010
sorenson_082011$Year <- 2011
sorenson <- rbind(sorenson_082010, sorenson_082011)
colnames(sorenson)[3] <- "Chelae_L" 
clydeholliday_2016$Year <- 2016

clyde20102016 <- rbind(sorenson[,c('Sex', 'CL', 'Chelae_L', 'Weight','Year')],clydeholliday_2016[,c('Sex', 'CL', 'Chelae_L', 'Weight','Year')])

ggplot(clyde20102016, aes(x=CL,y=Chelae_L, color=factor(Year))) + geom_point() + geom_smooth(span=1)
ggplot(clyde20102016, aes(x=log(CL),y=log(Chelae_L), color=factor(Year))) + geom_point() + geom_smooth()

ggplot(clyde20102016, aes(x=CL,y=Weight, color=factor(Year))) + geom_point() + geom_smooth(span=1)
ggplot(clyde20102016, aes(x=log(CL),y=log(Weight), color=factor(Year))) + geom_point() + geom_smooth(method='lm')

#Test differences in chela length
chela_CL_2010 <-lm(log(Chelae_L)~log(CL), data = clyde20102016[!is.na(clyde20102016$Chelae_L) & clyde20102016$Year==2010,])
summary(chela_CL_2010)
chela_CL_2011 <-lm(log(Chelae_L)~log(CL), data = clyde20102016[!is.na(clyde20102016$Chelae_L) & clyde20102016$Year==2011,])
summary(chela_CL_2011)
chela_CL_2016 <-lm(log(Chelae_L)~log(CL), data = clyde20102016[!is.na(clyde20102016$Chelae_L) & clyde20102016$Year==2016,])
summary(chela_CL_2016)

bla <- lmer(log(Chelae_L)~log(CL) + (log(CL)|Year), data=clyde20102016, REML=F)
summary(bla)

chela_CL <-lm(log(Chelae_L)~log(CL), data = clyde20102016[!is.na(clyde20102016$CL & !is.na(clyde20102016$Chelae_L)),])
summary(chela_CL)
clyde20102016[!is.na(clyde20102016$CL) & !is.na(clyde20102016$Chelae_L),'residchela'] <- residuals(chela_CL)
ggplot(clyde20102016, aes(x=factor(Year), y=residchela)) + geom_boxplot() +
  scale_y_continuous(name='Relative chela length', breaks=c(-0.25,0,0.25,0.5,0.75)) + 
  scale_x_discrete(name='Year') +
  theme_bw()

t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2010,'residchela'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2011,'residchela'])
t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2011,'residchela'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2016,'residchela'])
t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2010,'residchela'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2016,'residchela'])

#Test differences in weight
weight_CL_2010 <-lm(log(Weight)~log(CL), data = clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2010,])
summary(weight_CL_2010)
weight_CL_2011 <-lm(log(Weight)~log(CL), data = clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2011,])
summary(weight_CL_2011)
weight_CL_2016 <-lm(log(Weight)~log(CL), data = clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2016,])
summary(weight_CL_2016)

bla <- lmer(log(Weight)~log(CL) + (log(CL)|Year), data=clyde20102016, REML=F)
summary(bla)

weight_CL <-lm(log(Weight)~log(CL), data = clyde20102016[!is.na(clyde20102016$CL & !is.na(clyde20102016$Weight)),])
summary(weight_CL)
clyde20102016[!is.na(clyde20102016$CL) & !is.na(clyde20102016$Weight),'resid'] <- residuals(weight_CL)
ggplot(clyde20102016[!is.na(clyde20102016$CL) & !is.na(clyde20102016$Weight),], aes(x=resid)) + geom_histogram(aes(fill=factor(Year)))

ggplot(clyde20102016, aes(x=factor(Year), y=resid)) + geom_boxplot()+
  scale_y_continuous(name='Relative weight', limits=c(-0.4,0.75),breaks=c(-0.25,0,0.25,0.5,0.75)) + 
  scale_x_discrete(name='Year') +
  theme_bw()


t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2010,'resid'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2011,'resid'])
t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2011,'resid'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2016,'resid'])
t.test(clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2010,'resid'], 
       clyde20102016[!is.na(clyde20102016$Weight) & clyde20102016$Year==2016,'resid'])
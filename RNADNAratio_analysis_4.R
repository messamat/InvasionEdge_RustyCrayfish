library(readxl)
library(WriteXLS)
library(reshape2)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(abind)

rootdir <- "F:/Chapter3_ecoevo/" #UPDATE!
datadir <- file.path(rootdir, 'data/RNADNA_analysis/JDR_OR')
resdir <- file.path(rootdir, 'results/RNADNA_analysis')

############################FUNCTIONS TO IMPORT, FORMAT, AND PROCESS DATA ##############################
plotstandard <- function(data, standard, params, plateref, ylim) {
  ggplot(data[grepl(paste(standard), data$inst),], aes(x=con, y=rep1)) + 
    geom_point() + 
    geom_point(aes(y=rep2)) + 
    #stat_smooth(method='lm') +
    geom_abline(intercept=params[1], slope=params[2]) +
    annotate(geom='text', label=paste(plateref, standard, sep=' '), x=600, y = 5) + 
    scale_y_continuous(limits=c(0,ylim))+
    theme_bw()
} 

RNADNAcomp <- function(date, plate) {
  print(paste0('Analyzing RNA/DNA ratio data from ', date, ', plate ', plate))
  
  ref <- paste(date, '0', plate,sep="")
  #Read in tables after making sure they exist. Otherwise, throw error
  root <- paste(file.path(datadir,'Mathis_2017'),date,'_RNADNAratio_Plate',as.character(plate),sep="")
  if (file.exists(paste(root,'_format','.xls',sep=''))) {
    fluo <- read_xls(paste(root,'_format','.xls',sep=''))
  } else {
    if (file.exists(paste(root,'.xls',sep=''))) {
      fluo <- read_xls(paste(root,'.xls',sep=''))
    } else {
      stop(paste('The file ',root," doesn't exist.",sep=''))
    }
  }
  
  #Read layout
  layout <- read_xls(paste(root,'_layout','.xls',sep=''))
  #Make sure that non-NA values are all unique. Otherwise, throw error
  if (length(unique(unlist(layout)[which(!is.na(layout))]))/length(unlist(layout)[which(!is.na(layout))]) != 1) {
    stop(paste('The file ',root," layout table has duplicate values",sep=''))
  }
  
  #Format data
  setDT(fluo)
  setDT(layout)
  fluo_melt <- melt(fluo, measure.vars=2:13)
  layout_melt <- melt(layout, measure.vars=2:13)
  fluo_format <- data.frame(sample=layout_melt$value,plate1=fluo_melt$value)
  fluo_format$plate1 <- as.numeric(as.character(fluo_format$plate1))
  fluo_format$sample <- as.character(fluo_format$sample)
  fluo_format <- fluo_format[!is.na(fluo_format$sample),]
  fluo_format[grepl("RD", fluo_format$sample)|
                  grepl("RR",fluo_format$sample)|
                  grepl("DD",fluo_format$sample)|
                  grepl("DR",fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("RD", fluo_format$sample)|
                                                                                      grepl("RR",fluo_format$sample)|
                                                                                      grepl("DD",fluo_format$sample)|
                                                                                      grepl("DR",fluo_format$sample), 'sample'], 1,3)
  fluo_format$rep <- paste('rep', substr(fluo_format$sample, (nchar(fluo_format$sample)),nchar(fluo_format$sample)),sep="")
  fluo_format[grepl("OR", fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 1,8)
  fluo_format[grepl("OR", fluo_format$sample), 'tube'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 10,10)
  
  #Format curves
  curves_sub <- fluo_format[!grepl("OR", fluo_format$sample),]
  curves_formatcast <- dcast(curves_sub[,2:4], inst ~ rep, value.var = 'plate1')
  curves_formatcast$repmean <- rowMeans(curves_formatcast[,c("rep1", "rep2")])
  refcon <- data.frame(inst=c('DD1','DD2','DD3','DD4','DD5','DD6','DR1','DR2','DR3','DR4','DR5','DR6',
                                             'RD1','RD2','RD3','RD4','RD5','RD6','RR1','RR2','RR3','RR4','RR5','RR6'),
                                      con=c(0,40,80,160,400,800,0,40,80,160,400,800,
                                            0,50,100,200,500,1000,0,50,100,200,500,1000))
  curvesformat <- merge(curves_formatcast, refcon, by='inst',all.y=F)
  
  #Compute regression coefficients
  lmDD <- lm(repmean~con, data=curvesformat[grepl("DD", curvesformat$inst),])
  #print(summary(lmDD))
  lmDR <- lm(repmean~con, data=curvesformat[grepl("DR", curvesformat$inst),])
  #print(summary(lmDR))
  lmRD <- lm(repmean~con, data=curvesformat[grepl("RD", curvesformat$inst),])
  #print(summary(lmRD))
  lmRR <- lm(repmean~con, data=curvesformat[grepl("RR", curvesformat$inst),])
  #print(summary(lmRR))
  
  mDD <- coef(lmDD)[[2]]
  mDR <- coef(lmDR)[[2]]
  mRD <- coef(lmRD)[[2]]
  mRR <- coef(lmRR)[[2]]
  
  delta <- mDD/mDR
  rho <- mRR/mRD
  
  summary(lmDD)$adj.r.squared
  
  params <- data.frame(ref, 
              mRD,mRD_R2=summary(lmRD)$adj.r.squared,
              mRR,mRR_R2=summary(lmRR)$adj.r.squared,
              mDD,mDD_R2=summary(lmDD)$adj.r.squared, 
              mDR,mDR_R2=summary(lmDR)$adj.r.squared,
              delta,rho)
  
  #Plot standards
  plotRD <- plotstandard(curvesformat, "RD", coef(lmRD), plateref=ref, ylim=4000)
  plotRR <- plotstandard(curvesformat, "RR", coef(lmRR), plateref=ref, ylim=150)
  plotDR <- plotstandard(curvesformat, "DR", coef(lmDR), plateref=ref, ylim=6000)
  plotDD <- plotstandard(curvesformat, "DD", coef(lmDD), plateref=ref, ylim=150)
  
  plotcurves <- grid.arrange(plotRD,plotRR,plotDR,plotDD, ncol=2)
  #grid.draw(plotcurves)
  
  #Format samples fluorescence
  fluo_formatsamples <- fluo_format[!is.na(fluo_format$tube),]
  fluo_formatsamples$Site <- as.numeric(substr(fluo_formatsamples$inst, 3,5))
  fluo_formatsamples$Cray_ID <- as.numeric(substr(fluo_formatsamples$inst, 7,8))
  
  fluo_formatcast <- dcast(fluo_formatsamples, Site+Cray_ID+tube ~ rep, value.var = 'plate1')
  fluo_formatcast$repmean <- rowMeans(fluo_formatcast[,c("rep1", "rep2")], na.rm=T)
  fluo_formatcast <- dcast(fluo_formatcast, Site+Cray_ID~tube, value.var='repmean')
  datformat <- fluo_formatcast
  datformat$ref <- ref
  
  ##Substract fluorescence of tube C from that of tube A
  FaFc <- with(datformat, A-C)
  #Calculate the total RNA concentration
  datformat$RNAcon <- (FaFc/(1-rho))/mRD
  #Subtract the fluorescence of the tube C from that of tube B
  FbFc <- with(datformat, B-C)
  #Calculate the total DNA concentration
  datformat$DNAcon <- (FbFc/(1-delta))/mDR
  
  #Calculate RNA/DNA ratio
  datformat$RNADNAratio <- with(datformat, RNAcon/DNAcon)
  
  return(list(plots=plotcurves, data=datformat, params=params))
}

#Compute average RD/DR ratio for plates for which RNA standard was still good
mean(RNADNA_params_format$RD_DRratio[1:13])
#Function that uses RNA standard
RNADNAcomp_RNAreplace <- function(date, plate) {
  print(paste0('Analyzing RNA/DNA ratio data from ', date, ', plate ', plate))
  
  ref <- paste(date, '0', plate,sep="")
  #Read in tables after making sure they exist. Otherwise, throw error
  root <- paste(file.path(datadir,'Mathis_2017'),date,'_RNADNAratio_Plate',as.character(plate),sep="")
  if (file.exists(paste(root,'_format','.xls',sep=''))) {
    fluo <- read_xls(paste(root,'_format','.xls',sep=''))
  } else {
    if (file.exists(paste(root,'.xls',sep=''))) {
      fluo <- read_xls(paste(root,'.xls',sep=''))
    } else {
      stop(paste('The file ',root," doesn't exist.",sep=''))
    }
  }
  
  #Read layout
  layout <- read_xls(paste(root,'_layout','.xls',sep=''))
  #Make sure that non-NA values are all unique. Otherwise, throw error
  if (length(unique(unlist(layout)[which(!is.na(layout))]))/length(unlist(layout)[which(!is.na(layout))]) != 1) {
    stop(paste('The file ',root," layout table has duplicate values",sep=''))
  }
  
  #Format data
  setDT(fluo)
  setDT(layout)
  fluo_melt <- melt(fluo, measure.vars=2:13)
  layout_melt <- melt(layout, measure.vars=2:13)
  fluo_format <- data.frame(sample=layout_melt$value,plate1=fluo_melt$value)
  fluo_format$plate1 <- as.numeric(as.character(fluo_format$plate1))
  fluo_format$sample <- as.character(fluo_format$sample)
  fluo_format <- fluo_format[!is.na(fluo_format$sample),]
  fluo_format[grepl("RD", fluo_format$sample)|
                grepl("RR",fluo_format$sample)|
                grepl("DD",fluo_format$sample)|
                grepl("DR",fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("RD", fluo_format$sample)|
                                                                                grepl("RR",fluo_format$sample)|
                                                                                grepl("DD",fluo_format$sample)|
                                                                                grepl("DR",fluo_format$sample), 'sample'], 1,3)
  fluo_format$rep <- paste('rep', substr(fluo_format$sample, (nchar(fluo_format$sample)),nchar(fluo_format$sample)),sep="")
  fluo_format[grepl("OR", fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 1,8)
  fluo_format[grepl("OR", fluo_format$sample), 'tube'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 10,10)
  
  #Format curves
  curves_sub <- fluo_format[!grepl("OR", fluo_format$sample),]
  curves_formatcast <- dcast(curves_sub[,2:4], inst ~ rep, value.var = 'plate1')
  curves_formatcast$repmean <- rowMeans(curves_formatcast[,c("rep1", "rep2")])
  refcon <- data.frame(inst=c('DD1','DD2','DD3','DD4','DD5','DD6','DR1','DR2','DR3','DR4','DR5','DR6',
                              'RD1','RD2','RD3','RD4','RD5','RD6','RR1','RR2','RR3','RR4','RR5','RR6'),
                       con=c(0,40,80,160,400,800,0,40,80,160,400,800,
                             0,50,100,200,500,1000,0,50,100,200,500,1000))
  curvesformat <- merge(curves_formatcast, refcon, by='inst',all.y=F)
  
  #Compute regression coefficients
  lmDD <- lm(repmean~con, data=curvesformat[grepl("DD", curvesformat$inst),])
  #print(summary(lmDD))
  lmDR <- lm(repmean~con, data=curvesformat[grepl("DR", curvesformat$inst),])
  #print(summary(lmDR))
  lmRD <- lm(repmean~con, data=curvesformat[grepl("RD", curvesformat$inst),])
  #print(summary(lmRD))
  lmRR <- lm(repmean~con, data=curvesformat[grepl("RR", curvesformat$inst),])
  #print(summary(lmRR))
  
  mDD <- coef(lmDD)[[2]]
  mDR <- coef(lmDR)[[2]]
  #Adjust RNA with DNA standard curve and RNA/DNA ratio from good standard curves to include variability from ribogreen concentration and potential photobleaching
  mRD <- 0.51*coef(lmDR)[[2]]
  mRR <- coef(lmRR)[[2]]
  
  delta <- mDD/mDR
  #Keep original Rho to reflect completeness of RNase digestion
  rho <- mRR/coef(lmRD)[[2]]
  
  summary(lmDD)$adj.r.squared
  
  params <- data.frame(ref, 
                       mRD,mRD_R2=NA,
                       mRR,mRR_R2=summary(lmRR)$adj.r.squared,
                       mDD,mDD_R2=summary(lmDD)$adj.r.squared, 
                       mDR,mDR_R2=summary(lmDR)$adj.r.squared,
                       delta,rho)
  
  #Plot standards
  plotRD <- plotstandard(curvesformat, "RD", coef(lmRD), plateref=ref, ylim=4000)
  plotRR <- plotstandard(curvesformat, "RR", coef(lmRR), plateref=ref, ylim=150)
  plotDR <- plotstandard(curvesformat, "DR", coef(lmDR), plateref=ref, ylim=6000)
  plotDD <- plotstandard(curvesformat, "DD", coef(lmDD), plateref=ref, ylim=150)
  
  plotcurves <- grid.arrange(plotRD,plotRR,plotDR,plotDD, ncol=2)
  #grid.draw(plotcurves)
  
  #Format samples fluorescence
  fluo_formatsamples <- fluo_format[!is.na(fluo_format$tube),]
  fluo_formatsamples$Site <- as.numeric(substr(fluo_formatsamples$inst, 3,5))
  fluo_formatsamples$Cray_ID <- as.numeric(substr(fluo_formatsamples$inst, 7,8))
  
  fluo_formatcast <- dcast(fluo_formatsamples, Site+Cray_ID+tube ~ rep, value.var = 'plate1')
  fluo_formatcast$repmean <- rowMeans(fluo_formatcast[,c("rep1", "rep2")], na.rm=T)
  fluo_formatcast <- dcast(fluo_formatcast, Site+Cray_ID~tube, value.var='repmean')
  datformat <- fluo_formatcast
  datformat$ref <- ref
  
  ##Substract fluorescence of tube C from that of tube A
  FaFc <- with(datformat, A-C)
  #Calculate the total RNA concentration
  datformat$RNAcon <- (FaFc/(1-rho))/mRD
  #Subtract the fluorescence of the tube C from that of tube B
  FbFc <- with(datformat, B-C)
  #Calculate the total DNA concentration
  datformat$DNAcon <- (FbFc/(1-delta))/mDR
  
  #Calculate RNA/DNA ratio
  datformat$RNADNAratio <- with(datformat, RNAcon/DNAcon)
  
  return(list(plots=plotcurves, data=datformat, params=params))
}

##############################RUN ON ALL PLATES ####################
RNADNA_090101 <- RNADNAcomp('0901',1)
RNADNA_dat <- as.data.frame(setNames(replicate(length(colnames(RNADNA_090101$data)),numeric(0), simplify = F),colnames(RNADNA_090101$data)))
RNADNA_params <-  as.data.frame(setNames(replicate(length(colnames(RNADNA_090101$params)),numeric(0), simplify = F),colnames(RNADNA_090101$params)))

labdates <- unique(substr(list.files(datadir), 12,15))

#Run function for plates with the right RNA standard curves
for (date_loop in labdates[1:5]) {
  for (plate_loop in 1:4) {
    tryCatch({
      #print(paste(file.path(datadir,'Mathis_2017'),date_loop,'_RNADNAratio_Plate',as.character(plate_loop),'.xls',sep=""))
      RNADNA_comped <- RNADNAcomp(date_loop,plate_loop)
      RNADNA_dat <- rbind(RNADNA_dat,RNADNA_comped$data)
      RNADNA_params <- rbind(RNADNA_params,RNADNA_comped$params)
      ggsave(filename= paste0('curves',date_loop,'0',plate_loop,'.png'), plot=RNADNA_comped$plots,path= file.path(resdir, 'QA_QC'))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

#Run function for plates with an adjusted RNA standard cruve
for (date_loop in labdates[6:12]) {
  for (plate_loop in 1:4) {
    tryCatch({
      #print(paste(file.path(datadir,'Mathis_2017'),date_loop,'_RNADNAratio_Plate',as.character(plate_loop),'.xls',sep=""))
      RNADNA_comped <- RNADNAcomp_RNAreplace(date_loop,plate_loop)
      RNADNA_dat <- rbind(RNADNA_dat,RNADNA_comped$data)
      RNADNA_params <- rbind(RNADNA_params,RNADNA_comped$params)
      ggsave(filename= paste0('curves',date_loop,'0',plate_loop,'.png'), plot=RNADNA_comped$plots,path=file.path(resdir, 'QA_QC'))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

############################# Compute 091404 that has no standard curves #########################################################
date='0914'
plate=4
print(paste0('Analyzing RNA/DNA ratio data from ', date, ', plate ', plate))

ref <- paste(date, '0', plate,sep="")
#Read in tables after making sure they exist. Otherwise, throw error
root <- paste(file.path(datadir,'Mathis_2017'),date,'_RNADNAratio_Plate',as.character(plate),sep="")
if (file.exists(paste(root,'_format','.xls',sep=''))) {
  fluo <- read_xls(paste(root,'_format','.xls',sep=''))
} else {
  if (file.exists(paste(root,'.xls',sep=''))) {
    fluo <- read_xls(paste(root,'.xls',sep=''))
  } else {
    stop(paste('The file ',root," doesn't exist.",sep=''))
  }
}

#Read layout
layout <- read_xls(paste(root,'_layout','.xls',sep=''))
#Make sure that non-NA values are all unique. Otherwise, throw error
if (length(unique(unlist(layout)[which(!is.na(layout))]))/length(unlist(layout)[which(!is.na(layout))]) != 1) {
  stop(paste('The file ',root," layout table has duplicate values",sep=''))
}

#Format data
setDT(fluo)
setDT(layout)
fluo_melt <- melt(fluo, measure.vars=2:13)
layout_melt <- melt(layout, measure.vars=2:13)
fluo_format <- data.frame(sample=layout_melt$value,plate1=fluo_melt$value)
fluo_format$plate1 <- as.numeric(as.character(fluo_format$plate1))
fluo_format$sample <- as.character(fluo_format$sample)
fluo_format <- fluo_format[!is.na(fluo_format$sample),]
fluo_format[grepl("RD", fluo_format$sample)|
              grepl("RR",fluo_format$sample)|
              grepl("DD",fluo_format$sample)|
              grepl("DR",fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("RD", fluo_format$sample)|
                                                                              grepl("RR",fluo_format$sample)|
                                                                              grepl("DD",fluo_format$sample)|
                                                                              grepl("DR",fluo_format$sample), 'sample'], 1,3)
fluo_format$rep <- paste('rep', substr(fluo_format$sample, (nchar(fluo_format$sample)),nchar(fluo_format$sample)),sep="")
fluo_format[grepl("OR", fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 1,8)
fluo_format[grepl("OR", fluo_format$sample), 'tube'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 10,10)

#Get regression coefficients from all plates done that day and compute average
RNADNA_0914params <- RNADNA_params[RNADNA_params$ref %in% c('091401','091402','091403'),]

mDD <- mean(RNADNA_0914params$mDD)
mDR <- mean(RNADNA_0914params$mDR)
mRD <- mean(RNADNA_0914params$mRD)
mRR <- mean(RNADNA_0914params$mRR)

delta <- mDD/mDR
rho <- mRR/mRD

params <- data.frame(ref, 
                     mRD,mRD_R2=NA,
                     mRR,mRR_R2=NA,
                     mDD,mDD_R2=NA, 
                     mDR,mDR_R2=NA,
                     delta,rho)

#Format samples fluorescence
fluo_formatsamples <- fluo_format[!is.na(fluo_format$tube),]
fluo_formatsamples$Site <- as.numeric(substr(fluo_formatsamples$inst, 3,5))
fluo_formatsamples$Cray_ID <- as.numeric(substr(fluo_formatsamples$inst, 7,8))

fluo_formatcast <- dcast(fluo_formatsamples, Site+Cray_ID+tube ~ rep, value.var = 'plate1')
fluo_formatcast$repmean <- rowMeans(fluo_formatcast[,c("rep1", "rep2")], na.rm=T)
fluo_formatcast <- dcast(fluo_formatcast, Site+Cray_ID~tube, value.var='repmean')
datformat <- fluo_formatcast
datformat$ref <- ref

##Substract fluorescence of tube C from that of tube A
FaFc <- with(datformat, A-C)
#Calculate the total RNA concentration
datformat$RNAcon <- (FaFc/(1-rho))/mRD
#Subtract the fluorescence of the tube C from that of tube B
FbFc <- with(datformat, B-C)
#Calculate the total DNA concentration
datformat$DNAcon <- (FbFc/(1-delta))/mDR

#Calculate RNA/DNA ratio
datformat$ratio <- with(datformat, RNAcon/DNAcon)

#Combine results with main dataset
RNADNA_dat <- rbind(RNADNA_dat,datformat)
RNADNA_params <- rbind(RNADNA_params,params)

############################# Compute 092604 that has no standard curves aside from DR#########################################################
#First made sure that the DNA curve was very similar to otherplates that day
date='0926'
plate=4
print(paste0('Analyzing RNA/DNA ratio data from ', date, ', plate ', plate))

ref <- paste(date, '0', plate,sep="")
#Read in tables after making sure they exist. Otherwise, throw error
root <- paste(file.path(datadir,'Mathis_2017'),date,'_RNADNAratio_Plate',as.character(plate),sep="")
if (file.exists(paste(root,'_format','.xls',sep=''))) {
  fluo <- read_xls(paste(root,'_format','.xls',sep=''))
} else {
  if (file.exists(paste(root,'.xls',sep=''))) {
    fluo <- read_xls(paste(root,'.xls',sep=''))
  } else {
    stop(paste('The file ',root," doesn't exist.",sep=''))
  }
}

#Read layout
layout <- read_xls(paste(root,'_layout','.xls',sep=''))
#Make sure that non-NA values are all unique. Otherwise, throw error
if (length(unique(unlist(layout)[which(!is.na(layout))]))/length(unlist(layout)[which(!is.na(layout))]) != 1) {
  stop(paste('The file ',root," layout table has duplicate values",sep=''))
}

#Format data
setDT(fluo)
setDT(layout)
fluo_melt <- melt(fluo, measure.vars=2:13)
layout_melt <- melt(layout, measure.vars=2:13)
fluo_format <- data.frame(sample=layout_melt$value,plate1=fluo_melt$value)
fluo_format$plate1 <- as.numeric(as.character(fluo_format$plate1))
fluo_format$sample <- as.character(fluo_format$sample)
fluo_format <- fluo_format[!is.na(fluo_format$sample),]
fluo_format[grepl("RD", fluo_format$sample)|
              grepl("RR",fluo_format$sample)|
              grepl("DD",fluo_format$sample)|
              grepl("DR",fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("RD", fluo_format$sample)|
                                                                              grepl("RR",fluo_format$sample)|
                                                                              grepl("DD",fluo_format$sample)|
                                                                              grepl("DR",fluo_format$sample), 'sample'], 1,3)
fluo_format$rep <- paste('rep', substr(fluo_format$sample, (nchar(fluo_format$sample)),nchar(fluo_format$sample)),sep="")
fluo_format[grepl("OR", fluo_format$sample), 'inst'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 1,8)
fluo_format[grepl("OR", fluo_format$sample), 'tube'] <- substr(fluo_format[grepl("OR", fluo_format$sample), 'sample'], 10,10)

#Get regression coefficients from all plates done that day and compute average
RNADNA_0926params <- RNADNA_params_format[RNADNA_params$ref %in% c('092601','092602','092603'),]

mDD <- mean(RNADNA_0926params$mDD)
mDR <- mean(RNADNA_0926params$mDR)
mRD <- mean(RNADNA_0926params$mRD)
mRR <- mean(RNADNA_0926params$mRR)

delta <- mDD/mDR
rho <- mRR/mRD

params <- data.frame(ref, 
                     mRD,mRD_R2=NA,
                     mRR,mRR_R2=NA,
                     mDD,mDD_R2=NA, 
                     mDR,mDR_R2=NA,
                     delta,rho)

#Format samples fluorescence
fluo_formatsamples <- fluo_format[!is.na(fluo_format$tube),]
fluo_formatsamples$Site <- as.numeric(substr(fluo_formatsamples$inst, 3,5))
fluo_formatsamples$Cray_ID <- as.numeric(substr(fluo_formatsamples$inst, 7,8))

fluo_formatcast <- dcast(fluo_formatsamples, Site+Cray_ID+tube ~ rep, value.var = 'plate1')
fluo_formatcast$repmean <- rowMeans(fluo_formatcast[,c("rep1", "rep2")], na.rm=T)
fluo_formatcast <- dcast(fluo_formatcast, Site+Cray_ID~tube, value.var='repmean')
datformat <- fluo_formatcast
datformat$ref <- ref

##Substract fluorescence of tube C from that of tube A
FaFc <- with(datformat, A-C)
#Calculate the total RNA concentration
datformat$RNAcon <- (FaFc/(1-rho))/mRD
#Subtract the fluorescence of the tube C from that of tube B
FbFc <- with(datformat, B-C)
#Calculate the total DNA concentration
datformat$DNAcon <- (FbFc/(1-delta))/mDR

#Calculate RNA/DNA ratio
datformat$RNADNAratio <- with(datformat, RNAcon/DNAcon)

#Combine results with main dataset
RNADNA_dat <- rbind(RNADNA_dat,datformat)
RNADNA_params <- rbind(RNADNA_params,params)

#######################################################################################################################################
#######Look at distribution of parameters#########
qplot(RNADNA_params$mRD)
qplot(RNADNA_params$mRR)
qplot(RNADNA_params$mDD)
qplot(RNADNA_params$mDR)

params_relative<- apply(RNADNA_params[,c(2,4,6,8)], MARGIN=2, FUN=function(x) {x/mean(x)})
colnames(params_relative) <- paste(colnames(params_relative),'_rel',sep='')
RNADNA_params_format <- cbind(RNADNA_params, params_relative)
RNADNA_params_format <- RNADNA_params_format[,c(1,order(colnames(RNADNA_params_format[,(2:ncol(RNADNA_params_format))]))+1)]
RNADNA_params_format$RD_DRratio <- RNADNA_params_format$mRD/RNADNA_params_format$mDR
RNADNA_params_format$RD_RRratio <- RNADNA_params_format$mRD/RNADNA_params_format$mRR

###################################################### INSPECT AND CORRECT STANDARD CURVES #######################################################################
##090102 Corrections: replace 500ng/ml and 1000ng/ml values as higher than other plates that day, which increases RD/DR ratio 
# RNADNA_090101 <- read_xls(paste0(datadir, 'Mathis_20170901_RNADNAratio_Plate1_format.xls'))
# RNADNA_090102 <- read_xls(paste0(datadir, 'Mathis_20170901_RNADNAratio_Plate2.xls'))
# RNADNA_090103 <- read_xls(paste0(datadir, 'Mathis_20170901_RNADNAratio_Plate3.xls'))
# RNADNA_090102[1,12:13] <- rowMeans(abind(RNADNA_090101[1,12:13], RNADNA_090103[1,12:13], along=3), dims=2, na.rm=T)
# WriteXLS(RNADNA_090102, paste0(datadir,'Mathis_20170901_RNADNAratio_Plate2_format.xls'),row.names = FALSE)

# #090402 & 090403 Corrections: replace the first five points of the DR curve with the average of values from a day for which the max value for DR is similar
# RNADNA_090402 <- read_xls(paste0(datadir, 'Mathis_20170904_RNADNAratio_Plate2.xls'))
# RNADNA_090403 <- read_xls(paste0(datadir, 'Mathis_20170904_RNADNAratio_Plate3.xls'))
# 
# RNADNA_092401 <- read_xls(paste0(datadir, 'Mathis_20170924_RNADNAratio_Plate1.xls'))
# RNADNA_092402 <- read_xls(paste0(datadir, 'Mathis_20170924_RNADNAratio_Plate2.xls'))
# RNADNA_092403 <- read_xls(paste0(datadir, 'Mathis_20170924_RNADNAratio_Plate3.xls'))
# RNADNA_092404 <- read_xls(paste0(datadir, 'Mathis_20170924_RNADNAratio_Plate4.xls'))
# 
# RNADNA_090402[4,2:11] <- rowMeans(abind(RNADNA_092401[4,2:11],RNADNA_092402[4,2:11],RNADNA_092403[4,2:11],RNADNA_092404[4,2:11], along=3), dims=2, na.rm=T)
# RNADNA_090403[4,2:11] <- rowMeans(abind(RNADNA_092401[4,2:11],RNADNA_092402[4,2:11],RNADNA_092403[4,2:11],RNADNA_092404[4,2:11], along=3), dims=2, na.rm=T)
# 
# WriteXLS(RNADNA_090402, paste0(datadir,'Mathis_20170904_RNADNAratio_Plate2_format.xls'),row.names = FALSE)
# WriteXLS(RNADNA_090403, paste0(datadir,'Mathis_20170904_RNADNAratio_Plate3_format.xls'),row.names = FALSE)
# 
# #090601, 090602, 090603 Corrections: DD curve has an errroneous value for 400 ng/ml
# RNADNA_090601 <- read_xls(paste0(datadir, 'Mathis_20170906_RNADNAratio_Plate1.xls'))
# RNADNA_090602 <- read_xls(paste0(datadir, 'Mathis_20170906_RNADNAratio_Plate2.xls'))
# RNADNA_090603 <- read_xls(paste0(datadir, 'Mathis_20170906_RNADNAratio_Plate3.xls'))
# 
# RNADNA_090601[3,10:11] <- NA
# RNADNA_090602[3,10:11] <- NA
# RNADNA_090603[3,10:11] <- NA
# 
# WriteXLS(RNADNA_090601, paste0(datadir,'Mathis_20170906_RNADNAratio_Plate1_format.xls'),row.names = FALSE)
# WriteXLS(RNADNA_090602, paste0(datadir,'Mathis_20170906_RNADNAratio_Plate2_format.xls'),row.names = FALSE)
# WriteXLS(RNADNA_090603, paste0(datadir,'Mathis_20170906_RNADNAratio_Plate3_format.xls'),row.names = FALSE)
# 
# #091402 Corrections: max RR concentration, nothing was in the plate
# RNADNA_091402 <- read_xls(paste0(datadir, 'Mathis_20170914_RNADNAratio_Plate2.xls'))
# RNADNA_091402[2,12] <- NA
# 
# WriteXLS(RNADNA_091402, paste0(datadir,'Mathis_20170914_RNADNAratio_Plate2_format.xls'),row.names = FALSE)
# 
# ############################################# INSPECT AND CORRECT CRAYFISH RNA/DNA DATA #######################################################
# RNADNA_090702 <- read_xls(paste0(datadir, 'Mathis_20170907_RNADNAratio_Plate2.xls'))
# RNADNA_090702[7,2:7] <- NA
# WriteXLS(RNADNA_090702, paste0(datadir,'Mathis_20170907_RNADNAratio_Plate2_format.xls'),row.names = FALSE)
# 
# RNADNA_092404 <- read_xls(paste0(datadir, 'Mathis_20170924_RNADNAratio_Plate4.xls'))
# RNADNA_092404[6,2:7] <- NA
# WriteXLS(RNADNA_092404, paste0(datadir,'Mathis_20170924_RNADNAratio_Plate4_format.xls'),row.names = FALSE)

RNADNA_dat$day <- substr(RNADNA_dat$ref,1,4)
ggplot(RNADNA_dat, aes(x=ref, y=ratio, color=factor(Site))) + geom_point(size=3)
ggplot(RNADNA_dat, aes(x=Site, y=ratio, color=factor(ref))) + geom_point(size=3)
ggplot(RNADNA_dat, aes(x=Site, y=ratio, color=factor(day))) + geom_point(size=3)

RNADNA_dat[duplicated(RNADNA_dat[,c('Site','Cray_ID')]),]

craydat <- read.csv(file.path(rootdir, 'data/Fieldworkdata/Crayfish_4.csv'))
craydat_RNADNA <- merge(craydat, RNADNA_dat, by=c('Site',"Cray_ID"), all.x=T, all.y=T)
#2 missing ones 27/27, and 2/6

write.csv(RNADNA_dat, file.path(resdir,"RNADNA_datformat.csv"), row.names=F)
write.csv(RNADNA_params_format, file.path(resdir,"RNADNA_paramsformat.csv"), row.names=F)

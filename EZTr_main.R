##################
# EZTr_main #

# This framework is developed as part of the research article:
# Automated Discretization of ‘Transpiration Restriction to Increasing VPD’ Features from Outdoors High-Throughput Phenotyping Data
#
# Soumyashree Kar, Ryokei Tanaka, Lijalem Balcha Korbu, Jana Kholová, Hiroyoshi Iwata, Surya S. Durbha, J. Adinarayana, Vincent Vadez
##################

#' HTP data processing
#'
#' Convert  gravimetric sensors or load cells (LC) time series into
#' Evapotranspiration (ETr) and Transpiration (Tr) time series.
#'
#' The function iterates over the different time points (at 15 min interval) of
#' each sector or pot to generate a matrix of Tr values. Fifteen biologically
#' relevant features are extracted from the Tr time series for each day,
#' followed by computation of feature heritability.
#'
#' The entire process comprises four steps: 1st - Conversion of LC data to ETr;
#' 2nd - Calculation of reference i.e. Penman Monteith ET for the given weather
#' condition, and filtering ETr based on reference ET; 3rd - Generating smooth
#' ETr time series; 4th - Extraction of Tr from smooth ETr, Tr features and
#' each feature's broad-sense heritability estimate.
#' 
#' Input Files required:
#' 1.Loadcells²
#' 2.Loadcells metadata
#' 3.Sensor climate data
#' 4.Sensor unit map
#' 5.Plant eye
#'
#' #' Functions used for each step:
#' Step1 - extractRawLCmatrix(), curateRawLCgetET(), filterETrExtremeCols();
#' Step2 - extractWthrVar(), prepcsWthr(), calculateETref(),generateETrRatio(),
#'         genThreshVal(), dataPART(), threshETr(), ordFiltETr();
#' Step3 - smoothETr();
#' Step4 - calculateTr(), getFeatures(), getFeatureHe().
#'
#'
#' The function includes seven inputs, 'data' includes LC data,
#' experimental design and genotype-replicate information in metadata, weather
#' data and genotype-replicate specific leaf area data; first (Date1) and
#' last (Date2) dates of the experiment; 'irrg.dts' is a vector of dates when
#' plants were irrigated; 'opPATH' and opPATH.smth to store intermediate
#' feature-specific files.
#'
#' @param HTP_data /code{list} of 5 dataframes composed of input files
#' a) loadcell data with unit, genotype, g_alias, treatment, timestamp
#' and Mass(g) columns
#' b) metadata with unit, old_unit, Experiment, Treatment, Species
#' Genotype, G. Alias and Replicates columns
#' c) weather data with sensor, variable, timestamp, value columns
#' d) solar radiation data is stored in separate dataframe in format c)
#' e) leaf area data with Sector, Experiment, Treatment, Species,
#' Genotype, G..Alias, Replicates, timestamp and leaf area columns
#'
#' @param lastDate /code{string} containing last date of experiment
#' in YYYY-MM-DD format.
#'
#' @param irrg.dts /code{vector} of irrigation dates, each date mentioned as
#' string in YYYY-MM-DD format. If no such date is needed, same as the last date.
#'
#' @param Date1 /code{string} value specifying first date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param Date2 /code{string} value specifying last date of time series in
#' YYYY-MM-DD hh:mm:ss format.
#'
#' @param opPATH /code{string} value specifying the directory path for
#' intermediate results to be stored.
#'
#' @param opPATH.smth /code{string} value specifying the directory path for
#' feature-specific results to be stored.
#'
#' @return Return:
#'
#' The function returns a list of fifteen outcomes:
#' LCraw_tsmeta = Matrix of time stamp values present in raw load cell data,
#' LCraw = Matrix of raw load cell data,
#' ETrmm_Obs = Matrix of raw or observed ETr in mm
#' ETr_Obs_ERR.SEC.nms = Matrix of sector names with very high proportion of extreme values,
#' ETmm_Obs_FINAL = Matrix of ETr in mm after removing the erroneous sectors,
#' Wthr_agg15min = Matrix of all weather variables aggregated for 15 minutes interval,
#' Wthr.ETref.ETobs = Matrix of combined time series of Reference ET, ETr, weather values,
#' ETrRatio_TW1_ThreshVALS = Matrix of filtered values of Day-time ref ET/ETr ratios,
#' ETrRatio_TW2_ThreshVALS = Matrix of filtered values of Night-time ref ET/ETr ratios,
#' IrrgFilt_ETr = Matrix of Irrigation filtered ETr matrix,
#' IrrgFilt_ETr_Imptd = Matrix of irrigation filtered ETr matrix imputed,
#' ETr_smth = Matrix of smooth ETr time series,
#' Tr = Matrix of Tr time series,
#' featureH2 = Matrix of each feature heritability estimate on each day,
#' eachFeature_TS = List of each feature's time series for all genotypes.
#'
#' @author Soumyashree Kar email<ksoumya2301@gmail.com>
#'
#' @references
#'
#' Vadez, V., Kholová, J., Hummel, G., Zhokhavets, U., Gupta, S. K., &
#' Hash, C. T. (2015). LeasyScan: a novel concept combining 3D imaging and
#' lysimetry for high-throughput phenotyping of traits controlling plant water
#' budget. Journal of Experimental Botany, 66(18), 5581-5593.


####### Stage 0 : Libraries, packages, definition of parameters  ####### 

## Load libraries
library("easypackages")
libraries("readxl", "hms", "xts", "dplyr","mgcv","PerformanceAnalytics", "wavelets",  
          "signal", "tidyverse", "zoo", "h2o", "sqldf", "ggplot2", "plyr",
          "lubridate", "BioFTF", "plantecophys", "highfrequency", "stringr",
          "chron", "nonlinearTseries", "tsfeatures", "splitstackshape", "psych", "fixr", "dplyr")

library(SpATS)
library(reshape2)       
library("tidyr")
library(tidyr)     # data management
library(nlme)      # nonlinear model
library(locfit)    # local regression
library(gss)       # smoothing splines
library(factoextra)       
library(gridExtra)       
library(mgcv)
library(statgenHTP)
library(magick)
library(segmented)
library(MASS)
library(tidyverse)
library(rstatix)
library(stringr)
library(statgenSTA)
library(berryFunctions )


## Working directory
setwd(dir = "C:/Users/2021lg003/Documents/UPDATE_KAR_PIPELINE_LEASYSCAN_EZTr/")

# Load all the functions #
# Functions needed for processing Stage-I: LC to ETr extraction
source('./functions/extractRawLCmatrix.R') # Included 'Treatment col in meta.d'
source('./functions/curateRawLC.R')
source('./functions/filterLCExtremeCols.R')
source('./functions/getETr.R')
source('./functions/convETr.R')
source('./functions/filterETrExtremeCols.R')

# Functions needed for processing Stage-II: ETref extraction and ETr thresholding
source('./functions/extractWthrVar.R')
source('./functions/prepcsWthr.R')
source('./functions/calculateETref.R')
source('./functions/generateETrRatio.R')
source('./functions/genThreshVal.R')
source('./functions/dataPART.R')
source('./functions/threshETr.R')
source('./functions/ordFiltETr.R')

# Functions needed for processing Stage-III: raw and smooth Tr, Tr-features and feature-H2 extraction
# source('./functions/calculateTr.R')
source('./functions/getFeatures.R')
source('./functions/getFeatureHe.R')
source('./functions/smoothETr.R')
colMax <- function(data) sapply(data, max, na.rm = TRUE) ##function to find max values on the entire dataset
colMin <- function(data) sapply(data, min, na.rm = TRUE) ##function to find min values on the entire dataset


### Parameters need to be defined : 

## Last date
ldt <- readline(prompt = "Enter LAST DATE of experiment (YYYY-MM-DD): ")
lastDate=as.Date(as.character(ldt)) # follow this format  2023-10-18

## First date
# NOTE 1: Enter the day BEFORE the start of the experiment. Otherwise the next steps will "cut" the first day of measurement 
fdt <- readline(prompt = "Enter FIRST DATE of experiment (YYYY-MM-DD): ")
firstDate=as.Date(as.character(fdt)) # format 2023-09-27

## Time interval between two measurements
seq_by <- readline(prompt = "Enter desired time (in min) interval, e.g. 15/30/45/60: ")
seq_by=as.numeric(seq_by)

## Expected weight on the LC in grams 
avg_wgt <- readline(prompt = "Enter the expected weight (in grams) on the LC: ")
avg_wgt=as.numeric(avg_wgt) ## format 70000

## Irrigation dates (format required : "2023-09-30")
# NOTE 1: Enter date(s) of irrigation or when data is extremely noisy.
# NOTE 2: If no such date, enter LAST DATE of experiment
# NOTE 3: The entire date will be deleted : in the case of evening irrigation, it is not necessary to delete the entire date.
irrg.dts <- c("2023-09-30","2023-10-03","2023-10-06", "2023-10-10", "2023-10-13", "2023-10-17", "2023-10-18")

## Date 1 and 2 : format "YYYY-MM-DD HH:MM:SS"
Date1="2023-09-26 23:46:00" # Day before 'firstDate'(For 15min, "2021-04-19 23:46:00")
Date2="2023-10-18 23:45:00" #'lastDate'(For 15min, "2021-04-29 23:45:00")

### Define path to store results
opPATH.obj= "./saved_objects/"
opPATH <- "./results/"
opPATH.smth="./results/smthFeaturesTimeSeries/"
opPATH.raw="./results/rawFeaturesTimeSeries/"
opPATH.profile="./results/profile_example/"
opPATH.spats="./results/STAGENHTP/"

#### Prepare the LC data and meta data
# Load data
load("./data/Exp_60_LEASYSCAN_LAURA.RData")

allData <- IRD_data_Exp_NEW

# Get load cells data i.e. weights of sector
m.lc <- allData$m.lc

# Get Genotype and Exp design metadata
meta.d <- allData$meta.d
meta.d <- distinct(meta.d) #remove duplicate entries
meta.d <- meta.d[,1:8]

# Get the list of species from the metadata
species.nm <- unique(meta.d$Species)

# Include the species ID e.g. 1 for 'Sorghum' as per the data
meta.d.sp <- meta.d[meta.d$Species==species.nm[1], ]

# Find sectors with missing metadata
noEntrySecs <- which(!unique(m.lc$unit) %in% unique(meta.d.sp$unit))
noEntrySecNms <- unique(m.lc$unit)[noEntrySecs]

# Remove sectors-without-metadata from original loadcells file
m.lc <- m.lc[! m.lc$unit %in% noEntrySecNms, ]

####### Stage-I: Process LC data and generate ETr matrix #######
st.time <- Sys.time()

# Extract matrix of loadcell data #
LC.MAT.OP <- extractRawLCmatrix(x = m.lc, y = meta.d.sp, 
                                s=firstDate, z = lastDate, 
                                inter = seq_by)
# save(LC.MAT.OP, file = paste(opPATH.obj, "LC.MAT.OP.RData", sep = ""))
load( file = paste(opPATH.obj, "LC.MAT.OP.RData", sep = ""))

LC.MAT.f <- LC.MAT.OP$LC.MAT.f
LC.MAT.TSinfo <- LC.MAT.OP$LC_tsmeta

# Outliers detection
check_for_negative_values(LC.MAT.f) ### LC.MAT.f can contain negative weights
max_df=as.data.frame(colMax(LC.MAT.f)) ### LC.MAT.f can contain outliers weights 
max(max_df[6:nrow(max_df),]) #is this a possible weight?

# Draw the weight profile of the first LC to follow the data cleaning process (0-LC.MAT.f).
data=(LC.MAT.f) 
pdf(file = paste(opPATH.profile,"1-WEIGHT_RAW.pdf", sep = ""), width = 4,height = 4) 
plot<-
  ggplot(data =data, aes(data$TS, data[,2]))+
  ylab("Weight (g)")+
  xlab("Timestamp")+
  labs(title = paste("RAW DATA",colnames(data)[2]))+
  geom_point()
print(plot)
dev.off()


# Add metadata to LC matrix : 
# select from Metadata: "unit", "old_unit", 'Trtmt', "Genotype, "G..Alias", "Replicates"
meta.d.LCmat <- meta.d.sp[meta.d.sp$unit %in% colnames(LC.MAT.f)[-1], 
                          c(1,2,4,6,7,8)] 

LC.MAT.f.t <- as.data.frame(t(LC.MAT.f))

colnames(LC.MAT.f.t) <- LC.MAT.f$TS
# save(TS, file = paste(opPATH.obj, "TS.RData", sep = ""))

LC.MAT.f.t <- LC.MAT.f.t[-1,]

# reorder rows of 'meta.d.sp' according to rownames/unit of LC.MAT.f.t
meta.LCDF <- meta.d.LCmat[order(match(meta.d.LCmat$unit, rownames(LC.MAT.f.t))), ]

LC.MAT.raw <- as.data.frame(cbind(meta.LCDF, LC.MAT.f.t))
save(LC.MAT.raw, file = paste(opPATH.obj, "LC.MAT.raw.RData", sep = ""))
write.table(LC.MAT.raw, paste0(opPATH, "OP-1","_LCraw_wNA.csv"), sep = ";", row.names = F)


# Start outlier detection, removal and imputation of LC Matric to generate ETr profiles #
## STEP 1. Replace outliers by NA. Outliers if : negative weight + weight out the +/- 30% average weight + outlier defined boxplot of each LC,
## STEP 2. Linear imputation ,
## STEP 3. Identify remain columns with maximum missing, replace by  NA
## STEP 4. Keep if 30% of values

imputed.DF.final <- curateRawLC(x = LC.MAT.f, y = meta.LCDF, z= avg_wgt) ## +/- 30 % DONE
# save(imputed.DF.final, file = paste(opPATH.obj, "imputed.DF.final.RData", sep = ""))
load(file = paste(opPATH.obj, "imputed.DF.final.RData", sep = ""))

# Outliers detection post curateRawLC
check_for_negative_values(imputed.DF.final) ##imputed.DF.final should not contain negative values
max_df=as.data.frame(colMax(imputed.DF.final[7:ncol(imputed.DF.final)]))
max(max_df[1:nrow(max_df),])  ##max of imputed.DF.final should not be an outlier value

write.table(imputed.DF.final, paste0(opPATH,"OP-2","_LC_olrm_imputed.csv"), sep = ";", row.names = F)


# Identify the highly extreme valued sectors #
err.sec.info <- filterLCExtremeCols(x = imputed.DF.final, y = meta.LCDF)

err.sec.nm <- err.sec.info$err.sec.NM

err.sec.meta <- err.sec.info$err.sec.META

write.table(err.sec.meta, paste0(opPATH, "OP-3","_LCimp_errorUnits.csv"), sep = ";", row.names = F)


# Remove the err.cols i.e. sectors with extreme values
impData.errSEC.rmvd <- imputed.DF.final[!imputed.DF.final$unit %in% err.sec.nm, ]
# save(impData.errSEC.rmvd, file = paste(opPATH.obj, "impData.errSEC.rmvd.RData", sep = ""))

# Now, the weight profile should be cleaned. Verify on the first LC 
data=t(impData.errSEC.rmvd)
colnames(data)=data[1,]
data=data[-c(1:6),]
data<- as.data.frame(apply(data, 2, as.numeric))
data$TS= colnames(impData.errSEC.rmvd)[7:ncol(impData.errSEC.rmvd)]
data$TS =ymd_hms(data$TS)

pdf(file = paste(opPATH.profile,"2-WEIGHT_FILTERED.pdf", sep = ""), width = 4,height = 4)
plot<-
  ggplot(data =data, aes(TS, data[,1]))+
  ylab("Weight (g)")+
  xlab("Timestamp")+
  labs(title = paste("WEIGHT AFTER FILTERING",colnames(data)[1]))+
  geom_point()
print(plot)
dev.off()

# Generate ETr profiles from "impData.errSEC.rmvd" dataframe #
et.vals <- getETr(x = impData.errSEC.rmvd)

save(et.vals, file = paste(opPATH.obj, "et.vals.RData", sep = ""))
# load(file = paste(opPATH.obj, "et.vals.RData", sep = ""))

et.obs <- et.vals$obsETr_core

ETr_Meta <- et.vals$obsETr_meta


# Convert ETr in grams to mm (Y/N) # ## !! careful to put the Y or N in CAPITAL LETTER 
ETr_F <- convETr(x = ETr_Meta, y = et.obs) ## AJOUTER UNE PARTIE ADAPTEE pour L'IRD (surface LC differente avec ICRISAT)
# 
save(ETr_F , file = paste(opPATH.obj,  "ETr_F.RData", sep = ""))
# load(file = paste(opPATH.obj, "ETr_F.RData", sep = ""))


# #  Draw the profile of a LC  (weight/ETr) to follow the data cleaning process (2-ETr_F).
data=t(ETr_F)
colnames(data)=data[1,]
data=data[-c(1:6),]
TS=rownames(data)
data_num <- as.data.frame(apply(data, 2, as.numeric))
data=as.data.frame(data)
data=data_num
data$TS= TS
data$TS =ymd_hms(data$TS )

pdf(file = paste(opPATH.profile,"3-ET_RAW_VALUES.pdf", sep = ""), width = 4,height = 4)
plot<-
  ggplot(data =data, aes(TS, data[,1]))+
  ylab("ET (mm.15min-1)")+
  xlab("Timestamp")+
  labs(title = paste("RAW ET DATA",colnames(data)[1]))+
  geom_point()
print(plot)
dev.off()

write.table(ETr_F, paste0(opPATH, "OP-5","_ETr_Obs.csv"), col.names = TRUE, row.names = FALSE, sep = ";", dec = ".")


####### Stage-II: Process Weather data to obtain ETref + filter ETr based on ETref ####### 
wthr.DFagg15min <- IRD_data_Exp_NEW$climate
wthr.DFagg15min$Date <- lubridate::dmy(wthr.DFagg15min$Date) #Format here
wthr.DFagg15min <- wthr.DFagg15min[wthr.DFagg15min$Date>=firstDate & 
                                     wthr.DFagg15min$Date<=lastDate,]

wthr.DFagg15min$TS <- ymd_hms(paste0(as.character(wthr.DFagg15min$Date), 
                                     wthr.DFagg15min$TS))
wthr.DFagg15min <- wthr.DFagg15min[,-1]


wthr.DFagg15min= wthr.DFagg15min[-c(2061: nrow(wthr.DFagg15min)),]  ## there is empty row at the end of dataset .. why? need to be removed

# Create empty base matrix for weather variables
end.seq_hhmm  <- '23:45'

if (seq_by == '15'){
  end.seq_hhmm  <- '23:45'
} else if(seq_by == '30'){
  end.seq_hhmm  <- '23:30'
} else if(seq_by == '45'){
  end.seq_hhmm  <- '23:15'
} else if(seq_by == '60'){
  end.seq_hhmm  <- '23:00'}

TS_base<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
                                        ymd_hm(paste0(lastDate," ",end.seq_hhmm)), 
                                        by = paste0(seq_by, " ", "mins"))))
#  
# TS_base<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
#                                         ymd_hm(paste0(lastDate," ",'23:45')), by = '15 mins')))

names(TS_base)[1]<-c("int.val")

TS_base$time <- strftime(TS_base$int.val, format="%H:%M:%S", tz="UTC")

# Since, the hms values slightly differ in the orginal dataset than the ideal 15min interval values,
# replace them with the TS_base strftime (hms) values

hms.ts.base<-unique(TS_base$time)
names(hms.ts.base)<-c("time")

print("Climate matrix timestamp mapping status")

i<-nrow(wthr.DFagg15min)
pbar <- create_progress_bar('text')
pbar$init(i)

for(i in 1:nrow(wthr.DFagg15min))
{
  if(! wthr.DFagg15min$TS[i] %in% hms.ts.base)
  {
    j <- which.min(abs(chron(times=hms.ts.base) - 
                         chron(times=format(wthr.DFagg15min$TS[i],format = "%H:%M:%S")))) # find which value in ts-base-vector is nearest to each-DP
    
    # assign the nearest ts-base-vector value to that DP
    dd <- date(wthr.DFagg15min$TS[i])
    wthr.DFagg15min$TS[i] <- ymd_hms(paste0(dd, TS_base$time[j]))} 
  pbar$step()
}

# TS_ALL<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
#                                        ymd_hm(paste0(lastDate," ",'23:45')), by = '15 mins')))
TS_ALL<-as.data.frame(as.character(seq(ymd_hm(paste0(firstDate," ",'00:00')),
                                       ymd_hm(paste0(lastDate," ",end.seq_hhmm)), 
                                       by = paste0(seq_by, " mins"))))

names(TS_ALL)[1]<-c("TS.n") # MUST be the same as in original data set m.lc.df

TS_ALL$TS.n <- ymd_hms(TS_ALL$TS.n)

# create empty dataframe to store all values
Wthr.MAT <- as.data.frame(matrix(nrow = length(TS_ALL$TS.n), 
                                 ncol = ncol(wthr.DFagg15min)))

Wthr.MAT[ ,1] <- TS_ALL$TS.n; names(Wthr.MAT)[1] <- "TS"

colnames(Wthr.MAT)[2:ncol(Wthr.MAT)] <- names(wthr.DFagg15min)[2:ncol(wthr.DFagg15min)]

df <- merge(x = wthr.DFagg15min, y = Wthr.MAT, by = "TS", all = TRUE) # perform outer join to merge by id=TS.n
df.new <- df[,1:6]; names(df.new)[2:6]<-names(wthr.DFagg15min)[2:6]
df.new <- df.new[!duplicated(df.new[c("TS")]),]
# colnames(df.new)[2] <- names(wthr.DFagg15min)[1]


etrDTS <- as.data.frame(colnames(ETr_F)[7:ncol(ETr_F)])
names(etrDTS) <- "TS"
etrDTS$TS <- ymd_hms(etrDTS$TS)

df.new <- df.new[df.new$TS %in% etrDTS$TS, ] # subset after outer join to ensure that NAs don't add extra rows

# Pre-process weather
df.new$Temp <- if(sum(is.na(df.new$Temp))>0){df.new$Temp<-na.aggregate.default(df.new$Temp)}
df.new$RH <- if(sum(is.na(df.new$RH))>0){df.new$RH<-na.aggregate.default(df.new$RH)}
df.new$SR <- if(sum(is.na(df.new$SR))>0){df.new$SR<-na.aggregate.default(df.new$SR)}
df.new$WS <- if(sum(is.na(df.new$WS))>0){df.new$WS<-na.aggregate.default(df.new$WS)}
wthr.DFagg15min.filt <- df.new


# Compute VPD and insert into the weather DF #
SVP <- 610.7*(10^(7.5*wthr.DFagg15min.filt[ ,2]/(237.3+wthr.DFagg15min.filt[ ,2])))
VPD <- ((1 - (wthr.DFagg15min.filt[ ,3]/100))*SVP)/1000
wthr.DFagg15min.filt[ ,4] <- VPD

et.obs <- ETr_F
save(et.obs, file = paste(opPATH.obj,  "et.obs.RData", sep = ""))


# Calculate Penman Monteith ET #
wthr.df1 <- calculateETref(x=wthr.DFagg15min.filt)
max(wthr.df1$ETref) ## possible or outlier?
wthr.ETref.df <- as.data.frame(wthr.df1)
empty.MAT <- matrix(nrow = 8, 
                    ncol = (ncol(et.obs)-nrow(wthr.df1)))

# select columns "Temp"  "RH"    "VPD"   "SR"    "WS"    "Tmax"  "Tmin"  "ETref"
empty.MAT.wthr.ETref <- as.data.frame(cbind(empty.MAT, t(wthr.ETref.df[,c(2:6, 9:11)])))
colnames(empty.MAT.wthr.ETref) <- colnames(et.obs)
wthr.ETref.ETobs <- as.data.frame(rbind(empty.MAT.wthr.ETref, et.obs))
save(wthr.ETref.ETobs, file = paste(opPATH.obj,  "wthr.ETref.ETobs.RData", sep = ""))

# Remove irrigation dates #
file.colnms <- colnames(wthr.ETref.ETobs)
wthr.ETref.ETobs <- wthr.ETref.ETobs[ ,!substr(file.colnms,1,10) %in% irrg.dts]
save(wthr.ETref.ETobs, file = paste(opPATH.obj, "wthr.ETref.ETobs.RData", sep = "")) 

## save some parts of wthr.ETref.ETobs dataframe to reuse after (the format of wthr.ETref.ETobs will be use for the features extraction)
TS= colnames(wthr.ETref.ETobs[7:ncol(wthr.ETref.ETobs)])
save(TS,file= paste(opPATH.obj, "TS.RData", sep = "") )
LC= wthr.ETref.ETobs[9:nrow(wthr.ETref.ETobs),1]
save(LC,file= paste(opPATH.obj, "LC.RData", sep = "") )
metad_emptyrows= wthr.ETref.ETobs[,1:6]
weather= wthr.ETref.ETobs[1:8,]
save(weather, file =  paste(opPATH.obj,"weather.RData", sep = ""))


## Filter ETobs based on ETREF values on solar radiation values (SR)  instead of hours (in the initial version) 
### if SR = 0 -->  night period --> ET = 0
### if SR > 0 -->  day period --> -0.01<ET<ETref

ET_ratio_mat <- generateETrRatio(x = wthr.ETref.ETobs) 

save(ET_ratio_mat, file =  paste(opPATH.obj,"ET_ratio_mat.RData", sep = ""))

ET_ratio_mat[, 9:ncol(ET_ratio_mat) ]= round(ET_ratio_mat[, 9:ncol(ET_ratio_mat) ] , 3)
rownames(ET_ratio_mat)= TS
colnames(ET_ratio_mat)[9:ncol(ET_ratio_mat)]= LC


#  come back to wthr.ETref.ETobs format:  with metadata + a dataframe (TS in col, LC in row)
wthr.ETref.ETobs_ratio= t(ET_ratio_mat)
ETr_smth_metad= cbind(metad_emptyrows, wthr.ETref.ETobs_ratio) ## ET filtered non normalized, with the expected format (wthr.ETref.ETobs format)

save(ETr_smth_metad, file =  paste(opPATH.obj,"ETr_smth_metad.RData", sep = ""))
write.table(ETr_smth_metad, paste0(opPATH, "ETr_smth_metad.csv"), sep = ";", row.names = F, dec = ".")

# Outliers detection
check_for_negative_values(ETr_smth_metad[(9:nrow(ETr_smth_metad)),(7:ncol(ETr_smth_metad))]) 
max_df=as.data.frame(colMax(ETr_smth_metad[(9:nrow(ETr_smth_metad)),(7:ncol(ETr_smth_metad))]))  
max(max_df[1:nrow(max_df),]) 


# # Identify error plots from ETr values using the similar method as above #
ETr_err.sec.info <- filterETrExtremeCols(x = ETr_smth_metad, y = meta.LCDF) ## FILTRER APRES LA PARTIE ET ETREF !!!

err.sec.nm <- ETr_err.sec.info$ETr_err.sec.NM

err.sec.meta <- ETr_err.sec.info$ETr_err.sec.META

if(length(err.sec.nm)>0){write.table(err.sec.meta,
                                     paste0(opPATH, "OP-6","_ET_Obs_ERR.SEC.nms.csv"), sep = ";", row.names = F, dec = ".")}

# Remove the err.cols i.e. sectors with extreme values #
ETr_Meta_ERRsec.rmvd <- ETr_smth_metad

ETr_Meta_ERRsec.rmvd <- ETr_Meta_ERRsec.rmvd[!ETr_Meta_ERRsec.rmvd$unit %in% err.sec.nm, ]

write.table(ETr_Meta_ERRsec.rmvd, paste0(opPATH,
                                         "OP-7","_ETr_Obs_FINAL.csv"), sep = ";", row.names = F, dec = ".")
save(ETr_Meta_ERRsec.rmvd, file =  paste(opPATH.obj,"ETr_Meta_ERRsec.rmvd.RData", sep = ""))
load(paste(opPATH.obj,"ETr_Meta_ERRsec.rmvd.RData", sep = ""))


# Outliers detection
max_df=as.data.frame(colMax(ETr_Meta_ERRsec.rmvd[(9:nrow(ETr_Meta_ERRsec.rmvd)),(7:ncol(ETr_Meta_ERRsec.rmvd))]))  
max(max_df[1:nrow(max_df),]) 

## END of filtering process


# Now the ET profile should be cleaned. Draw the ET profile of the first LC
data=t(ETr_Meta_ERRsec.rmvd)
data=data[,-c(1:8)]
data=data[-c(1:6),]
data_num <- as.data.frame(apply(data, 2, as.numeric))
data=data_num
data= cbind(TS, data)
data$TS =ymd_hms(data$TS )
colnames(data)[2:ncol(data)]=LC

pdf(file = paste(opPATH.profile, "4-ET_PROFIL_FILTERED.pdf", sep = ""), width = 4,height = 4)
plot<-
  ggplot(data =data, aes(TS, data[,2]))+
  ylab("ET mm.15min-1")+
  xlab("Timestamp")+
  labs(title = paste("ET profile filtered",colnames(data)[1]))+
  geom_point()+
  geom_line()+
  ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), date_breaks = "1 days")+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
print(plot)
dev.off()



####### STAGE III: Process Plant Eye data with spatial correction ####### 
library("easypackages")

libraries("dplyr","LoadCellDataProcessing", "statgenHTP", "platformDataAnalysis" , 
          "ggplot2", "lubridate", "SpATS", "reshape2", "tidyr", "nlme", "locfit",
          "gss","factoextra", "gridExtra","mgcv", "magick", "segmented", "MASS",
          "tidyverse","rstatix")


##Input pe data and spatial correction
pe_data <- read.csv(file = "./data/Exp60 Sorghum Ref set IRD Trial Sep 2023-477 mm_20231027_planteye.csv", sep=";")
## !! note that pe_data have been prepared with position line and 
colnames(pe_data) ## 1] "barcode","n","r", "old_unit", "unit",  "Line" ,  "Plot", "ID","genotype", "timestamp","timeNumber",  "Digital_biomass",...  

##SI BESOIN DE RETOURNER au design 
# LABEL <- read.csv("C:/Users/2021lg003/Documents/LEASYSCAN_EXP_INDE_2023/DATASET/FULL_DATASET/Label_SPATS.csv", header=T, sep=";", dec = ".")
# LABEL=na.omit(LABEL)
# LABEL= LABEL[-3] ## removed Variety ID
# Label_geno <- read.csv("C:/Users/2021lg003/Documents/LEASYSCAN_EXP_INDE_2023/DATASET/FULL_DATASET/Label_Genotype.csv", header=T, sep=";")
#save for later
d_exp= read.csv2( "./data/d_exp_template.csv")
save(d_exp,file = paste(opPATH.obj,"d_exp.RData", sep = ""))


pe_data$genotype<-as.factor(pe_data$genotype) ##321 genotypes
pe_data$treatment<-as.factor(pe_data$treatment)
pe_data=na.omit(pe_data)
pe_data$timestamp= dmy_hm(pe_data$timestamp)
pe_data$date= date(pe_data$timestamp)
pe_data<-pe_data %>% arrange(timestamp) ##trie par ordre croissant

start=which((pe_data$date) == "2023-09-28")
end=which((pe_data$date) == "2023-10-17") 
pe_data=pe_data[(start[1]):(end[length(end)]),] 

pe_data$timeNumber= as.character(pe_data$timeNumber)
unique(pe_data$timeNumber)
for (i in  1:length(unique(pe_data$timeNumber))) {
  pe_data$timeNumber[pe_data$timeNumber== unique(pe_data$timeNumber)[i]] <- i
  
}
pe_data$timeNumber = as.factor(pe_data$timeNumber)

for (i in levels(pe_data$timeNumber)){
  Data  <- pe_data %>%
    dplyr::filter(pe_data$timeNumber == i) 
  DF= Data[duplicated(Data$ID), ] ##
  Data$timestamp = Data$timestamp[1] ## ne mettre qu'un seul timestamp commun pour tout le numero de scan
  Data_unique= Data[!duplicated(Data$ID), ] ##
  
  write.table(DF,file = "./DOUBLONS.csv", append = TRUE, col.names = FALSE, row.names = FALSE, sep = ";")
  write.table(Data_unique,file = "./pe_data_unique.csv", append = TRUE, col.names = FALSE, row.names = FALSE, sep = ";")
  
}

pe_data_unique= read.csv2("./pe_data_unique.csv", header = FALSE)
colnames(pe_data_unique)=colnames(pe_data)
pe_data_unique$timeNumber=as.factor(pe_data_unique$timeNumber)
pe_data=pe_data_unique
rm(pe_data_unique)
rm(Data_unique)



colnames(pe_data) <- mdf_raw_pe_colnames(colnames = colnames(pe_data))

ref_trait_nm <- c("Digital_biomass", "Height", "Leaf_angle", "Leaf_area",
                  "Leaf_area_index", "Leaf_area_projected", "Leaf_inclination",
                  "Light_penetration_depth")

pe_data <- pe_data %>% select(timeNumber, timestamp, unit, genotype, treatment,
                              Digital_biomass, Height, Leaf_angle, Leaf_area, Leaf_area_index,
                              Leaf_area_projected, Leaf_inclination, Light_penetration_depth)

pe_data[, 6:ncol(pe_data)] <- sapply(pe_data[, 6:ncol(pe_data)],  as.numeric)
pe_data$timestamp= lubridate::ymd_hms(pe_data$timestamp)


colnames(pe_data) <-  c("timeNumber", "timePoint", "unit", "genotype", "treatment","Digital_biomass", "Height", "Leaf_angle", "Leaf_area",
                        "Leaf_area_index", "Leaf_area_projected", "Leaf_inclination",
                        "Light_penetration_depth")
# add extra columns from the experimental design
pe_data <- add_exp_des_col(data = pe_data, d_exp_des = d_exp,
                           data_unit = "unit",
                           d_exp_unit = "new_unit",
                           col_add = c("rowNum", "colNum", "block", "plotId"))


# arrange the columns in a certain order
pe_data <- pe_data %>% select(timeNumber, timePoint, block, rowNum, colNum, plotId,
                              genotype, Digital_biomass, Height,
                              Leaf_angle, Leaf_area, Leaf_area_index,
                              Leaf_area_projected, Leaf_inclination,
                              Light_penetration_depth) 


# PE pipeline: remove tp with too high missing values
pe_data$Leaf_area= as.numeric(pe_data$Leaf_area) 
pe_data$timeNumber=as.numeric(pe_data$timeNumber)
# pdf(file = paste(opPATH.spats, "Leaf_area_raw_values.pdf", sep = ""), width = 4,height = 4) 
# plot =  ggplot(data =pe_data, aes(timeNumber, Leaf_area, colour = plotId))+
#   ylab("Leaf area")+
#   xlab("TimeNumber")+
#   geom_point()+
#   geom_line()
# 
# dev.off()
#  
### no missing values, aleready removed !
# prop_non_miss <- timepoint_prop_non_missing(pe_data)
# tp_rem <- names(prop_non_miss[prop_non_miss < 0.3])
# if(length(tp_rem) > 0){
#   pe_data <- pe_data[!(pe_data$timePoint %in% tp_rem), ]
# } 



# PE pipeline: remove outliers
# 
# pe_data <- outlier_boxplot_detect(pe_data)
# 
# plot_trend(data = pe_data, trait = "Height")
plot_trend(data = pe_data, trait = "Leaf_area")
# plot_trend(data = pe_data, trait = "Digital_biomass")

# PE pipeline: trim the time series

# remove days progressively the end of the exp until obtain a linear regression 

sel_tp <- sort(unique(pe_data$timePoint)) # all timePoint in the TS
sel_tp<-sel_tp[order(sel_tp)]
sel_tp

sel_tp <- sel_tp[-c(29:length(sel_tp))] # remove last days one by one, We decide to keep only the linear part at the start of the exp, until the "2023-10-11 ", where we see a inflexion point

pe_data <- pe_data[pe_data$timePoint %in% sel_tp, ]

save(pe_data,file = paste(opPATH.obj,"pe_data.RData", sep = ""))
load(file = paste(opPATH.obj,"pe_data.RData", sep = ""))

plot_trend(data = pe_data, trait = "Leaf_area") ## visualize to keep only the linear part at the start


# PE pipeline: Creation of TP object

TP_PE <- createTimePoints(dat = pe_data,
                          experimentName = "EXP60",
                          genotype = "genotype",
                          timePoint = "timePoint",
                          plotId = "plotId",
                          rowNum = "rowNum", colNum = "colNum")
summary(TP_PE)

save(TP_PE, file=paste(opPATH.obj,"TP_PE.RData", sep = ""))
load(file=paste(opPATH.obj,"TP_PE.RData", sep = ""))

# PE pipeline: Spatial adjustment

Timepoint <- getTimePoints(TP_PE)
save(Timepoint, file =paste(opPATH.obj, "Timepoint.RData", sep = ""))

## Plot the layout for the each point.
for (i in 1:nrow(Timepoint)) {
  png(file = paste(opPATH.spats,i, "_PE_LAYOUT.png", sep = "")) 
  plot<-
    plot(TP_PE, 
         plotType = "layout",
         timePoints = i,  
         traits = "Leaf_area")
  print(plot)
  dev.off()
  
}

modPhenoSpGDP <- NULL

modPhenoSpGDP <- fitModels(TP = TP_PE,
                           trait = "Leaf_area",
                           timePoints = c(1:nrow(Timepoint)) )

save(modPhenoSpGDP, file = paste(opPATH.obj,"modPhenoSpGDP.RData", sep = ""))
load(file = paste(opPATH.obj,"modPhenoSpGDP.RData", sep = ""))

summary(modPhenoSpGDP)

### plot SPATS results over timePoint
for (i in 1:nrow(Timepoint)) {
  png(file = paste(opPATH.spats,i, "_SPATS.png", sep = "")) 
  plot(modPhenoSpGDP,
       timePoints =  i,
       plotType = "spatial",
       spaTrend = "percentage",
       colorBy = "repId")
  dev.off()
}

## Plot heritability over time
png(file = paste(opPATH.spats,i, "_H2_over_time.png", sep = "")) 
plot(modPhenoSpGDP,
     plotType = "herit",
     ylim = c(0.5,1))
getHerit(modPhenoSpGDP)
dev.off()

## Plot residual, column and raw variance over time
png(file = paste(opPATH.spats,i, "_variance_over_time.png", sep = "")) 
plot(modPhenoSpGDP,
     plotType = "variance")
dev.off()

## Extract the corrected values:
spatCorrSpGDP <- getCorrected(modPhenoSpGDP, timePoints = 1:nrow(Timepoint))
write.table(spatCorrSpGDP, paste(opPATH.spats,"Leaf_Area_Corrected_statgenHTP.csv", sep = ""))


save(spatCorrSpGDP, file = paste(opPATH.obj,"spatCorrSpGDP.RData", sep = "" ))
load(file = paste(opPATH.obj,"spatCorrSpGDP.RData", sep = "" ))

spatCorrSpGDP$genotype= as.factor(spatCorrSpGDP$genotype)
spatCorrSpGDP$timePoint= as.factor(spatCorrSpGDP$timePoint)

## plot 3D-LA raw and 3D-LA corrected over time for each LC 
# for (i in levels(spatCorrSpGDP$plotId)) {
i=levels(spatCorrSpGDP$plotId)[1]
exemple_i <- spatCorrSpGDP %>% 
  dplyr::filter(spatCorrSpGDP$plotId== i) 
png(file = paste(opPATH.profile,i,"_LA_corr_raw.png", sep = ""))

plot<-
  ggplot(data =exemple_i, aes(timeNumber))+
  ylab("LA mm²")+
  xlab("Timestamp")+
  labs(title = paste("3D-LA and 3D-LA corr"))+
  geom_point(aes(y = Leaf_area_corr), color = "steelblue")+
  geom_point(aes(y = Leaf_area), color = "darkred")+
  geom_line(aes(y = Leaf_area_corr), color="steelblue", linetype="twodash") +
  geom_line(aes(y = Leaf_area), color="darkred", linetype="twodash") 
print(plot)
dev.off()

print(i)
# }



LB_MEAN_plot <- spatCorrSpGDP %>%
  select(plotId, Leaf_area_corr) %>%
  group_by(plotId) %>%
  summarize_all(funs(mean_Leaf_area_corr = mean (Leaf_area_corr, na.rm=T), n = n(), sd = sd(Leaf_area_corr, na.rm=T), se = sd(.)/sqrt(n-1))) %>%
  filter(is.finite(mean_Leaf_area_corr)) %>%
  filter(mean_Leaf_area_corr > 0)


LB_MEAN_genotype <- spatCorrSpGDP %>%
  select(genotype, Leaf_area_corr) %>%
  group_by(genotype) %>%
  summarize_all(funs(mean_Leaf_area_corr = mean (Leaf_area_corr, na.rm=T), n = n(), sd = sd(Leaf_area_corr, na.rm=T), se = sd(.)/sqrt(n-1))) %>%
  filter(is.finite(mean_Leaf_area_corr)) %>%
  filter(mean_Leaf_area_corr > 0)

save(LB_MEAN_genotype, file = paste(opPATH.obj,"LB_MEAN_genotype.RData", sep = ""))
save(LB_MEAN_plot, file = paste(opPATH.obj,"LB_MEAN_plot.RData", sep = ""))
# 
# load(file = "LB_MEAN_genotype.RData")
# load(file = "LB_MEAN_plot.RData")
# 
# 
# LB_MEAN_plot_day<- spatCorrSpGDP %>%
#   select(plotId, Leaf_area_corr, date) %>%
#   group_by(date, plotId) %>%
#   summarize_all(funs(mean_Leaf_area_corr_day = mean (Leaf_area_corr, na.rm=T))) %>%
#   filter(is.finite(mean_Leaf_area_corr_day)) %>%
#   filter(mean_Leaf_area_corr_day > 0)
# save(LB_MEAN_plot_day, file= "mean_Leaf_area_corr_day.RData")

##make link between  pe_data (after spatial correction) and  pe.df.ETr dataset used in Kar et al. 2020
## "Sector", "Genotype", "Replicates", "LeafArea3D", "TS", "date"
spatCorrSpGDP$date= date(spatCorrSpGDP$timePoint)
spatCorrSpGDP$date = ymd(spatCorrSpGDP$date)

spatCorrSpGDP= spatCorrSpGDP %>% inner_join(d_exp, by = c('plotId'))
# Rename as the following columns: "Sector", "Genotype", "Replicates", "LeafArea3D", "TS", "date"
colnames(spatCorrSpGDP)[2]= "TS"
colnames(spatCorrSpGDP)[15]= "Genotype"
colnames(spatCorrSpGDP)[18]= "Sector"
colnames(spatCorrSpGDP)[13]= "Replicates"
colnames(spatCorrSpGDP)[3]= "LeafArea3D" ## PE CORRECTED VALUES
colnames(spatCorrSpGDP)[4]= "LeafArea3D_raw" ## PE RAW VALUES

spatCorrSpGDP <- spatCorrSpGDP %>% select(Sector, Genotype, Replicates, LeafArea3D, 
                                          LeafArea3D_raw, TS,date)

pe.df.ETr= spatCorrSpGDP

colnames(pe.df.ETr) <- c("old_unit", "Genotype", "Replicates", "LeafArea3D", "LeafArea3D_raw", "TS", "date")

save(pe.df.ETr, file = paste(opPATH.obj,"pe.df.ETr.RData", sep = ""))
load(file =  paste(opPATH.obj,"pe.df.ETr.RData", sep = ""))

pe.df.ETr$Genotype <- factor(pe.df.ETr$Genotype)
pe.df.ETr$Replicates <- factor(pe.df.ETr$Replicates)
pe.df.ETr$date=as.factor(pe.df.ETr$date)


# ## should we used 3D LA corrected or 3D-LA raw : ANOVA LA - repetition significant or not ?
# for (i in levels(pe.df.ETr$date)) {
#   
#   exemple_i <- pe.df.ETr %>% 
#     dplyr::filter(pe.df.ETr$date== i) 
#   
#   model=lm(LeafArea3D~Replicates*Genotype, data = pe.df.ETr)
#   aov<-Anova(model)
#   print(aov)
#   
# }


####### STAGE IV: Calculation of LAI per day, LAI per TS, observed LA, Transpiration and TRrate ####### 


# conversion3D-LA --> observed LA

# create files adapted for Kar et al pipeline LAI.mat, sel.secs, subs.d
load(file = "./saved_objects/ETr_Meta_ERRsec.rmvd.RData")
load(file = "./saved_objects/pe.df.ETr.RData")
load(file = "./saved_objects/TS.RData")
load(file = "./saved_objects/wthr.ETref.ETobs.RData")
ETr_smoothFILE = ETr_Meta_ERRsec.rmvd

sel.secs <- unique(na.omit(ETr_smoothFILE$old_unit)) 
length(sel.secs) 
save(sel.secs, file = "./saved_objects/sel.secs.RData")

length(TS)
TS= as.character(TS)

LAI.mat <- matrix(NA, nrow = length(sel.secs), ncol = length(TS))
rownames(LAI.mat) <- sel.secs
colnames(LAI.mat) = TS
save(LAI.mat, file = "./saved_objects/LAI.mat.RData")

unq.dts <- unique(substr(colnames(wthr.ETref.ETobs)[-c(1:6)], 1, 10))
length(unq.dts)

LAI.mat <- LAI.mat[order(match(rownames(LAI.mat), 
                               ETr_smoothFILE[(9:nrow(ETr_smoothFILE)), "old_unit"])),]

LAI.all.dates <- LAI.mat; unq.dts.copy <- unq.dts

LAI.all.dates <- LAI.all.dates[ , c(1:length(unq.dts))]
colnames(LAI.all.dates) <- c(as.character(unq.dts.copy))


# Calculate LAI value per day #
################################ 1. CALCUL LEAF AREA INDEX PER DAY    ################################ 
#   objectives : 
# - calculate the median of 3D PE, converted 3D-PE to observed LA with the relation of Vadez 2015  
# - insert planimeter value 
# - calculate the  relation of LA over time with  observed LA and planimeter value +  save the slope
# - replace NA value of LA (from last 3D PE scan to planimeter) with previous value +slope
# - plot LAI per day


# 1) LAI when we have PE value (from 28.09.23 to 11.10.23) 
for(i in 1:nrow(LAI.mat)){ ## for all unit of LAI.mat 
  
  la.df.tmp <- pe.df.ETr[pe.df.ETr$old_unit == rownames(LAI.mat)[i], ] ## select unit in common between pe.data and LAI
  
  la.df.med <- la.df.tmp %>% group_by(date) %>%  ## group by date and provide median of LeafArea3D
    dplyr::summarise(Med_LeafArea3D = median(LeafArea3D))
  
  date.mat <- data.frame(date = unq.dts, val = rep(NA, length(unq.dts))) ## create df with date, val (empty, length unique date)
  date.mat$date = ymd(date.mat$date)
  for(j in 1:nrow(date.mat)) ## for all row of this date.mat df
  {
    ifelse(date.mat$date[j] %in% la.df.med$date, 
           date.mat$val[j] <- la.df.med$Med_LeafArea3D[la.df.med$date == date.mat$date[j]],
           date.mat$val[j] <- NA)    ## fill in df with val= median per date
  }
  
  # Relation from Vadez  et al. 2015 LeasyScan: a novel concept ... )
  # y=  3DLA x= Observed LA. Relation species-dependant
  # y=0.36x +29.9             --> Observed LA(x) = (3DLA_mm²/0.36)-29.9 for pearl millet
  # y=0.36x +63.2             --> Observed LA(x) = (3DLA_mm²/0.36)-63.2  for peanut
  # y=0.39x + 302.2           --> Observed LA(x) = (3DLA_mm²/0.39)-302.2 for cowpea
  
  # * 0.26 m² to measure LAI  #!! careful; 0.26m² = surface of LC in ICRISAT !! not suitable for IRD platform
  
  sec.lai.dates <- rep(date.mat$val, each =1) 
  LAI.all.dates[i, ] <- ((((sec.lai.dates/100) - 29.9)/0.36)*(1/0.26)/10000) ##for pearl millet 
  
  
  
}
# LAI per day with NA value until limit PE

## 2) insert planimeter value and linked with PE value
### insert planimeter values on excel
LAI = read.csv(file = "./data/LAI_converted_with_NA_planimeter.csv", sep = ";") # (correspond à la date du 19/10/2023)
rownames(LAI)= LAI[,1]
LAI = LAI %>% arrange(LAI$X) ##trie par ordre croissant de nom de LC
LAI$X10.10.2023= NA 
LAI$X11.10.2023 = NA ## les données du 10 et 11 oct sont outliers dans certains cas de nombreux cas de  LC, mieux vaut les enlever pour améliorer la relation
colnames(LAI)[9:ncol(LAI)]= sub("X","",colnames(LAI)[9:ncol(LAI)])

## removed 38:10:01 and 38:10:02 (row 854 and 855) because there are just planimeter data (PE 3DLA are missing data) 
LAI= LAI[-c(854, 855),]

# add a day to unq.dts for planimeter data 
unq.dts=colnames(LAI)[9:ncol(LAI)]## there is still irrigation dates (22 days)
unq.dts= dmy(unq.dts)
length(unq.dts)
date.mat <- data.frame(date = unq.dts, val = rep(NA, length(unq.dts))) ## create df with date, val (empty, length unique date)
date.mat$date = dmy(date.mat$date)

### replace NA  with last measure +slope value of LAI over time
for (i in 1:nrow(LAI)) {
  data_LC=LAI[i,9:ncol(LAI)]
  date.mat$val= t(data_LC)
  colnames(date.mat)= c("date", "val")
  if (sum(is.na(date.mat$val))!= length(date.mat$val)) {
    mod<-lm(formula =val~date, data = date.mat)
    slope=as.numeric(mod$coefficients[2])
    
    for(j in 1:length(date.mat$val)){
      if(is.na(date.mat$val[j])== TRUE){
        date.mat$val[j]= date.mat$val[j-1]+slope} else {j=j+1}
    }
    LAI[i,9:ncol(LAI)]=date.mat$val
    print(i)
  }
}
LC= LAI[,1]
rownames(LAI) = LC 
TS= colnames(LAI[9:ncol(LAI)])
TS= dmy(TS)

save(LC, file = paste(opPATH.obj,"LC.RData", sep = ""))
save(unq.dts, file =  paste(opPATH.obj,"unq.dts_for_LAI_values.RData"))
save(LAI,  file = paste(opPATH.obj,"LAI_per_day.RData", sep = ""))
write.csv2(LAI, file = "LAI_observed_converted3D_plani_estim.csv")


max_df=as.data.frame(colMax(LAI[,9:ncol(LAI)])) 
max(max_df[1:nrow(max_df),]) 
min_df=as.data.frame(colMin(LAI[,9:ncol(LAI)]))
min(min_df[1:nrow(min_df),]) 


## removed negative LAI if present 
for (j in (9:ncol(LAI))){ 
  for (i in 1:nrow(LAI)) { 
    LAI[i,j][LAI[i,j]<0]<- NA 
  }
}

##LAI values per day with LAI converted until 9 oct + slope value + planimeter
data=as.data.frame(t(LAI))
data= data[-c(1:8),]
rownames(data)= TS
data= cbind(TS, data)
data_num <- as.data.frame(apply(data[,2:ncol(data)], 2, as.numeric))  # Convert all variable types to numeric
data= cbind(data[1], data_num)

i=2
# for ( i in 2: ncol(data)) {
  png(file = paste(opPATH.profile,i-1, "LAI_per_day.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=data[,1]))+ ## data[,1] = TS
    ylab("Leaf area index ")+
    xlab("Timestamp")+
    labs(title = paste("Leaf area index  ",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "darkred")+
    scale_x_date(date_labels = "%d-%m", date_breaks = "1 day")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }

max_df=as.data.frame(colMax(LAI[,9:ncol(LAI)])) 
max(max_df[1:nrow(max_df),]) 
min_df=as.data.frame(colMin(LAI[,9:ncol(LAI)]))
min(min_df[1:nrow(min_df),]) 
save(LAI, file =  paste(opPATH.obj,"LAI_rmd_out.RData", sep = ""))

################################ 2. CALCUL LEAF AREA INDEX PER INTERVAL of 15 min    ################################ 
val_per_interval= data.frame(matrix(NA, nrow = nrow(LAI), ncol = ncol(LAI[9:ncol(LAI)])*96)) ## on commence à 9 car debut de la meure /jour, 96 nb fof interval per day !! to adapt to other 
rownames(val_per_interval)=rownames(LAI)

## create TS all day
start_seqby= ymd_hms(paste(dmy(unq.dts[1]),as.character("00:00:00")))
end_seqby= ymd_hms(paste(dmy(unq.dts[length(unq.dts)]),as.character("23:45:00")))
start_seqby
end_seqby

TS_seqby= seq(start_seqby,end_seqby, by = '15 mins')
length(TS_seqby)

a=1
b=96 #= seq_by*60

for (i in 9:ncol(LAI)) { ## les 8 premières lignes sont des metad
  # for (i in 1:ncol(LAI.mat)) {
  sec.lai.tmp <- as.data.frame(LAI[,i])
  colnames(sec.lai.tmp) = colnames(LAI)[i]
  sec.lai.tmp = cbind(sec.lai.tmp, replicate(95,sec.lai.tmp)) #replicate the column LA per day * 96 interval of 15 min
  
  val_per_interval[a:b]= sec.lai.tmp
  # colnames(val_per_interval)[a:b]= colnames(sec.lai.tmp)
  a=a+96
  b=b+96
}
colnames(val_per_interval)=TS_seqby
save(val_per_interval, file =  paste(opPATH.obj,'LAI_per_15mins.RData', sep = ""))
write.table(val_per_interval, file = "LAI_per_15mins.csv", row.names = F, sep = ";")

## removed irrigation dates
load(file = "./saved_objects/unq.dts.RData")
# unq.dts = dmy(unq.dts)
unq.dts_irrg_rvd= unq.dts[unq.dts !=  c("2023-09-30" )]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-03")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-06")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-10")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-13")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-17")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-18")]
unq.dts_irrg_rvd= unq.dts_irrg_rvd[unq.dts_irrg_rvd !=  c("2023-10-19")] ## planimeter date
unq.dts_irrg_rvd
save(unq.dts_irrg_rvd, file=paste(opPATH.obj, "unq.dts_irrg_rvd.RData", sep = ""))


## create TS per 15 min on the entire trial except irrigation dates
i=1
## do the first day  i=1
start_day =  ymd_hms(paste(ymd( unq.dts_irrg_rvd[i]),as.character("00:00:00")))
start_day
end_day = ymd_hms(paste(ymd( unq.dts_irrg_rvd[i]),as.character("23:45:00")))
end_day
date = as.data.frame(seq(start_day,end_day, by = '15 mins'))
TS= date

for (i in 2:length(unq.dts_irrg_rvd)) {
  start_day =  ymd_hms(paste(ymd( unq.dts_irrg_rvd[i]),as.character("00:00:00")))
  end_day = ymd_hms(paste(ymd( unq.dts_irrg_rvd[i]),as.character("23:45:00")))
  date = as.data.frame(seq(start_day,end_day, by = '15 mins'))
  TS= rbind(TS, date)
}
TS= TS[1:nrow(TS),]
save(TS, file =paste(opPATH.obj,  "TS_entire_trial_rmd_irrig.RData", sep = ""))

save(val_per_interval, file=paste(opPATH.obj, "LAI_per_15min_entire_trial.RData", sep = ""))

##val_per_interval contains values of LAI for all timepoints ( length(TS_seqby) = 2112), we want to keep only Timestamps TS,  length(TS)= 1344, with irrigation dates removed
LAI.mat= val_per_interval[,as.character(TS)]
ncol(LAI.mat)

save(val_per_interval, file="LAI_per_15min_rmd_irrig.RData")



save(LAI.mat, file = paste (opPATH.obj, "LAI.mat.RData", sep = ""))
write.table(LAI.mat, file = "LAI.mat.csv", sep = ";", dec = ".")


##### plot LAI over time to check the profile
data=as.data.frame(t(LAI.mat))

TS=rownames(data)
data$TS= TS
data$TS =ymd_hms(data$TS)

i=1
# for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "LAI_over_time.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("Leaf area index ")+
    xlab("Timestamp")+
    labs(title = paste("Leaf area index  ",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "darkred")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }


max_df=as.data.frame(colMax(LAI.mat)) 
max(max_df[1:nrow(max_df),]) ## 
min_df=as.data.frame(colMin(LAI.mat))
min(min_df[1:nrow(min_df),]) ## 
# ## REMOVED OUTLIERS INF if Negative LAI value
# for  (j in 1:ncol(LAI.mat_rmd_out)){ 
#   for (i in 1:nrow(LAI.mat_rmd_out)) { 
#     LAI.mat_rmd_out[i,j][LAI.mat_rmd_out[i,j]<0]<- NA     }
#   print(j)
# }
# 
# LAI.mat_rmd_out= LAI.mat


################################ 3. CALCUL LEAF AREA ################################ 
# Leaf Area Index  = Leaf Area /surface de LC 
# Leaf Area =Leaf Area Index * surface LC 
# LA = LAI * 0.26 m²
LEAF_AREA= data.frame(matrix(NA, nrow = nrow(LAI.mat), ncol =  ncol(LAI.mat)))
for (i in 1:nrow(LAI.mat)){
  LEAF_AREA[i, ] <- LAI.mat[i, ]* 0.26 
  print(i)
}
colnames(LEAF_AREA) = TS 

## !! now all row in LC ascending order : should start with 131:07:02
LC = rownames(LAI.mat)
rownames(LEAF_AREA)= LC
colnames(LEAF_AREA)=TS
save(LEAF_AREA,  file = paste (opPATH.obj, "LEAF_AREA.RData", sep = ""))

write.table(LEAF_AREA, file = "LEAF_AREA.csv", sep = ";", dec = ".")


## removed outliers 
max_df=as.data.frame(colMax(LEAF_AREA)) 
max(max_df[1:nrow(max_df),]) #
min_df=as.data.frame(colMin(LEAF_AREA))
min(min_df[1:nrow(min_df),]) ##  

save(LEAF_AREA, file = paste(opPATH.obj,"LEAF_AREA_rmd_out.RData", sep = ""))
write.table(LEAF_AREA, file = "LEAF_AREA_rmd_out.csv", sep = ";", dec = ".")


data=as.data.frame(t(LEAF_AREA))
colnames(data)= LC
rownames(data)=TS
data$TS= TS
data$TS =ymd_hms(data$TS)

for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "LA_over_time.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("Leaf area (m²) ")+
    xlab("Timestamp")+
    labs(title = paste("Leaf area ",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "darkgreen")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
  print(i)
}


################################ 4. ETr NON NORMALIZED ############################### 
load(file = "ETr_filtered.RData")
x= ETr_filtered
ETr_filt_imputed_FILE <- x
ETr_smth.mat <- ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE), 
                                      9:ncol(ETr_filt_imputed_FILE)]

rownames(ETr_smth.mat)=   ETr_filt_imputed_FILE[9:nrow(ETr_filt_imputed_FILE),1]
ETr_smth.mat= cbind(rownames(ETr_smth.mat), ETr_smth.mat)
ETr_smth.mat<-ETr_smth.mat %>% arrange(`rownames(ETr_smth.mat)`) ##trie par ordre croissant de nom de LC
rownames(ETr_smth.mat)= ETr_smth.mat[,1]
ETr_smth.mat= ETr_smth.mat[,-1]
colnames(ETr_smth.mat)=TS

## add ETref values
ETref=x[8,9:ncol(x)]
colnames(ETref)= TS
ETr_smth.mat = rbind(ETref,ETr_smth.mat )
rownames(ETr_smth.mat)[1]= "ETref" 
rownames(ETr_smth.mat)=LC
# ETr_smth.mat=ETr_smth.mat[,-]

max_df=as.data.frame(colMax(ETr_smth.mat)) 
max(max_df[1:nrow(max_df),]) ##  0.122
min_df=as.data.frame(colMin(ETr_smth.mat))
min(min_df[1:nrow(min_df),]) ## -0.01


## removed negative values
for (i in 1:nrow(ETr_smth.mat)){
  for (j in 1:ncol(ETr_smth.mat)){
    ETr_smth.mat[i,j][ETr_smth.mat[i,j]<0]<- NA } }
load( file = "ETr_smth.mat_25_03_rmd_negative.RData,")

save(ETr_smth.mat, file = "ETr_smth.mat_25_03_rmd_negative.RData,")
write.table(ETr_smth.mat, file = "ETr_smth.mat_update_25_03_rmd_negative.csv", sep = ";", dec = ".")

ETref= t(x[8,9:ncol(x)])
data=as.data.frame(t(ETr_smth.mat))
colnames(data)= LC
data=  as.data.frame(apply(data, 2, as.numeric))
data = cbind(ETref, data)
colnames(data)= c("ETref",LC)

data$TS= TS
data$TS =ymd_hms(data$TS)
i=2
for ( i in 2: ncol(data)) {
  png(file = paste(opPATH.profile,i-1, "_ETr_non_normalized_ETref.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("ETr and ETref (mm.15min-1) ")+
    xlab("Timestamp")+
    labs(title = paste("ETr filtered (non normalized) ",colnames(data)[i]))+
    geom_line(aes(y = data[,1]), color = "darkolivegreen")+
    geom_point(aes(y = data[,i]), color = "darkblue")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
  print(i)
}


################################ 5. ETr  NORMALIZED BY LA ############################### 
ETr_smth.mat = ETr_smth.mat[-1,] ## enlever la premiere ligne ETref
## removed 38:10:01 and 38:10:02 (row 854 and 855) because there are just planimeter (PE 3DLA missing data) 
rownames(ETr_smth.mat)[854]
rownames(ETr_smth.mat)[855]

ETr_smth.mat= ETr_smth.mat[-c(854, 855),]
save(ETr_smth.mat, file = "ETr_smth.mat.RData")
LC=rownames(ETr_smth.mat)
## length(LC) 1204 LC ok !

### normalization by LA 
ETR_normalized= data.frame(matrix(NA, nrow = nrow(ETr_smth.mat), ncol =  ncol(ETr_smth.mat)))
for (i in 1:nrow(ETr_smth.mat)){
  for(j in 1:ncol(ETr_smth.mat)) {
    ETR_normalized[i, j] <- ETr_smth.mat[i,j ]/ LEAF_AREA[i,j]
    ETR_normalized[i,j][ETR_normalized[i,j]<0]<- NA ## if the ETr is negative
  }
  print(i)
}
colnames(ETR_normalized)= TS
rownames(ETR_normalized)= LC


max_df=as.data.frame(colMax(ETR_normalized))
max(max_df[1:nrow(max_df),]) 
min_df=as.data.frame(colMin(ETR_normalized)) ##is this maximum weight value is possible?  or outlier?
min(min_df[1:nrow(min_df),]) 

save(ETR_normalized, file = "ETR_normalized_data_normalized_LA.RData")
write.table(ETR_normalized, file = "ETR_normalized_by_LA.csv", sep = ";", dec = ".")



### plot ETr normalized by LA 

data=as.data.frame(t(ETR_normalized))
# data=as.data.frame(t(ETR_normalized))
data$TS= TS
data$TS =ymd_hms(data$TS)

for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "_ETR_normalized.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("ETr normalized (mm.15min-1.m-2) ")+
    xlab("Timestamp")+
    labs(title = paste("ETR_normalized by LA",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "chocolate4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
  print(i)
}
# 
## We observed high peak that does not seem logical. It is not always the same day (not related to weather condition).  not present for all genotypes
## Objective : removed these peaks based on max per day 


data=as.data.frame(t(ETR_normalized))
data$TS= TS
data$date= date(TS)
data$date= as.factor(data$date)

## prepare dataframe unique_date used in the loop to call date and n_day
unique_date= data.frame(matrix( nrow = length( unique(data$date)), ncol = 2))
colnames(unique_date)=c ("n_day", "day")
unique_date$n_day= (1:nlevels(data$date))
unique_date$day= levels(data$date)
unique_date$day= date( unique_date$day)

## prepare dataframe max_day to store value
max_day= data.frame(matrix(NA, nrow = length( unique(data$date)), ncol = 2)) ## create a empty dataframe to store data
colnames(max_day)=c ("date", "max")
max_day$max = as.numeric(max_day$max)
max_day$date = date(max_day$date)

dat_rmd= data.frame(matrix(NA, nrow = length(LC), ncol = 1)) ## create a empty dataframe to store data
rownames(dat_rmd)=LC


for (n in (1:(ncol(data)-2))) { # for each LC (2 last col = TS and date)
  print(n)
  subset_LC=data.frame(data$TS, data$date, data[n]) # create dataframe with only TS, date and LC n
  colnames(subset_LC)= c("TS", "date", colnames(data)[n])
  subset_LC$date = date(subset_LC$date)
  
  for (i in 1:nrow(unique_date) ){ # for each day
    
    subset_day<- subset_LC[subset_LC$date %in% c(unique_date[i,2]), ] ## create a dataframe with only one day
    max_day[i,1]= unique_date[i,2]
    max_day[i,2]= max(subset_day[,3])
  }
  ol <- boxplot(max_day$max, plot = FALSE)$out
  tf <- which(max_day$max %in% ol)
  date_to_rmd = max_day$date[tf]
  date_to_rmd
  data$date=date(subset_LC$date)
  
  if (length(date_to_rmd) == 0){
    print("all days of data have been kept")
    dat_rmd[n,] =0
    
  } else {
    for (date in 1:nrow(data)) {
      for (nday_to_rmd in 1:length(date_to_rmd)){
        if (data$date[date]  == date_to_rmd[nday_to_rmd]) {
          row_to_rmd=  which(data$date %in% date_to_rmd)
          data[row_to_rmd,n]=NA
          print(paste(length(date_to_rmd),"day(s) have been removed"))
          dat_rmd[n,] = length(date_to_rmd)}
      }
    }
  }
}


ETR_normalized= as.data.frame(t(data))
ETR_normalized= ETR_normalized[-nrow(ETR_normalized),] # removed TS nd date
ETR_normalized=  as.data.frame(apply(ETR_normalized, 2, as.numeric))
rownames(ETR_normalized)= LC

save(ETR_normalized, file = "ETR_normalized_rmd_out.RData")



for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "_ETR_normalized_rmd_out.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("ETr normalized rmd out (mm.15min-1.m-2) ")+
    xlab("Timestamp")+
    labs(title = paste("ETR_normalized by LA  rmd out",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "chocolate4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
  print(i)
}

################################ 6. TRANSPIRATION 1ST STEP Tr ############################### 
Trans.mat= data.frame(matrix(NA, nrow = nrow(ETr_smth.mat), ncol =  ncol(LAI.mat))) ### ncol 96*nombre jour 14 
rownames(Trans.mat)=LC
colnames(Trans.mat)=TS

for (i in 1:nrow(ETr_smth.mat)){
  if (rownames(ETr_smth.mat)[i] ==rownames(LAI.mat)[i]){
    Trans.mat[i, ] <- (1-exp(-0.463*LAI.mat[i, ]))*ETr_smth.mat[i, ] #!! why 1-(1-)
    print(i)
  }
}

colnames(Trans.mat)= TS
rownames(Trans.mat)=LC

for (j in (1:ncol(Trans.mat))){ ## i col LC
  for (i in 1:nrow(Trans.mat)) { ## j row TS
    Trans.mat[i,j][Trans.mat[i,j]<0]<- NA   }
  
}

save (Trans.mat, file = "Transpiration_Trans.mat.RData")
write.table(Trans.mat, file = "Trans.mat.csv", sep = ";", dec = ".")


max_df=as.data.frame(colMax(Trans.mat))
max(max_df[1:nrow(max_df),]) #   
min_df=as.data.frame(colMin(Trans.mat)) 
min(min_df[1:nrow(min_df),]) #   

# update 20.03 rmd outliers : 

max_df=as.data.frame(colMax(Trans.mat))
max(max_df[1:nrow(max_df),]) #   0.07207621
min_df=as.data.frame(colMin(Trans.mat)) 
min(min_df[1:nrow(min_df),]) # 0


Trans.mat= Trans.mat[-ncol(Trans.mat)]
save (Trans.mat, file = "Transpiration_Trans.mat_rmd_out.RData")
write.table(Trans.mat, file = "Trans.mat_rmd_out.csv", sep = ";", dec = ".")

## add VPD 
VPD=x[3,7:ncol(x)] # LENGTH + 1344 TS

data=as.data.frame(t(Trans.mat))
data=cbind(t(VPD), data)
colnames(data)[1]="VPD"

rownames(data)=TS
data <- as.data.frame(apply(data, 2, as.numeric))
data$TS= TS
data$TS =ymd_hms(data$TS)

i=2
# for ( i in 2: ncol(data)) {
  png(file = paste(opPATH.profile,i-1, "_Transpiration_VPD.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    geom_point(aes(y = data[,1]), color = "black")+
    geom_point(aes(y = data[,i]*100), color = "aquamarine4")+
    scale_y_continuous(name="VPD (kPa)", sec.axis=sec_axis(~./100, name="Transpiration (mm.15min-1)")) +
    xlab("Timestamp")+
    labs(title = paste("Tr (LAI~ETr ) & VPD",colnames(data)[i]))+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }
# 


i=2
# for ( i in 2: ncol(data)) {
  png(file = paste(opPATH.profile,i-1, "_Transpiration.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    geom_point(aes(y = data[,i]), color = "aquamarine4")+
    ylab("Transpiration (mm.15min-1)")+
    xlab("Timestamp")+
    labs(title = paste("Tr (LAI~ETr )",colnames(data)[i]))+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }


### plot ETr normalized and Transpiration 
# ETR_normalized= cbind(rownames(ETR_normalized), ETR_normalized)
# colnames(ETR_normalized)[1]="LC"
# rownames(ETR_normalized)= paste("ETr_norm",1:nrow(ETR_normalized), sep = "")

ETr_smth.mat= cbind(rownames(ETr_smth.mat), ETr_smth.mat)
colnames(ETr_smth.mat)[1]="LC"
rownames(ETr_smth.mat)= paste("ETr_",1:nrow(ETr_smth.mat), sep = "")

Trans.mat= cbind(rownames(Trans.mat), Trans.mat)
colnames(Trans.mat)[1]="LC"
rownames(Trans.mat)= paste("Trans_",1:nrow(Trans.mat), sep = "")


## !! pour merger TR rate and Tr, on ne peut pas garde rles LC comme rownames. Rownames doivent être unique
data_ETr_TR= rbind(ETr_smth.mat, Trans.mat)
data_ETr_TR$LC =as.factor(data_ETr_TR$LC)

i=levels(data_ETr_TR$LC)[1]
# for ( i in (levels(data_ETr_TR$LC)) ) {  ## for unq_LC = name LC 
  subset_LC<- data_ETr_TR[data_ETr_TR$LC %in% c(i), ] ## create a dataframe with data of one genotype only
  j=rownames(subset_LC)[1]
  j=sub("ETr", "",j)
  subset_LC= t(subset_LC)
  LC_name= subset_LC[1,1]
  subset_LC= subset_LC[-1,]
  subset_LC= as.data.frame(subset_LC)
  subset_LC <- as.data.frame(apply(subset_LC, 2, as.numeric))
  subset_LC$TS = ymd_hms(TS)
  colnames(subset_LC)= c("ETr", "Transp", "TS")
  
  x <- subset_LC$TS  # used to draw the plot 
  y <- subset_LC$Transp
  y2 = subset_LC$ETr
  
  
  png(file = paste(opPATH.profile,j, "_ETr_and_Transp.png", sep = ""))
  plot<-
    ggplot(data =subset_LC, aes(x))+
    ylab("ETr and Tr (mm.15min-1)")+
    xlab("Timestamp")+
    labs(title = paste("Tr (green) and ETr (brown)  ",i))+
    geom_point(aes(y = y2), color = "chocolate4")+
    geom_point(aes(y = y), color = "aquamarine4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }

### plot ETr  and Transpiration 
ETr_smth.mat= cbind(rownames(ETr_smth.mat), ETr_smth.mat)
colnames(ETr_smth.mat)[1]="LC"
rownames(ETr_smth.mat)= paste("ETr",1:nrow(ETr_smth.mat), sep = "")

Trans.mat= cbind(rownames(Trans.mat), ETr_smth.mat)
colnames(Trans.mat)[1]="LC"
rownames(Trans.mat)= paste("Trans",1:nrow(Trans.mat), sep = "")


## !! pour merger TR rate and Tr, on ne peut pas garde rles LC comme rownames. Rownames doivent être unique
data_ETr_TR= rbind(ETr_smth.mat, Trans.mat)
data_ETr_TR$LC =as.factor(data_ETr_TR$LC)

i=levels(data_ETr_TR$LC)[1]
# for ( i in (levels(data_ETr_TR$LC)) ) {  ## for unq_LC = name LC 
  subset_LC<- data_ETr_TR[data_ETr_TR$LC %in% c(i), ] ## create a dataframe with data of one genotype only
  j=rownames(subset_LC)[1]
  j=sub("ETr", "",j)
  subset_LC= t(subset_LC)
  LC_name= subset_LC[1,1]
  subset_LC= subset_LC[-1,]
  subset_LC= as.data.frame(subset_LC)
  subset_LC <- as.data.frame(apply(subset_LC, 2, as.numeric))
  subset_LC$TS = ymd_hms(TS)
  colnames(subset_LC)= c("ETr", "Transp", "TS")
  
  x <- subset_LC$TS  # used to draw the plot 
  y <- subset_LC$Transp
  y2 = subset_LC$ETr
  
  
  png(file = paste(opPATH.profile,j, "_ETr_and_Transp.png", sep = ""))
  plot<-
    ggplot(data =subset_LC, aes(x))+
    ylab("ETr  and Tr (mm.15min-1)")+
    xlab("Timestamp")+
    labs(title = paste("Tr (green) and ETr norm (brown)  ",i))+
    geom_point(aes(y = y), color = "aquamarine4")+
    geom_point(aes(y = y2), color = "chocolate4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%Y-%m-%d"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }

################################ 7. TRANSPIRATION RATE 2ND STEP  ############################### 
rownames(Trans.mat)= Trans.mat[,1]
Trans.mat= Trans.mat[,-1]

Trans.mat <- as.data.frame(apply(Trans.mat, 2, as.numeric))


TR.mat <- data.frame(matrix(NA, nrow = nrow(Trans.mat), ncol =  ncol(Trans.mat))) #Empty Trans Matrix
rownames(TR.mat)= LC
colnames(TR.mat)= TS

# Calculate Transpiration Rate #
#Kar 2020 : "The Tr values thus obtained were converted to Transpiration Rate (TR , i.e. the transpiration per unit of leaf area) by dividing Tr of each day
with the corresponding Observed LA (Eq.6 : TR = Tr/ ObservedLA)
for (i in 1:nrow(TR.mat)){
  # for (j in 1:ncol(TR.mat)) {
  TR.mat[i, ] <- (Trans.mat[i, ]/ (LEAF_AREA[i, ]))
  print(i)
}

for (j in (1:ncol(TR.mat))){ ## i col LC
  for (i in 1:nrow(TR.mat)) { ## j row TS
    TR.mat[i,j][TR.mat[i,j]<0]<- NA   }
  
}

save(TR.mat, file = "TR.mat_Transpiration_rate.RData")
write.table(TR.mat, file = "TR.mat.csv", dec = ".", sep = ";")

max_df=as.data.frame(colMax(TR.mat))
max(max_df[1:nrow(max_df),])
min_df=as.data.frame(colMin(TR.mat)) 
min(min_df[1:nrow(min_df),])  

data=as.data.frame(t(TR.mat))
data=data[-nrow(data),]
rownames(data)=TS
data$TS= TS
data$TS =ymd_hms(data$TS)

i=1
# for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "_TRrate.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("TR rate (mm.15min-1.m-2)")+
    xlab("Timestamp")+
    labs(title = paste("Transpiration rate ",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "deeppink4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }


## still high peaks to removed
data=as.data.frame(t(TR.mat))
data  = data[-nrow(data),]
rownames(data) = TS
colnames(data)= LC
data$TS= TS
data$date= date(TS)
data$date= as.factor(data$date)
nlevels(data$date)


## prepare dataframe unique_date used in the loop to call date and n_day
unique_date= data.frame(matrix( nrow = length( unique(data$date)), ncol = 2))
colnames(unique_date)=c ("n_day", "day")
unique_date$n_day= 1:nlevels(data$date)
unique_date$day= levels(data$date)
unique_date$day= date( unique_date$day)

## prepare dataframe max_dayto store value
max_day= data.frame(matrix(NA, nrow = length( unique(data$date)), ncol = 2)) ## create a empty dataframe to store data
colnames(max_day)=c ("date", "max")
max_day$max = as.numeric(max_day$max)
max_day$date = date(max_day$date)

dat_rmd= data.frame(matrix(NA, nrow = length(LC), ncol = 1)) ## create a empty dataframe to store data
rownames(dat_rmd)=LC


for (n in (1:(ncol(data)-2))) { # for each LC (2 last col = TS and date)
  print(n)
  subset_LC=data.frame(data$TS, data$date, data[n]) # create dataframe with only TS, date and LC n 
  colnames(subset_LC)= c("TS", "date", colnames(data)[n])
  subset_LC$date = date(subset_LC$date)
  
  for (i in 1:nrow(unique_date) ){ # for each day  
    
    subset_day<- subset_LC[subset_LC$date %in% c(unique_date[i,2]), ] ## create a dataframe with only one day
    max_day[i,1]= unique_date[i,2]
    max_day[i,2]= max(subset_day[,3])
  }
  ol <- boxplot(max_day$max, plot = FALSE)$out
  tf <- which(max_day$max %in% ol)
  date_to_rmd = max_day$date[tf]
  date_to_rmd
  data$date=date(subset_LC$date)
  
  if (length(date_to_rmd) == 0){
    print("all days of data have been kept")
    dat_rmd[n,] =0
    
  } else {
    for (date in 1:nrow(data)) {
      for (nday_to_rmd in 1:length(date_to_rmd)){
        if (data$date[date]  == date_to_rmd[nday_to_rmd]) {
          row_to_rmd=  which(data$date %in% date_to_rmd)
          data[row_to_rmd,n]=NA 
          print(paste(length(date_to_rmd),"day(s) have been removed"))
          dat_rmd[n,] = length(date_to_rmd)}
      }
    }
  }
}



i=1
# for ( i in 1: ncol(data)) {
  png(file = paste(opPATH.profile,i, "_TR.mat_rmd_out.png", sep = ""))
  plot<-
    ggplot(data =data, aes(x=TS))+
    ylab("TR rate rmd out (mm.15min-1.m-2) ")+
    xlab("Timestamp")+
    labs(title = paste("TR rate rmd out",colnames(data)[i]))+
    geom_point(aes(y = data[,i]), color = "chocolate4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%d-%m"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }





TR.mat_rmd= as.data.frame(t(data))
TR.mat_rmd= TR.mat_rmd[-nrow(TR.mat_rmd),]
TR.mat_rmd= TR.mat_rmd[-nrow(TR.mat_rmd),] # removed TS nd date 
TR.mat_rmd=  as.data.frame(apply(TR.mat_rmd, 2, as.numeric))
rownames(TR.mat_rmd)= LC

save(TR.mat_rmd, file = "TR.mat_rmd_out.RData,")
write.table(TR.mat_rmd, file = "TR.mat_rmd_out.csv", sep = ";", dec = ".")

## plot TRrate and Tr on the 

Trans.mat= cbind(rownames(Trans.mat), Trans.mat)
TR.mat= cbind(rownames(TR.mat), TR.mat)
colnames(Trans.mat)[1]= "LC"
colnames(TR.mat)[1]= "LC"
rownames(Trans.mat)= paste("Tr",1:nrow(Trans.mat), sep = "")
rownames(TR.mat)= paste("TR_rate",1:nrow(TR.mat), sep = "")

## !! pour merger TR rate and Tr, on ne peut pas garde rles LC comme rownames. Rownames doivent être unique
data_TRrate_TR= rbind(Trans.mat, TR.mat)
data_TRrate_TR$LC =as.factor(data_TRrate_TR$LC)

i=levels(data_TRrate_TR$LC)[1]
# for ( i in (levels(data_TRrate_TR$LC)) ) {  ## for unq_LC = name LC 
  subset_LC<- data_TRrate_TR[data_TRrate_TR$LC %in% c(i), ] ## create a dataframe with data of one genotype only
  j=rownames(subset_LC)[1]
  subset_LC= t(subset_LC)
  LC_name= subset_LC[1,1]
  subset_LC= subset_LC[-1,]
  subset_LC= as.data.frame(subset_LC)
  subset_LC <- as.data.frame(apply(subset_LC, 2, as.numeric))
  subset_LC$TS = ymd_hms(TS)
  colnames(subset_LC)= c("Transpiration", "TRrate", "TS")
  
  x <- subset_LC$TS  # used to draw the plot 
  y <- subset_LC$Transpiration
  y2 = subset_LC$TRrate
  
  
  png(file = paste(opPATH.profile,j, "_TRrate.png", sep = ""))
  plot<-
    ggplot(data =subset_LC, aes(x))+
    ylab("TRrate and Tr (mm.15min-1.m-²)")+
    xlab("Timestamp")+
    labs(title = paste("Tr (green) and TR rate (pink)  ",i))+
    geom_point(aes(y = y2), color = "deeppink4")+
    # geom_point(aes(y = y), color = "aquamarine4")+
    ggplot2::scale_x_datetime(labels = scales::date_format(format = "%Y-%m-%d"), breaks  = "1 days")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
  print(plot)
  dev.off()
#   print(i)
# }





####### STAGE V : Feature Extraction on TRraate profile ####### 


## remettre les dataframe  TR.mat adapté pour l'extraction de features
load(file = "metad.RData")
metad = metad %>% arrange(metad$unit) ##trie par ordre croissant de nom de LC

## enlever les LC manquantes :    ## removed 38:10:01 and 38:10:02 (row 854 and 855) because there are just planimeter (PE 3DLA missing data) 
metad= metad[-855,]
metad= metad[-854,]

## add wther
load(file = "wthr.ETref.ETobs.RData")
weather= wthr.ETref.ETobs[1:8,]
rm(wthr.ETref.ETobs)

## select weather on TS saved
weather= weather[,as.character(TS)]

# TR.mat= TR.mat[,-1]
# TR.mat= rbind(weather,TR.mat) ## weather and Trans merged
# 
# ### add 9 rows with NA values
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# metad <- insertRows(metad, 1 , new = NA)
# 
# TR.mat_metad = cbind(metad, TR.mat)
featuresRES <- getFeatures(x = TR.mat_rmd_out_metad)

allFeatures <- featuresRES$allFeatures

# create H2 dataframe to store H2 est. of each feature for each day
F.He <- as.data.frame(matrix(NA, nrow = length(unq.dts), ncol = 15)) # Date-ROW, feature-COL
colnames(F.He) <- c("maxET", "slope.maxET-6", "slope.07-maxET", "slope.00-07", "slope.19-23:45", "curvmaxET", 
                    "total.auc","auc.10-15", "sd.10-15", "auc.prop.10-15", "auc.07-19", "sd.07-19",  
                    "auc.prop.07-19", "auc.night", "cos.sim.index")
rownames(F.He) <- unq.dts

featureHeRES <- getFeatureHe(x = allFeatures, y = TR.mat_rmd_out_metad, d = unq.dts, p = opPATH.smth)
write.csv(featureHeRES, paste0(opPATH, "TR.mat_rmd_out_metad_featureH2.csv"))


## Prepare data for 'each feature'
maxET <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
slope.maxET.6 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
slope.07maxET <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
slope.00.07 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
slope.19.2345 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
curvmaxET <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
total.auc <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
auc.10.15 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
sd.10.15 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
auc.prop.10.15 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
auc.07.19 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
sd.07.19 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
auc.prop.07.19 <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
auc.night <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))
cos.sim.index <- as.data.frame(matrix(nr = (nrow(TR.mat_rmd_out_metad)-8), nc = length(unq.dts)))

for (j in 1:(nrow(TR.mat_rmd_out_metad)-8)){
  
  for(d in 1:length(unq.dts))
  {maxET[j, d] <- data.frame(allFeatures[[j]][d, 1])
  slope.maxET.6[j, d] <- data.frame(allFeatures[[j]][d, 2])
  slope.07maxET[j, d] <- data.frame(allFeatures[[j]][d, 3])
  slope.00.07 [j, d] <- data.frame(allFeatures[[j]][d, 4])
  slope.19.2345[j, d] <- data.frame(allFeatures[[j]][d, 5])
  
  curvmaxET[j, d] <- data.frame(allFeatures[[j]][d, 6])
  total.auc[j, d] <- data.frame(allFeatures[[j]][d, 7])
  auc.10.15[j, d] <- data.frame(allFeatures[[j]][d, 8])
  sd.10.15 [j, d] <- data.frame(allFeatures[[j]][d, 9])
  auc.prop.10.15[j, d] <- data.frame(allFeatures[[j]][d, 10])
  
  auc.07.19[j, d] <- data.frame(allFeatures[[j]][d, 11])
  sd.07.19[j, d] <- data.frame(allFeatures[[j]][d, 12])
  auc.prop.07.19[j, d] <- data.frame(allFeatures[[j]][d, 13])
  auc.night [j, d] <- data.frame(allFeatures[[j]][d, 14])
  cos.sim.index[j, d] <- data.frame(allFeatures[[j]][d, 15])
  }
  
} 

names(maxET)=names(slope.maxET.6)=names(slope.07maxET)=names(slope.00.07)=names(slope.19.2345)<-unq.dts
names(curvmaxET)=names(total.auc)=names(auc.10.15)=names(sd.10.15)=names(auc.prop.10.15)<-unq.dts
names(auc.07.19)=names(sd.07.19)=names(auc.prop.07.19)=names(auc.night)=names(cos.sim.index)<-unq.dts


write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], maxET)), 
            file= paste0(opPATH.smth, "maxET.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], slope.maxET.6)), 
            paste0(opPATH.smth, "slope.maxET.6.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], slope.00.07)), 
            paste0(opPATH.smth, "slope.00.07.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], slope.07maxET)), 
            paste0(opPATH.smth, "slope.07maxET.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], slope.19.2345)), 
            paste0(opPATH.smth, "slope.19.2345.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], curvmaxET)), 
            paste0(opPATH.smth, "curvmaxET.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], total.auc)), 
            paste0(opPATH.smth, "total.auc.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], auc.10.15)), 
            paste0(opPATH.smth, "auc.10.15.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], sd.10.15)), 
            paste0(opPATH.smth, "sd.10.15.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], auc.prop.10.15)), 
            paste0(opPATH.smth, "auc.prop.10.15.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], auc.07.19)), 
            paste0(opPATH.smth, "auc.07.19.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], sd.07.19)), 
            paste0(opPATH.smth, "sd.07.19.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], auc.prop.07.19)), 
            paste0(opPATH.smth, "auc.prop.07.19.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], auc.night)), 
            paste0(opPATH.smth, "auc.night.csv"), dec = ".", sep = ";")

write.table(as.data.frame(cbind(TR.mat_rmd_out_metad[9:nrow(TR.mat_rmd_out_metad), 1:5], cos.sim.index)), 
            paste0(opPATH.smth, "cos.sim.index.csv"), dec = ".", sep = ";")

end.time <- Sys.time()

print(paste0("Complete processing executed in: ", round((end.time-st.time), 2), "minutes"))



write.table(as.data.frame(TR.mat_rmd_out_metad), 
            paste0(opPATH.smth, "TR_rate_rmd_out.csv"), dec = ".", sep = ";")


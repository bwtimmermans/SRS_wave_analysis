# source("/home/ben/research/NOC/SRS_wave_analysis/comparison_GRL/buoy_trend.R")

# Script to load some SRS data from different satellites and compare with buoy NDBC 46066.
# BT 09/2019

# Libraries.
   library(ncdf4)
   library(TTR)
   library(abind)
   require(zyp)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")

# ================================================================= #
# Buoy data specification.
# ================================================================= #
# Loops for multiple plots (QC etc).
#for (AA in c("H","L")) {
for (AA in c(TRUE,FALSE)) {
#   for (BB in c("1hr","24hr")) {

# Buoy 46066, 52.785N 155.047W
# Pacific.
   #buoy_list <- c("46006")
   #buoy_list <- c("46066")
   #buoy_list <- c("46035")
# Hawaii.
   #buoy_list <- c("51001")
   buoy_list <- c("51003")
   #buoy_list <- c("51004")
# Atlantic.
   #buoy_list <- c("41001")
   #buoy_list <- c("41002")
   #buoy_list <- c("41010")

# Averaging flag.
   #flag_av <- "6hr"
   #flag_av <- "1hr"
   #flag_av <- BB
   flag_av <- "24hr"

# Quality control threshold flag.
# 85%
   #flag_QC_thresh <- "H"
# 75%
   flag_QC_thresh <- "L"
   #flag_QC_thresh <- AA
# Title label and threshold.
   if (flag_QC_thresh == "H") {
      #lab_QC <- "<85% rejection"
      lab_QC <- "11/12 mths, 12 days"
      QC_mth_thresh <- 10
      QC_day_thresh <- 12
      QC_hour_thresh <- 12
   } else {
      #lab_QC <- "<75% rejection"
      lab_QC <- "10/12 mths, 8 days"
      QC_mth_thresh <- 9
      QC_day_thresh <- 8
      QC_hour_thresh <- 6
   }

# Regression covariate.
   lab_reg <- "none"

# Regression on / off.
   #flag_buoy <- FALSE
   flag_buoy <- AA

# Plot buoy monthly averages.
   flag_buoy_monthly <- flag_buoy

# Regression on / off.
   flag_reg_lines <- TRUE

#=======================================================================#
# Buoy data.
#=======================================================================#
# Locations.
# Pacific.
   buoy_46066 <- data.frame(lat=52.785, lon=-155.047, name="46066", nudge=0, lab_x_nudge=1)
   buoy_46006 <- data.frame(lat=40.782, lon=-137.397, name="46006", nudge=0, lab_x_nudge=1)
   buoy_46005 <- data.frame(lat=46.134, lon=-131.079, name="46005", nudge=0, lab_x_nudge=1)
   buoy_46059 <- data.frame(lat=38.094, lon=-129.951, name="46059", nudge=0, lab_x_nudge=1)
   buoy_46035 <- data.frame(lat=57.018, lon=-177.708, name="46035", nudge=0, lab_x_nudge=1)
# Atlantic.
   buoy_42060 <- data.frame(lat=16.406, lon=-63.188, name="42060", nudge=0)
   buoy_41001 <- data.frame(lat=34.502, lon=-72.522, name ="41001", nudge=0)
   buoy_41002 <- data.frame(lat=31.760, lon=-74.840, name ="41002", nudge=0)
   buoy_41009 <- data.frame(lat=28.501, lon=-80.184, name ="41009", nudge=0)
   buoy_41010 <- data.frame(lat=28.878, lon=-78.485, name ="41010", nudge=0)
   buoy_41040 <- data.frame(lat=14.559, lon=-53.073, name ="41040", nudge=0)
   buoy_41043 <- data.frame(lat=21.132, lon=-64.856, name ="41043", nudge=0)
   buoy_41044 <- data.frame(lat=21.575, lon=58.625, name ="41044", nudge=0)
   buoy_41046 <- data.frame(lat=23.832, lon=-68.417, name ="41046", nudge=0)
   buoy_41047 <- data.frame(lat=27.520, lon=-71.53, name ="41047", nudge=0)
# Gulf (Ike).
   buoy_42001 <- data.frame(lat=25.897, lon=-89.668, name = "42001", nudge=0)
   buoy_42002 <- data.frame(lat=26.091, lon=-93.758, name = "42002", nudge=0)
   buoy_42003 <- data.frame(lat=26.007, lon=-85.648, name = "42003", nudge=0)
   buoy_42019 <- data.frame(lat=27.907, lon=-95.352, name = "42019", nudge=0)
   buoy_42036 <- data.frame(lat=28.501, lon=-84.516, name = "42036", nudge=0)
   buoy_42035 <- data.frame(lat=29.232, lon=-94.413, name = "42035", nudge=0)
   buoy_42039 <- data.frame(lat=28.788, lon=-86.008, name = "42039", nudge=0.5)
   buoy_42040 <- data.frame(lat=29.208, lon=-88.226, name = "42040", nudge=1)
   buoy_42055 <- data.frame(lat=22.120, lon=-93.96, name = "42055", nudge=0)
# Hawaii.
   buoy_51004 <- data.frame(lat=17.604, lon=-152.364, name = "51004", nudge=0)
   buoy_51003 <- data.frame(lat=19.172, lon=-160.662, name = "51003", nudge=0)
   buoy_51001 <- data.frame(lat=24.453, lon=-162.0, name = "51001", nudge=0)
#-----------------------------------------------------------------------#

# Load historical data for NDBC buoys.
   mat_buoy_obs <- matrix(0,nrow=500000,ncol=16)
   mat_data_dims <- matrix(0,ncol=2,nrow=length(buoy_list))
# Loop over buoys.
   for (b.idx in 1:length(buoy_list)) {
      #buoy_name <- c("46066")
      buoy_name <- buoy_list[b.idx]
# File path.
      buoy_data_dir <- paste("/home/ben/research/waves/buoy_data/NDBC",buoy_name,sep="")
      data_files <- list.files(path = buoy_data_dir, pattern = paste("^",buoy_name,"*",sep="") )
      no_files <- length(data_files)

      file_path <- 0
      data_dims <- matrix(0,nrow=no_files,ncol=2)
      mat_buoy <- 0

      for (i in 1:no_files) {
         file_path[i] <- paste(buoy_data_dir,"/",data_files[i],sep='')
         temp_table <- read.table(file_path[i],skip=1)
         data_dims[i,] <- dim(temp_table)
         mat_buoy <- rbind(mat_buoy,temp_table[,c(1:16)])
      }

      mat_data_dims[b.idx,] <- dim(mat_buoy)
      #array_buoy_obs[(1:dim(mat_buoy)[1]),,b.idx] <- as.matrix(mat_buoy)
      mat_buoy_obs[(1:dim(mat_buoy)[1]),] <- as.matrix(mat_buoy)
   }

# Remove missing data.
   mat_buoy_obs[mat_buoy_obs==99.00] <- NA
   mat_buoy_obs[mat_buoy_obs==0.00] <- NA

# Adjust data for the change in data formats in 2005.
   hist_data2 <- seq(which(mat_buoy_obs[,1] == 2005)[1],dim(mat_buoy_obs)[1])
   hist_data1 <- seq(1,hist_data2[1]-1)
   hist_buoy_hs <- c(mat_buoy_obs[hist_data1,8],mat_buoy_obs[hist_data2,9])

# Quantiles for 
   q_plot <- quantile(hist_buoy_hs,probs=c(0.5,0.9),na.rm=T)

# Bin data by year.
   if (mat_buoy_obs[3,1] < 100) {
      start_year <- mat_buoy_obs[3,1] + 1900
   } else {
      start_year <- mat_buoy_obs[3,1]
   }
   seq_years <- start_year:2018

#-----------------------------------------------------------------------#
# Gridded data.
#-----------------------------------------------------------------------#
# Labels.
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   list_names <- list("CCI","SRS","ERA","ERA","NOC","GOW")
   list_var_lab <- list("swh_mean","swh_mean","swh_mean","swh_avg","hs_mean","swh_mean")
   vec_var_col <- c(2,5,1,1,1,1)
# Adjust longitude iand latitude for comparison with buoy.
   vec_lon_adjust <- c(0,360,360,360,0,360)
   vec_lat_adjust <- c(0,0,0,0,0,0)

   list_data_path <- list(6)
# Data path CCI.
   #list_data_path[[1]] <- "/home/ben/research/NOC/SRS_wave_analysis/CCI/analysis/output/4deg/list_stats_1992-2018_annual_n_summer_none.Robj"
   list_data_path[[1]] <- "/home/ben/research/NOC/SRS_wave_analysis/CCI/analysis/output/1deg/list_stats_1992-2018_annual_n_summer_none.Robj"

# Data path Young [ALL].
   #list_data_path[[2]] <- paste("/home/ben/research/NOC/SRS_wave_analysis/analysis_NOC/output/",res,"deg/list_stats_",data_type,"_",paste(vec_datasets[mission_idx],collapse='_'),"_",lab_years,"_",lab_months,"_global_",flag_reg,".Robj",sep="")
   #list_data_path[[2]] <- "/home/ben/research/NOC/SRS_wave_analysis/analysis_NOC/output/4deg/list_stats_KU-all_ERS-1_TOPEX_ERS-2_GFO_JASON-1_ENVISAT_JASON-2_CRYOSAT-2_HY-2_1992-2017_annual_global_none.Robj"
   list_data_path[[2]] <- "/home/ben/research/NOC/SRS_wave_analysis/analysis_NOC/output/1deg/list_stats_KU-all_ERS-1_TOPEX_ERS-2_GFO_JASON-1_ENVISAT_JASON-2_CRYOSAT-2_HY-2_1992-2017_annual_global_none.Robj"

# Data path ERA5.
   #list_data_path[[3]] <- "/home/ben/research/NOC/SRS_wave_analysis/ERA5/analysis/output/4deg/list_stats_1980-2018_annual_n_summer_none_mean.Robj"
   list_data_path[[3]] <- "/home/ben/research/NOC/SRS_wave_analysis/ERA5/analysis/output/1deg/list_stats_1979-2018_annual_n_summer_none_mean.Robj"

# Data path ECMWF
   #list_data_path[[4]] <- "/home/ben/research/NOC/SRS_wave_analysis/ECMWF/output/4deg/list_stats_1979-2017_annual_n_summer_none_mean.Robj"
   list_data_path[[4]] <- "/home/ben/research/NOC/SRS_wave_analysis/ECMWF/output/1deg/list_stats_1979-2017_annual_n_summer_none_mean.Robj"

# Data path NOC_WW3.
   list_data_path[[5]] <- "/home/ben/research/NOC/SRS_wave_analysis/NOC_WW3/analysis/output/4deg/list_stats_1992-2015_annual_n_summer_none_1deg.Robj"

# Data path GOW.
   list_data_path[[6]] <- "/home/ben/research/NOC/SRS_wave_analysis/GOW/analysis/output/4deg/list_stats_1990-2018_annual_n_summer_none.Robj"

# Labels.
   plot_labels <- c("(A) CCI [1992-2017]","(B) Young [1992-2017]","(C) ERA5 [1980-2018]","(D) ERA5 WAM (Bidlot) [1979-2017]","(E) NOC WW3 [1992-2015]","(F) GOW WW3 [1990-2018]")

#-----------------------------------------------------------------------#
# Product selection.
   vec_prod <- c(1,2,3,4,5,6)
# CCI, Young, ERA, ECMWF
   vec_prod <- c(1,2,3,4)

# Data structures.   
   list_lab_dataset <- list(6)
   list_buoy_series <- list(6)

# Load data sets in a loop.
   for (d_idx in vec_prod) {
# Matrix to hold datasets.
      attach(list_data_path[[d_idx]])
      attached_data <- eval(parse(text=paste("list_",list_names[[d_idx]],"_stats",sep="")))
      detach(pos=2)

# Set up data structures.
#   mat_plot_data_1 <- matrix(NA,nrow=length(attached_data[[1]]$lat_cell),ncol=length(attached_data[[1]]$lon_cell))
      list_lab_dataset[[d_idx]] <- attached_data[[1]]$dataset_name
      lat_mid <- attached_data[[1]]$lat_mid
      #print(paste("Dataset",list_lab_dataset[[d_idx]]))
      #print(lat_mid)
      lon_mid <- attached_data[[1]]$lon_mid
      #print(lon_mid)
      #print(" ")
      lab_stats <- attached_data[[1]]$trend_stats

# Assign data.
      list_data <- attached_data[[2]]
      rm(attached_data)

# Find longitude and latitude.
      lat_cell_temp <- abs(lat_mid - eval(parse(text=paste("buoy_",buoy_list[1],"$lat",sep=""))))
      buoy_lat_cell_idx <- which( lat_cell_temp == min(lat_cell_temp) )
      lon_cell_temp <- abs(lon_mid - (eval(parse(text=paste("buoy_",buoy_list[1],"$lon",sep=""))) + vec_lon_adjust[d_idx]))
      buoy_lon_cell_idx <- which( lon_cell_temp == min(lon_cell_temp) )[1]
# Fix for Young.
      if ( d_idx == 2 ) {
         series_temp <- list_data[[buoy_lat_cell_idx,buoy_lon_cell_idx]][,,2]
         list_buoy_series[[d_idx]] <- series_temp
      } else {
# Extract stats.
         list_buoy_series[[d_idx]] <- list_data[[buoy_lat_cell_idx,buoy_lon_cell_idx]]
      }
   }

#=======================================================================#
# Load index data (ONI, NAO, PDO).
#=======================================================================#
# ONI.
   df_ONI <- read.csv("/home/ben/research/NOC/SRS_wave_analysis/datasets/indices/ONI.csv")
# Matrix containing mean, var and "min or max".
   mat_ONI <- matrix(NA,nrow=dim(df_ONI)[1],ncol=3)
   mat_ONI[,1] <- apply(X=df_ONI[,2:13],MAR=1,FUN=mean)
   mat_ONI[,2] <- apply(X=df_ONI[,2:13],MAR=1,FUN=var)
# Min or max based upon whether the mean is positive or negative.
   #for (i in 1:dim(df_ONI)[1]) {
   for (i in 1:69) {
      if (mat_ONI[i,1] < 0) {
         mat_ONI[i,3] <- min(df_ONI[i,2:13])
      } else {
         mat_ONI[i,3] <- max(df_ONI[i,2:13])
      }
   }
# Column and row names.
   colnames(mat_ONI) <- c("ONI_mean","ONI_var","ONI_minmax")
   rownames(mat_ONI) <- df_ONI[,1]
# ONI range based upon 2018.
   ONI_years <- which( as.numeric(rownames(mat_ONI)) == start_year ):which( as.numeric(rownames(mat_ONI)) == 2018 )

# NAO.
   df_NAO <- read.csv("/home/ben/research/NOC/SRS_wave_analysis/datasets/indices/norm.nao.monthly.b5001.current.ascii.table",header=FALSE, sep=',')
# Matrix containing mean, var and "min or max".
   mat_NAO <- matrix(NA,nrow=dim(df_NAO)[1],ncol=3)
   mat_NAO[,1] <- apply(X=df_NAO[,2:13],MAR=1,FUN=mean)
   mat_NAO[,2] <- apply(X=df_NAO[,2:13],MAR=1,FUN=var)
# Min or max based upon whether the mean is positive or negative.
   #for (i in 1:dim(df_NAO)[1]) {
   for (i in 1:69) {
      if (mat_NAO[i,1] < 0) {
         mat_NAO[i,3] <- min(df_NAO[i,2:13])
      } else {
         mat_NAO[i,3] <- max(df_NAO[i,2:13])
      }
   }
# Column and row names.
   colnames(mat_NAO) <- c("NAO_mean","NAO_var","NAO_minmax")
   rownames(mat_NAO) <- df_ONI[,1]
# NAO range based upon 2018.
   NAO_years <- which( as.numeric(rownames(mat_NAO)) == start_year ):which( as.numeric(rownames(mat_NAO)) == 2018 )

# PDO.
   df_PDO <- read.csv("/home/ben/research/NOC/SRS_wave_analysis/datasets/indices/PDO_NOAA.csv",header=TRUE, sep=',')
# Matrix containing mean, var and "min or max".
   mat_PDO_temp <- t(matrix(df_PDO$Value[-c(1981:1985)],nrow=12))

   mat_PDO <- matrix(NA,nrow=dim(mat_PDO_temp)[1],ncol=3)
   mat_PDO[,1] <- apply(X=mat_PDO_temp,MAR=1,FUN=mean)
   mat_PDO[,2] <- apply(X=mat_PDO_temp,MAR=1,FUN=var)
# Min or max based upon whether the mean is positive or negative.
   #for (i in 1:dim(df_PDO)[1]) {
   for (i in 1:dim(mat_PDO_temp)[1]) {
      if (mat_PDO[i,1] < 0) {
         mat_PDO[i,3] <- min(mat_PDO_temp[i,])
      } else {
         mat_PDO[i,3] <- max(mat_PDO_temp[i,])
      }
   }
# Column and row names.
   colnames(mat_PDO) <- c("PDO_mean","PDO_var","PDO_minmax")
   rownames(mat_PDO) <- 1854:2018
# PDO range based upon 2018.
   PDO_years <- which( as.numeric(rownames(mat_PDO)) == start_year ):which( as.numeric(rownames(mat_PDO)) == 2018 )

#=======================================================================#
# Data processing and QC for buoy data.
#=======================================================================#
# Data structures for monthly.
   mat_monthly_hs <- mat_monthly_hs1 <- matrix(list(),length(seq_years),12)
# Loop over years to get monthly raw data.
   for (yy in 1:length(seq_years)) {
      if ( seq_years[yy] <= 1998 ) {
         for (mm in 1:12) {
            #mat_monthly_hs[[yy,mm]] <- hist_buoy_hs[which(mat_buoy_obs[,1] == (seq_years[yy]-1900) & mat_buoy_obs[,2] %in% mm)]
            month_hs_idx <- which(mat_buoy_obs[,1] == (seq_years[yy]-1900) & mat_buoy_obs[,2] %in% mm)
            mat_monthly_hs[[yy,mm]] <- cbind( hist_buoy_hs[month_hs_idx],mat_buoy_obs[month_hs_idx,3] )
	 }
      } else {
         for (mm in 1:12) {
            month_hs_idx <- which(mat_buoy_obs[,1] == seq_years[yy] & mat_buoy_obs[,2] %in% mm)
            mat_monthly_hs[[yy,mm]] <- cbind( hist_buoy_hs[month_hs_idx],mat_buoy_obs[month_hs_idx,3] )
	 }
      }
   }

# Data structure for annual.
   list_annual_hs <- list_annual_hs1 <- list()

# Quality control monthly.
   mat_month_q_flag <- matrix(TRUE,length(seq_years),12)
# Quality control annual.
   year_q_flag <- logical(length(seq_years))
   year_q_flag[1:length(seq_years)] <- TRUE

# Loop over years for quality control.
   for (yy in 1:length(seq_years)) {
      for (mm in 1:12) {
         mat_day_hs <- mat_monthly_hs[[yy,mm]]
         #if ( (sum(!is.na(mat_day_hs[,1])) < 96) | (length(unique(mat_day_hs[,2])) < 13) ) {
# Require at least 13 days, each containing at least 13 measurements (typically one per hour).
# Assume equally weighted averages.
         if ( ! sum( sapply(X=unique(mat_day_hs[,2]),FUN=function(x) { sum(!is.na(mat_day_hs[mat_day_hs[,2]==x,1])) }) > QC_hour_thresh ) > QC_day_thresh ) {
            mat_month_q_flag[yy,mm] <- FALSE
         } else {
            mat_monthly_hs1[[yy,mm]] <- sapply(X=unique(mat_day_hs[,2]),FUN=function(x) { mean(mat_day_hs[mat_day_hs[,2]==x,1],na.rm=T) })
         }
      }
      #print(paste("Year:",seq_years[yy],"SUM:",sum(unlist(sapply(X=c(1,2,3,10,11,12),FUN=function(x) { mat_month_q_flag[[yy,x]] })))))
      if ( sum(unlist(sapply(X=1:12,FUN=function(x) { mat_month_q_flag[[yy,x]] }))) < QC_mth_thresh ) {
         year_q_flag[yy] <- FALSE
         list_annual_hs[[yy]] <- NA
         list_annual_hs1[[yy]] <- NA
      } else {
# Raw data.
         hs_temp <- unlist(sapply(X=1:12,FUN=function(x) { mat_monthly_hs[[yy,x]][,1] }))
         list_annual_hs[[yy]] <- hs_temp[!is.na(hs_temp)]
# Daily sampling.
         hs_temp1 <- unlist(sapply(X=1:12,FUN=function(x) { mat_monthly_hs1[[yy,x]] }))
         list_annual_hs1[[yy]] <- hs_temp1[!is.na(hs_temp1)]
      }
   }
   list_annual_hs <- list_annual_hs1

#=======================================================================#
# Annual / seasonal data for data buoy.
#=======================================================================#
# Find mean, quantiles and uncertainty.
   vec_q <- c(0.5,0.9,0.95)
   mat_b_mean <- sapply(X=list_annual_hs,FUN=mean,na.rm=T)
   mat_b_mean[is.nan(mat_b_mean)] <- NA
   mat_b_q <- cbind(t(sapply(X=list_annual_hs,FUN=quantile,probs=vec_q,na.rm=T)),mat_b_mean)
   colnames(mat_b_q) <- c(paste("Q",100*c(0.5,0.9,0.95),sep=""),"mean")
   df_Q_temp <- data.frame(cbind(year=(seq_years+0.5),mat_b_q,mat_ONI[ONI_years,],mat_NAO[NAO_years,],mat_PDO[PDO_years,]))
   rownames(df_Q_temp) <- seq_years+0.5
   df_Q <- df_Q_temp[year_q_flag,]
# Estimate CIs.
   array_q_CI_temp <- array(0,dim=c(length(seq_years),2,(length(vec_q)+1)))
   for (qq in 1:length(vec_q)) {
      array_q_CI_temp[,,qq] <- t(sapply(X=list_annual_hs,FUN=function(x) { sort(x)[quantile.CI(length(x),q=vec_q[qq])$Interval] }))
   }
# Mean and variance.
   mat_b_var <- sqrt(sapply(X=list_annual_hs,FUN=var,na.rm=T))
   mat_b_var <- cbind(mat_b_mean-mat_b_var,mat_b_mean+mat_b_var)
   array_q_CI_temp[,,4] <- mat_b_var
   array_q_CI <- array_q_CI_temp[year_q_flag,,]
# Weights for weighted least squares.
   vec_CI_w <- ( max(array_q_CI[,1,qq],na.rm=T)/array_q_CI[,1,qq] + max(array_q_CI[,2,qq],na.rm=T)/array_q_CI[,2,qq] ) / 2
   vec_CI_w[is.na(vec_CI_w)] <- 1

#=======================================================================#
# Create the figure.
#=======================================================================#
# Line colours.
   line_cols <- c("orangered2","magenta4","gold2","green4","springgreen3","cornflowerblue","blue4","")

# File generation.
   lab_av <- "24 hr mean"
   file_lab_av <- paste(strsplit(lab_av,split=' ')[[1]][1:2],collapse='')
   if ( flag_buoy ) { file_lab_prod <- paste("B",paste(vec_prod,collapse=''),sep='') } else { file_lab_prod <- paste(vec_prod,collapse='') }
   fig_file_name <- paste("./figures/buoy_trends/melbourne/",buoy_name,"_QA_annual_",file_lab_av,"_",lab_reg,"_",flag_QC_thresh,"_",file_lab_prod,".png",sep="")
   #X11()
   png(filename = fig_file_name, width = 3200, height = 1800)
   par(mfrow=c(1,1),oma=c(1.5,1.5,6,1),mar=c(7.0,7.0,5.0,3),mgp=c(5,2,0))

#-----------------------------------------------------------------------#
# Regression (lm & Sen slope).
#-----------------------------------------------------------------------#
   #for (qq in 1:length(vec_q)) {
   for (qq in 4) {
      if ( lab_reg == "PDO" ) {
         #lm_B <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q,weights=vec_CI_w)
         lm_B <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year + PDO_mean,data=df_Q)
         top_title <- paste("Trend in Q50, Q90, Q95 at NDBC",buoy_name," (",lab_av,")\nMultiple regression (Hs ~ year + ",lab_reg,")",sep='')
      } else if ( lab_reg == "NAO" ) {
         lm_B <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year + NAO_mean,data=df_Q)
         top_title <- paste("Trend in Q50, Q90, Q95 at NDBC",buoy_name," (",lab_av,")\nMultiple regression (Hs ~ year + ",lab_reg,")",sep='')
      } else {
         lm_B <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q)
         #sen_B <- zyp.sen(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q)
         sen_B <- zyp.sen(mean ~ year,data=df_Q)
# CCI regression.
         mat_grid <- list_buoy_series[[1]]
         mat_grid_plot <- data.frame(year=(as.numeric(rownames(mat_grid))+0.5),mean=mat_grid[,vec_var_col[1]])
         lm_C <- lm(mean ~ year,data=mat_grid_plot)
         sen_C <- zyp.sen(mean ~ year,data=mat_grid_plot)
# Young.
         mat_grid <- list_buoy_series[[2]]
         mat_grid_plot <- data.frame(year=(as.numeric(rownames(mat_grid))+0.5),mean=mat_grid[,vec_var_col[2]])
         lm_Y <- lm(mean ~ year,data=mat_grid_plot)
         sen_Y <- zyp.sen(mean ~ year,data=mat_grid_plot)
# ERA5.
         mat_grid <- list_buoy_series[[3]]
         mat_grid_plot <- data.frame(year=(as.numeric(rownames(mat_grid))+0.5),mean=mat_grid[,vec_var_col[3]])
         lm_E <- lm(mean ~ year,data=mat_grid_plot)
         sen_E <- zyp.sen(mean ~ year,data=mat_grid_plot)
# ERA5 [1995 - 2005].
         mat_grid <- list_buoy_series[[3]]
         mat_grid_plot <- data.frame(year=(as.numeric(rownames(mat_grid))+0.5),mean=mat_grid[,vec_var_col[3]])
         lm_Es <- lm(mean ~ year,data=mat_grid_plot[14:39,])
         sen_Es <- zyp.sen(mean ~ year,data=mat_grid_plot[14:39,])

         #top_title <- paste("Trend in Q50, Q90, Q95 at NDBC",buoy_name," (",lab_av,")\nLinear regression (Hs ~ year)",sep='')
         top_title <- paste("Trend in mean (annual) at NDBC",buoy_name," (",lab_av,", ",lab_QC,")\nSen regression (Hs ~ year)",sep='')
      }
      sum_lm_B <- summary(lm_B)
      print(summary(lm_B))
#-----------------------------------------------------------------------#
# Find plot range.
#-----------------------------------------------------------------------#
      grid_range <- range(df_Q[,(qq+1)])
      for (d_idx in vec_prod) {
         mat_grid <- list_buoy_series[[d_idx]]
         grid_range[1] <- min(grid_range[1],min(mat_grid[,vec_var_col[d_idx]],na.rm=T),na.rm=T)
         grid_range[2] <- max(grid_range[2],max(mat_grid[,vec_var_col[d_idx]],na.rm=T),na.rm=T)
      }
      y_lim <- c( floor(grid_range[1]/0.25)*0.25, ceiling(grid_range[2]/0.25)*0.25 )
      #if (flag_buoy_monthly) {
      #   y_lim <- y_lim + c(-1.0,1.0)
      #}
      y_lim <- c(1.25,3.75)
#-----------------------------------------------------------------------#
# Plot buoy and axes.
#-----------------------------------------------------------------------#
      mat_buoy_plot <- df_Q[,c(1,(qq+1))]
      if ( flag_buoy ) {
         plot(mat_buoy_plot,xlim=c(1980,2020),ylim=y_lim,cex.main=3.0,cex.lab=3.0,xlab="Year",ylab="Mean (m)",pch=19,cex=4,axes=FALSE)
      } else {
         plot(NULL,xlim=c(1980,2020),ylim=y_lim,cex.main=3.0,cex.lab=3.0,xlab="Year",ylab="Mean (m)",pch=19,cex=4,axes=FALSE)
      }
      axis(side=1,at=seq(1980,2020,5),cex.lab=3.0,cex.axis=3.0,)
      axis(side=2,at=seq(0,7.5,0.25),cex.lab=3.0,cex.axis=3.0)
#-----------------------------------------------------------------------#
# Plot monthly means.
#-----------------------------------------------------------------------#
      if ( flag_buoy_monthly ) {
         seq_month <- seq(start_year+1/24,2019-1/24,1/12)
         vec_month_temp <- unlist(suppressWarnings(lapply(X=t(mat_monthly_hs1),FUN=mean,na.rm=T)))
         ## Regression on raw.
         #df_month_temp <- data.frame(time=seq(start_year+1/24,2019-1/24,1/12),mean=vec_month_temp)
         #df_month_temp <- df_month_temp[!is.na(df_month_temp[,2]),]
         #sen_mth1 <- zyp.sen(mean ~ time,data=df_month_temp)
         #lines(cbind(sen_mth1$x,(sen_mth1$coefficients[2]*sen_mth1$x+sen_mth1$coefficients[1])),lty=2,lwd=10,col="red")
         ## Regression on QC.
         vec_month_temp[!t(mat_month_q_flag)] <- NA
         #df_month_temp <- data.frame(time=seq(start_year+1/24,2019-1/24,1/12),mean=vec_month_temp)
         #df_month_temp <- df_month_temp[!is.na(df_month_temp[,2]),]
         #sen_mth1 <- zyp.sen(mean ~ time,data=df_month_temp)
         #lines(cbind(sen_mth1$x,(sen_mth1$coefficients[2]*sen_mth1$x+sen_mth1$coefficients[1])),lty=3,lwd=10,col="blue")
	 # Trace.
         points(cbind(seq_month,vec_month_temp),pch=19,cex=2,col="grey30")
         lines(cbind(seq_month,vec_month_temp),lwd=2,col="grey30")
      }
#-----------------------------------------------------------------------#
# Plot background grid lines.
#-----------------------------------------------------------------------#
      abline(v=seq(1980,2020,5),col="darkgrey",lwd=3,lty=3)
      abline(h=seq(0,7.5,0.25),col="darkgrey",lwd=3,lty=3)
#-----------------------------------------------------------------------#
# Plot gridded products.
#-----------------------------------------------------------------------#
      for (d_idx in vec_prod) {
         mat_grid <- list_buoy_series[[d_idx]]
         mat_grid_plot <- cbind((as.numeric(rownames(mat_grid))+0.5),mat_grid[,vec_var_col[d_idx]])
         #points(mat_grid_plot,pch=(d_idx+1),lwd=4,cex=6,col=line_cols[d_idx])
         lines(mat_grid_plot,lwd=9,col=line_cols[d_idx])
      }
# Buoy lines, uncertainty and trend.
#      buoy_seg_plot <- (mat_buoy_plot[2:dim(mat_buoy_plot)[1],1] - mat_buoy_plot[(1:(dim(mat_buoy_plot)[1])-1),1]) == 1
#      for (y_idx in 2:dim(mat_buoy_plot)[1]) {
#         if (buoy_seg_plot[y_idx-1]) {
#            lines(mat_buoy_plot[c((y_idx-1):y_idx),],lwd=6)
#         }
#      }

      #lines(df_Q[,c(1,(qq+1))],lwd=6)
      #arrows(df_Q$year, array_q_CI[,1,qq], df_Q$year, array_q_CI[,2,qq], length=0.05, angle=90, code=3, lwd=2)
#-----------------------------------------------------------------------#
# Plot trends.
#-----------------------------------------------------------------------#
      if ( flag_buoy ) {
# Buoy.
# Linear trend.
         newdata <- data.frame(year=as.numeric(names(lm_B$fitted.values)))
         pred_lm_B <- predict(lm_B,newdata)
         lines(cbind(newdata,pred_lm_B),lty=5,lwd=6)
# Sen slope.
         lines(cbind(sen_B$x,(sen_B$coefficients[2]*sen_B$x+sen_B$coefficients[1])),lty=1,lwd=6)
      }
      if ( flag_reg_lines ) {
# Young.
#         newdata <- data.frame(year=as.numeric(names(lm_Y$fitted.values)))
#         pred_lm_Y <- predict(lm_Y,newdata)
#         lines(cbind(newdata,pred_lm_Y),col=line_cols[2],lty=2,lwd=9)
# Sen slope.
         lines(cbind(sen_Y$x,(sen_Y$coefficients[2]*sen_Y$x+sen_Y$coefficients[1])),col=line_cols[2],lty=1,lwd=6)
# ERA5.
         #newdata <- data.frame(year=lm_E$model$year)
         #pred_lm_E <- predict(lm_E,newdata)
         #lines(cbind(newdata,pred_lm_E),col=line_cols[3],lty=1,lwd=6)
# Sen slope.
         lines(cbind(sen_E$x,(sen_E$coefficients[2]*sen_E$x+sen_E$coefficients[1])),col=line_cols[3],lty=1,lwd=6)
# ERA 1995 - 2015.
#         newdata <- data.frame(year=lm_Es$model$year)
#         pred_lm_Es <- predict(lm_Es,newdata)
#         lines(cbind(newdata,pred_lm_Es),col=line_cols[3],lty=2,lwd=9)
# Sen slope.
         lines(cbind(sen_E$x,(sen_E$coefficients[2]*sen_E$x+sen_E$coefficients[1])),col="black",lty=1,lwd=6)
         lines(cbind(sen_E$x,(sen_E$coefficients[2]*sen_E$x+sen_E$coefficients[1])),col=line_cols[3],lty=5,lwd=6)
         lines(cbind(sen_Es$x,(sen_Es$coefficients[2]*sen_Es$x+sen_Es$coefficients[1])),col="black",lty=1,lwd=6)
         lines(cbind(sen_Es$x,(sen_Es$coefficients[2]*sen_Es$x+sen_Es$coefficients[1])),col=line_cols[3],lty=5,lwd=6)
# CCI.
         #newdata <- data.frame(year=lm_C$model$year)
         #pred_lm_C <- predict(lm_C,newdata)
         #lines(cbind(newdata,pred_lm_C),col=line_cols[1],lty=1,lwd=6)
# Sen slope.
         lines(cbind(sen_C$x,(sen_C$coefficients[2]*sen_C$x+sen_C$coefficients[1])),col=line_cols[1],lty=1,lwd=6)
      }
#-----------------------------------------------------------------------#
# Legend.
#-----------------------------------------------------------------------#
      if ( lab_reg == "PDO" ) {
         legend(x=2000,y=3.8,legend=c("Reg. line","PDO","ONI"),lty=c(1,1,1),lwd=c(5,4,4),col=c("red","darkblue","darkorange"),cex=3)
      } else if ( flag_buoy ) {
         legend(x=1980,y=y_lim[2],legend=c("Buoy","Lin. reg. line (buoy)","Sen reg. line (buoy)",plot_labels[vec_prod]),lty=c(NA,5,1,1,1,1,1,1,2,2,3),lwd=c(5),pch=c(19,rep(NA,10)),col=c("black","black","black",line_cols[vec_prod],line_cols[vec_prod[c(1:3,3)]]),cex=3)
         #legend(x=1980,y=y_lim[2],legend=c("Buoy","Lin. reg. line (Annual QC)","Sen reg. line (Annual QC)","Sen reg. line (Monthly ALL)","Sen reg. line (Monthly QC)"),lty=c(NA,5,1,2,3),lwd=c(5),pch=c(19,rep(NA,10)),col=c("black","black","black","red","blue"),cex=3)
      } else if ( !flag_buoy ) {
         legend(x=1980,y=y_lim[2],legend=c(plot_labels[vec_prod]),lty=c(1,1,1,1,1,1,2,2,3),lwd=c(5),pch=c(rep(NA,10)),col=c(line_cols[vec_prod],line_cols[vec_prod[c(1:3,3)]]),cex=3)
      }
#-----------------------------------------------------------------------#
# Plot index for multiple regression.
#-----------------------------------------------------------------------#
      if ( lab_reg == "PDO" ) {
         par(new=TRUE)
         plot(df_Q$year,df_Q$PDO_mean,ylim=c(-2,9),xlab="",ylab="",col="darkblue",axes=FALSE)
	 axis(side=4,at=c(-2:2),cex.axis=2.0)
         lines(df_Q$year,df_Q$PDO_mean,col="darkblue",lwd=4)
         par(new=TRUE)
         plot(df_Q$year,df_Q$ONI_mean,ylim=c(-2,9),xlab="",ylab="",col="darkorange",axes=FALSE)
         lines(df_Q$year,df_Q$ONI_mean,col="darkorange",lwd=4)
      #} else {
      #   par(new=TRUE)
      #   plot(df_Q$year,df_Q$NAO_mean,ylim=c(-4,4),xlab="",ylab="",axes=FALSE)
      #   lines(df_Q$year,df_Q$NAO_mean)
      }

      if ( flag_buoy ) {
         if ( sum_lm_B$coefficients[2,4] < 0.1 ) {
            trend_sig <- paste("Buoy data, year (trend) coef.: ",format(sum_lm_B$coefficients[2,1],digits=2),"m per year.\nSignificant at 10%, Pr(>|t|) = ",format(sum_lm_B$coefficients[2,4],digits=2),sep="")
         } else {
            trend_sig <- paste("Buoy data, year (trend) coef.: ",format(sum_lm_B$coefficients[2,1],digits=2),"m per year.\nNot significant, Pr(>|t|) = ",format(sum_lm_B$coefficients[2,4],digits=2),sep="")
         }
         mtext(text = trend_sig, side = 3, line = -8, cex = 3)
      }
   }
   #mtext(text = paste("Trend in Q50, Q90, Q95 at NDBC",buoy_name," (",lab_av,")\nMultiple regression (year,",lab_reg,")",sep=''), outer = TRUE, side = 3, line = -3, cex = 3.5)
   mtext(text = top_title, outer = TRUE, side = 3, line = -3, cex = 3.5)
   dev.off()

   system(paste("okular",fig_file_name,"&> /dev/null &"))
#   }
}
   

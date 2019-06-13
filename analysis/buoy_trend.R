# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/buoy_trend.R")

# Script to load some SRS data from different satellites and compare with buoy NDBC 46066.
# BT 04/2019

# Libraries.
   library(ncdf4)
   library(TTR)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")

# File base name.
   #IMOS_SRS-Surface-Waves_MW_JASON-1_FV02_052N-205E-DM00.nc

# Buoy 46066, 52.785N 155.047W
   buoy_list <- c("46066")
   #buoy_list <- c("51004")
   #buoy_list <- c("41002")
   #buoy_list <- c("46006")

# Averaging flag.
   #flag_av <- "24hr"
   flag_av <- "6hr"

# Regression covariate.
   lab_reg <- "NAO"

#-----------------------------------------------------------------------#
# Buoy data.
#-----------------------------------------------------------------------#

# Load historical data for NDBC buoys.
   #array_buoy_obs <- array(0,dim=c(100000,16,length(buoy_list)))
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
# Load index data (ONI, NAO).
#-----------------------------------------------------------------------#
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

#-----------------------------------------------------------------------#
# Trend based upon annual data.
#-----------------------------------------------------------------------#
   list_years <- list(length(seq_years))
   for (yy in 1:length(seq_years)) {
      if ( seq_years[yy] <= 1998 ) {
         #hs_temp <- hist_buoy_hs[which(mat_buoy_obs[,1] == (seq_years[yy]-1900))]
         #list_years[[yy]] <- hs_temp[!is.na(hs_temp)]
         list_years[[yy]] <- hist_buoy_hs[which(mat_buoy_obs[,1] == (seq_years[yy]-1900))]
      } else {
         #hs_temp <- hist_buoy_hs[which(mat_buoy_obs[,1] == seq_years[yy])]
         #list_years[[yy]] <- hs_temp[!is.na(hs_temp)]
         list_years[[yy]] <- hist_buoy_hs[which(mat_buoy_obs[,1] == seq_years[yy])]
      }
   }

# Take fixed daily mean.
   #hs_temp_resid <- NULL
   year_q_flag <- logical(length(seq_years))
   year_q_flag[1:length(seq_years)] <- TRUE

# Averaging.
   if (flag_av == "24hr") {
# 24 hour averaging.
      mat_day_idx <- cbind(seq(1,,24,400),seq(24,,24,400))
      mat_day_idx_2016 <- cbind(seq(1,,144,2400),seq(144,,144,2400))
      qual_thresh <- 274
      lab_av <- "24 hr mean"
   } else if (flag_av == "6hr") {
# 6 hour averaging.
      mat_day_idx <- cbind(seq(1,,6,1600),seq(6,,6,1600))
      mat_day_idx_2016 <- cbind(seq(1,,36,2400),seq(36,,36,2400))
      qual_thresh <- 4*274
      lab_av <- "6 hr mean"
   }

   for (yy in 1:length(seq_years)) {
      if (buoy_name == 46006 & seq_years[yy] >= 2016) {
         hs_temp <- apply(X=mat_day_idx_2016,MAR=1,FUN=function(x) { mean(list_years[[yy]][x[1]:x[2]],na.rm=T) })
         list_years[[yy]] <- hs_temp[!is.nan(hs_temp)]
         if (length(list_years[[yy]]) < qual_thresh) {
            year_q_flag[yy] <- FALSE
         }
      } else if (buoy_name == 41002 & seq_years[yy] >= 2018) {
         hs_temp <- apply(X=mat_day_idx_2016,MAR=1,FUN=function(x) { mean(list_years[[yy]][x[1]:x[2]],na.rm=T) })
         list_years[[yy]] <- hs_temp[!is.nan(hs_temp)]
         if (length(list_years[[yy]]) < qual_thresh) {
            year_q_flag[yy] <- FALSE
         }
      } else if (buoy_name == 51004 & seq_years[yy] >= 2018) {
         hs_temp <- apply(X=mat_day_idx_2016,MAR=1,FUN=function(x) { mean(list_years[[yy]][x[1]:x[2]],na.rm=T) })
         list_years[[yy]] <- hs_temp[!is.nan(hs_temp)]
         if (length(list_years[[yy]]) < qual_thresh) {
            year_q_flag[yy] <- FALSE
         }
      } else {
         hs_temp <- apply(X=mat_day_idx,MAR=1,FUN=function(x) { mean(list_years[[yy]][x[1]:x[2]],na.rm=T) })
         #hs_temp_resid_temp <- apply(X=mat_day_idx,MAR=1,FUN=function(x) { list_years[[yy]][x[1]:x[2]] - median(list_years[[yy]][x[1]:x[2]]) })
         #hs_temp_resid <- c(hs_temp_resid,hs_temp_resid_temp[!is.na(hs_temp_resid_temp)])
         list_years[[yy]] <- hs_temp[!is.na(hs_temp)]
         if (length(list_years[[yy]]) < qual_thresh) {
            year_q_flag[yy] <- FALSE
         }
      }
   }

## Take moving daily average.
#   for (yy in 1:length(seq_years)) {
#      if ( length(list_years[[yy]]) > 100 ) {
#         hs_temp <- SMA(x=list_years[[yy]],n=48)
#         list_years[[yy]] <- hs_temp[!is.na(hs_temp)]
#      } else {
#         list_years[[yy]] <- NA
#      }
#   }

# Find quantiles and uncertainty.
   vec_q <- c(0.5,0.9,0.95)
   mat_b_q <- t(sapply(X=list_years,FUN=quantile,probs=vec_q,na.rm=T))
   colnames(mat_b_q) <- paste("Q",100*c(0.5,0.9,0.95),sep="")
   df_Q_temp <- data.frame(cbind(year=seq_years,mat_b_q,mat_ONI[ONI_years,],mat_NAO[NAO_years,],mat_PDO[PDO_years,]))
   df_Q <- df_Q_temp[year_q_flag,]
# Estimate CIs.
   array_q_CI_temp <- array(0,dim=c(length(seq_years),2,length(vec_q)))
   for (qq in 1:length(vec_q)) {
      array_q_CI_temp[,,qq] <- t(sapply(X=list_years,FUN=function(x) { sort(x)[quantile.CI(length(x),q=vec_q[qq])$Interval] }))
   }
   array_q_CI <- array_q_CI_temp[year_q_flag,,]
# Weights for weighted least squares.
   vec_CI_w <- ( max(array_q_CI[,1,qq],na.rm=T)/array_q_CI[,1,qq] + max(array_q_CI[,2,qq],na.rm=T)/array_q_CI[,2,qq] ) / 2
   vec_CI_w[is.na(vec_CI_w)] <- 1

   file_lab_av <- paste(strsplit(lab_av,split=' ')[[1]][1:2],collapse='')
   fig_file_name <- paste("./figures/buoy_trends/",buoy_name,"_QA_annual_",file_lab_av,"_",lab_reg,".png",sep="")
   #X11()
   png(filename = fig_file_name, width = 2200, height = 2200)
   par(mfrow=c(2,2),oma=c(1.5,1.5,3,1),mar=c(7.0,7.0,5.0,3),mgp=c(5,2,0))

   for (qq in 1:length(vec_q)) {
      if ( lab_reg == "PDO" ) {
         #lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q,weights=vec_CI_w)
         lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year + PDO_mean,data=df_Q)
      } else if ( lab_reg == "NAO" ) {
         lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year + NAO_mean,data=df_Q)
      } else {
         lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q)
      }
      sum_lm_Q <- summary(lm_Q)
      print(summary(lm_Q))
# Plot data, error bars and trend.
      plot(df_Q[,c(1,(qq+1))],ylim=c(0,7),cex.main=3.0,cex.lab=3.0,cex.axis=3.0,lwd=3.0)
      arrows(df_Q$year, array_q_CI[,1,qq], df_Q$year, array_q_CI[,2,qq], length=0.05, angle=90, code=3, lwd=2)
      abline(lm_Q,lwd=3)
# Plot index.
      if ( lab_reg == "PDO" ) {
         par(new=TRUE)
         plot(df_Q$year,df_Q$PDO_mean,ylim=c(-4,4),xlab="",ylab="",axes=FALSE)
         lines(df_Q$year,df_Q$PDO_mean)
      } else {
         par(new=TRUE)
         plot(df_Q$year,df_Q$NAO_mean,ylim=c(-4,4),xlab="",ylab="",axes=FALSE)
         lines(df_Q$year,df_Q$NAO_mean)
      }

      if ( sum_lm_Q$coefficients[2,4] < 0.1 ) {
         trend_sig <- paste("Trend: ",format(sum_lm_Q$coefficients[2,1],digits=2),"m per year.\nSignificant at 10%, Pr(>|t|) = ",format(sum_lm_Q$coefficients[2,4],digits=2),sep="")
      } else {
         trend_sig <- paste("Trend: ",format(sum_lm_Q$coefficients[2,1],digits=2),"m per year.\nNot significant, Pr(>|t|) = ",format(sum_lm_Q$coefficients[2,4],digits=2),sep="")
      }
      mtext(text = trend_sig, side = 3, line = -8, cex = 3)
   }
   mtext(text = paste("Trend in Q50, Q90, Q95 at NDBC",buoy_name," (",lab_av,")",sep=''), outer = TRUE, side = 3, line = -2, cex = 4)
   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))
   
##-----------------------------------------------------------------------#
## Satellite data.
##-----------------------------------------------------------------------#
## Start figures.
#   #X11()
#   png(filename = "./figures/46066_5678_annual.png",width=2000,height=2000)
#   par(mfrow=c(2,2),oma=c(1.5,1.5,2,1),mar=c(6.0,7.0,5.0,3),mgp=c(5,2,0))
#
## Loop over different data sets.
#   SHW_KU <- NULL
#
#   #for (k in 1:length(vec_datasets) {
#   for (k in c(5,6,7,8)) {
#      #vec_lat_sat1 <- paste("0",50:55,"N",sep="")
#      #vec_lon_sat1 <- paste(200:210,"E",sep="")
#      nc1_filename <- paste(vec_datasets[k],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[k],"_FV02_052N-205E-DM00.nc",sep="")
#
#      nc1 = nc_open(paste(data_path,nc1_filename,sep=""))
#   
## SWH.
#      nc1_SHW_KU <- ncvar_get(nc1,"SWH_KU")
#      SHW_KU <- c(SHW_KU,nc1_SHW_KU)
#      nc1_SHW_C <- ncvar_get(nc1,"SWH_C")
#
## Load data.
## Time.
#      nc1_time_idx <- ncvar_get(nc1,"TIME")
#      nc1_date <- as.POSIXct(nc1_time_idx*3600*24, origin = '1985-01-01', tz='GMT')
#      nc1_start <- nc1_date[1]
#
#      nc1_lon <- ncvar_get(nc1,"LONGITUDE")
#      nc1_lat <- ncvar_get(nc1,"LATITUDE")
#
## Look for breaks.
#      nc1_breaks <- NULL
#      for (i in 2:length(nc1_time_idx)) {
#         AA <- nc1_time_idx[i] - nc1_time_idx[i-1]
#         if ( AA > 1e-4 ) {
#            nc1_breaks <- c(nc1_breaks,i)
#            #print(paste(" Break before:",i))
#         }
#      }
#   
## Matrix of 'block' indices for contiguous readings.
#      mat_nc1_breaks <- cbind(c(1,nc1_breaks),c(nc1_breaks-1,length(nc1_time_idx)))
#      vec_nc1_break_mid <- apply(X=mat_nc1_breaks,MAR=1,FUN=function(x) { floor(mean(x[1]:x[2])) })
#
## ================================================================= #
## Analysis:
## Wave height in each block. Matrix dim should not be the size of the entire record!
#      mat_nc1_block_first <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
#      mat_nc1_block_max <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
#      mat_nc1_block_mean <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
#      nc1_block_cor <- 0
#
#      for (i in 1:dim(mat_nc1_breaks)[1]) {
#         CC <- mat_nc1_breaks[i,1]:mat_nc1_breaks[i,2]
#         mat_nc1_block_first[i,] <- c( nc1_SHW_C[CC[1]],nc1_SHW_KU[CC[1]] )
#         mat_nc1_block_max[i,] <- c( max(nc1_SHW_C[CC]),max(nc1_SHW_KU[CC]) )
#         mat_nc1_block_mean[i,] <- c( mean(nc1_SHW_C[CC]),mean(nc1_SHW_KU[CC]) )
#         nc1_block_cor[i] <- cor(nc1_SHW_C[CC],nc1_SHW_KU[CC])
#      }
## Seasonal info.
#      nc1_month <- sapply(X=nc1_date[vec_nc1_break_mid],FUN=substr,start=6,stop=7)
#      nc1_ONDJFM <- which(nc1_month == "01" | nc1_month == "02" | nc1_month == "03" | nc1_month == "10" | nc1_month == "11" | nc1_month == "12")
#
#      #cor(mat_nc1_block_mean)
#
## Close files.
#      nc_close(nc1)
#   #}
#
#   #X11()
#   #qqplot(hist_buoy_hs,mat_nc1_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main="Comparison at NDBC 46066")
#   func_qqplot(hist_buoy_hs,SHW_KU,
#          #xlim=c(0,12),ylim=c(0,12),
#          xlim=c(0,4),ylim=c(0,4),regress = FALSE,
#          xlab="Buoy",ylab="Satellite",
#          main=paste(vec_datasets[k],": QQ at NDBC 46066",sep=""),cex.main=3.0,cex.lab=3.0,cex.axis=3.0,lwd.reg=5.0)
## Quantiles.
#   for (kk in 1:2) {
#      segments(x0=q_plot[kk], y0=0, x1=q_plot[kk], y1=q_plot[kk], lty=3, lwd=3, col="red")
#      segments(x0=0, y0=q_plot[kk], x1=q_plot[kk], y1=q_plot[kk], lty=3, lwd=3, col="red")
#   }
## Seasonal.
#   #buoy_month <- mat_buoy_obs[,2]
#   #buoy_ONDJFM <- which(buoy_month == "1" | buoy_month == "2" | buoy_month == "3" | buoy_month == "10" | buoy_month == "11" | buoy_month == "12")
#   #qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc1_block_mean[nc1_ONDJFM,],xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main="Seasonal (ONDJFM) at NDBC 46066")
#
#   }
#
#   #X11()
#   #qqplot(hist_buoy_hs,mat_nc2_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main="Comparison at NDBC 46066")
## Seasonal.
#   #qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc2_block_mean[nc2_ONDJFM,],xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main="Seasonal (ONDJFM) at NDBC 46066")
#   #abline(a=0,b=1)
#
#   dev.off()
#

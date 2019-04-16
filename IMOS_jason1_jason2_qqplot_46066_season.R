# source("/home/ben/research/NOC/waves/IMOS_jason1_jason2_qqplot_46066_season.R")

# Script to load some IMOS data from two different satellites and poke around.
# BT 03/2019

# Libraries.
   library(ncdf4)

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/waves/datasets/IMOS/"

# File base name.
   #IMOS_SRS-Surface-Waves_MW_JASON-1_FV02_052N-205E-DM00.nc
   nc_sat1 <- "IMOS_SRS-Surface-Waves_MW_JASON-1_FV02_"
   vec_lat_sat1 <- paste("0",50:55,"N",sep="")
   vec_lon_sat1 <- paste(200:210,"E",sep="")
   vec_lonlat_sat1 <- paste(rep(vec_lat_sat1,times=length(vec_lon_sat1)),"-",rep(vec_lon_sat1,each=length(vec_lat_sat1)),sep="")
   vec_filenames_sat1 <- paste(nc_sat1,vec_lonlat_sat1,"-DM00.nc",sep="")

   nc_sat2 <- "IMOS_SRS-Surface-Waves_MW_JASON-2_FV02_"
   vec_lat_sat2 <- paste("0",50:55,"N",sep="")
   vec_lon_sat2 <- paste(200:210,"E",sep="")
   vec_lonlat_sat2 <- paste(rep(vec_lat_sat2,times=length(vec_lon_sat2)),"-",rep(vec_lon_sat2,each=length(vec_lat_sat2)),sep="")
   vec_filenames_sat2 <- paste(nc_sat2,vec_lonlat_sat2,"-DM00.nc",sep="")

# Time slice.
   time_lab <- "winter"
   file_name <- paste("./figures/J1_J2_46066_",time_lab,".pdf",sep="")

# ================================================================= #
# Loop over different locations.

# Which grid cell.
   file_list <- c(26,32,38,27,33,39,28,34,40)
   #file_list <- c(33)

# Variables.
   mat_nc1_plot <- NULL
   mat_nc2_plot <- NULL

   pdf(file = file_name,width=10,height=3)
   #if ( length(file_list) == 9) {
   #   par(mfrow=c(3,3))
   #} else {
      par(mfrow=c(1,3))
   #}

   for (k in 1:length(file_list)) {
# Buoy 46066.
      nc1 = nc_open(paste(data_path,vec_filenames_sat1[file_list[k]],sep=""))
      nc2 = nc_open(paste(data_path,vec_filenames_sat2[file_list[k]],sep=""))
   
# SWH.
      nc1_SHW_KU <- ncvar_get(nc1,"SWH_KU")
      nc1_SHW_C <- ncvar_get(nc1,"SWH_C")
      nc2_SHW_KU <- ncvar_get(nc2,"SWH_KU")
      nc2_SHW_C <- ncvar_get(nc2,"SWH_C")

# Load data.
# Time.
      nc1_time_idx <- ncvar_get(nc1,"TIME")
      nc1_date <- as.POSIXct(nc1_time_idx*3600*24, origin = '1985-01-01', tz='GMT')
      nc1_start <- nc1_date[1]

      nc1_lon <- ncvar_get(nc1,"LONGITUDE")
      nc1_lat <- ncvar_get(nc1,"LATITUDE")

# Look for breaks.
      nc1_breaks <- NULL
      for (i in 2:length(nc1_time_idx)) {
         AA <- nc1_time_idx[i] - nc1_time_idx[i-1]
         if ( AA > 1e-4 ) {
            nc1_breaks <- c(nc1_breaks,i)
            #print(paste(" Break before:",i))
         }
      }
   
# Matrix of 'block' indices for contiguous readings.
      mat_nc1_breaks <- cbind(c(1,nc1_breaks),c(nc1_breaks-1,length(nc1_time_idx)))
      vec_nc1_break_mid <- apply(X=mat_nc1_breaks,MAR=1,FUN=function(x) { floor(mean(x[1]:x[2])) })

# ----------------------------------------------------------------- #

# Load data.
# Time.
      nc2_time_idx <- ncvar_get(nc2,"TIME")
      nc2_date <- as.POSIXct(nc2_time_idx*3600*24, origin = '1985-01-01', tz='GMT')

      nc2_lon <- ncvar_get(nc2,"LONGITUDE")
      nc2_lat <- ncvar_get(nc2,"LATITUDE")

# Look for breaks.
      nc2_breaks <- NULL
      for (i in 2:length(nc2_time_idx)) {
         AA <- nc2_time_idx[i] - nc2_time_idx[i-1]
         if ( AA > 1e-4 ) {
            nc2_breaks <- c(nc2_breaks,i)
            #print(paste(" Break before:",i))
         }
      }
   
# Matrix of 'block' indices for contiguous readings.
      mat_nc2_breaks <- cbind(c(1,nc2_breaks),c(nc2_breaks-1,length(nc2_time_idx)))
      vec_nc2_break_mid <- apply(X=mat_nc2_breaks,MAR=1,FUN=function(x) { floor(mean(x[1]:x[2])) })

# ================================================================= #
# Analysis:
# Wave height in each block. Matrix dim should not be the size of the entire record!
      mat_nc1_block_first <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      mat_nc1_block_max <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      mat_nc1_block_mean <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      nc1_block_cor <- 0

      for (i in 1:dim(mat_nc1_breaks)[1]) {
         CC <- mat_nc1_breaks[i,1]:mat_nc1_breaks[i,2]
         mat_nc1_block_first[i,] <- c( nc1_SHW_C[CC[1]],nc1_SHW_KU[CC[1]] )
         mat_nc1_block_max[i,] <- c( max(nc1_SHW_C[CC]),max(nc1_SHW_KU[CC]) )
         mat_nc1_block_mean[i,] <- c( mean(nc1_SHW_C[CC]),mean(nc1_SHW_KU[CC]) )
         nc1_block_cor[i] <- cor(nc1_SHW_C[CC],nc1_SHW_KU[CC])
      }
# Seasonal info.
      nc1_month <- sapply(X=nc1_date[vec_nc1_break_mid],FUN=substr,start=6,stop=7)
      nc1_ONDJFM <- which(nc1_month == "01" | nc1_month == "02" | nc1_month == "03" | nc1_month == "10" | nc1_month == "11" | nc1_month == "12")
      #nc1_ONDJFM <- which(nc1_month == "04" | nc1_month == "05" | nc1_month == "06" | nc1_month == "07" | nc1_month == "08" | nc1_month == "09")
      #nc1_ONDJFM <- which(nc1_month == "01")

# Aggregate over grid cells.
      mat_nc1_plot <- c( mat_nc1_plot,as.vector(mat_nc1_block_mean[nc1_ONDJFM,]) )

# Wave height in each block.
      mat_nc2_block_first <- matrix(0,nrow=dim(mat_nc2_breaks)[1],ncol=2)
      mat_nc2_block_max <- matrix(0,nrow=dim(mat_nc2_breaks)[1],ncol=2)
      mat_nc2_block_mean <- matrix(0,nrow=dim(mat_nc2_breaks)[1],ncol=2)
      nc2_block_cor <- 0

      for (i in 1:dim(mat_nc2_breaks)[1]) {
         CC <- mat_nc2_breaks[i,1]:mat_nc2_breaks[i,2]
         mat_nc2_block_first[i,] <- c( nc2_SHW_C[CC[1]],nc2_SHW_KU[CC[1]] )
         mat_nc2_block_max[i,] <- c( max(nc2_SHW_C[CC]),max(nc2_SHW_KU[CC]) )
         mat_nc2_block_mean[i,] <- c( mean(nc2_SHW_C[CC]),mean(nc2_SHW_KU[CC]) )
         nc2_block_cor[i] <- cor(nc2_SHW_C[CC],nc2_SHW_KU[CC])
      }
# Seasonal info.
      nc2_month <- sapply(X=nc2_date[vec_nc2_break_mid],FUN=substr,start=6,stop=7)
      nc2_ONDJFM <- which(nc2_month == "01" | nc2_month == "02" | nc2_month == "03" | nc2_month == "10" | nc2_month == "11" | nc2_month == "12")
      #nc2_ONDJFM <- which(nc2_month == "04" | nc2_month == "05" | nc2_month == "06" | nc2_month == "07" | nc2_month == "08" | nc2_month == "09")
      #nc2_ONDJFM <- which(nc2_month == "01")

# Aggregate over grid cells.
      mat_nc2_plot <- c( mat_nc2_plot,as.vector(mat_nc2_block_mean[nc2_ONDJFM,]) )

# Close files.
      nc_close(nc1)
      nc_close(nc2)
   }

# Seasonal.
      qqplot(mat_nc1_plot,mat_nc2_plot,xlim=c(0,12),ylim=c(0,12),
             xlab="Jason 1",
             ylab="Jason 2",
             main=paste("Hs: J 1 vs J 2 at NDBC 46066 (",time_lab,")",sep=""))
      abline(a=0,b=1)

#-----------------------------------------------------------------------#
# Observed data.
#-----------------------------------------------------------------------#

   buoy_list <- c("46066")

# Load historical data for NDBC buoys.
   #array_buoy_obs <- array(0,dim=c(100000,16,length(buoy_list)))
   mat_buoy_obs <- matrix(0,nrow=100000,ncol=16)
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

   #X11()
   #qqplot(hist_buoy_hs,mat_nc1_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main="Comparison at NDBC 46066")
# Seasonal.
   buoy_month <- mat_buoy_obs[,2]
   #buoy_ONDJFM <- which(buoy_month == "1")
   buoy_ONDJFM <- which(buoy_month == "1" | buoy_month == "2" | buoy_month == "3" | buoy_month == "10" | buoy_month == "11" | buoy_month == "12")
   #buoy_ONDJFM <- which(buoy_month == "4" | buoy_month == "5" | buoy_month == "6" | buoy_month == "7" | buoy_month == "8" | buoy_month == "9")
   qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc1_plot,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main=paste("Seasonal (",time_lab,") at NDBC 46066",sep=""))
   abline(a=0,b=1)

   #X11()
   #qqplot(hist_buoy_hs,mat_nc2_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main="Comparison at NDBC 46066")
# Seasonal.
   qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc2_plot,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main=paste("Seasonal (",time_lab,") at NDBC 46066",sep=""))
   abline(a=0,b=1)

   dev.off()

# Display.
   system(paste("okular",file_name,"&> /dev/null &"))


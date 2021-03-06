# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/qqplot_SRS_buoy.R")

# Script to load some SRS data from different satellites and compare with buoy NDBC 46066.
# BT 04/2019

# Libraries.
   library(ncdf4)
# QQ plot function.
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/func_qqplot.R")

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")

# File base name.
   #IMOS_SRS-Surface-Waves_MW_JASON-1_FV02_052N-205E-DM00.nc

# Buoy 46006, 40.782 N, 137.397 W
# Buoy 46066, 52.785 N, 155.047 W
# Buoy 51004, 17.604 N, 152.364 W
   buoy_list <- c("46006")

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
      mat_buoy <- NULL

      for (i in 1:no_files) {
         file_path[i] <- paste(buoy_data_dir,"/",data_files[i],sep='')
         temp_table1 <- read.table(file_path[i])
         temp_table <- read.table(file_path[i],skip=1,
                                  col.names=c("Y","M","D","H","m",rep("X",(dim(temp_table1)[2]-5))),
                                  colClasses=c("Y"="character","M"="character","D"="character","H"="character","m"="character"))
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
   hist_data2 <- seq(which(mat_buoy_obs[,1] == "2005")[1],dim(mat_buoy_obs)[1])
   hist_data1 <- seq(1,hist_data2[1]-1)

# Create a matrix of dates.
   date_temp1 <- mat_buoy_obs[hist_data1,]
   date_temp2 <- mat_buoy_obs[hist_data2,]
   mat_b_date <- rbind(cbind(date_temp1[!is.na(date_temp1[,9]),1:4],"00"),date_temp2[!is.na(date_temp2[,9]),1:5])
   colnames(mat_b_date) <- c("Y","M","D","H","m")
# Correct for missing centuries.
   AA <- mat_b_date[as.numeric(mat_b_date[,1]) <= 100,1]
   mat_b_date[as.numeric(mat_b_date[,1]) <= 100,1] <- as.character(as.numeric(AA)+1900)
# Function to write dates.
   func_buoy_date <- function(x) {
      #b_date <- paste(paste(mat_b_date[x,1:3],collapse='-'),paste(mat_b_date[x,4],mat_b_date[x,5],"00",sep=':'))
      b_date <- paste(paste(x[1:3],collapse='-'),paste(x[4],x[5],"00",sep=':'))
      return(b_date)
   }
   vec_buoy_dates <- apply(X=mat_b_date,MAR=1,FUN=func_buoy_date)
# Matrix of lists for date matching and searching.
   year_list <- c(mat_b_date[1,1]:mat_b_date[dim(mat_b_date)[1],1])
   mat_b_date_search <- matrix(list(),nrow=length(year_list),ncol=12)
   for (yy_idx in 1:dim(mat_b_date_search)[1]) {
      for (mm_idx in 1:12) {
         mat_b_date_search[[yy_idx,mm_idx]] <- which(as.numeric(mat_b_date[,1]) == year_list[yy_idx] & as.numeric(mat_b_date[,2]) == mm_idx)
      }
   }

# Hs data for histogram.
   hist_buoy_hs <- c(mat_buoy_obs[hist_data1,8],mat_buoy_obs[hist_data2,9])

# Quantiles for 
   q_plot <- quantile(hist_buoy_hs,probs=c(0.5,0.9),na.rm=T)

#-----------------------------------------------------------------------#
# Satellite data.
#-----------------------------------------------------------------------#
   fig_file_name <- paste("./figures/",buoy_name,"_5678_annual.png",sep="")
# Start figures.
   #X11()
   png(filename = fig_file_name,width=2200,height=2200)
   par(mfrow=c(2,2),oma=c(1.5,1.5,2,1),mar=c(7.0,7.0,5.0,3),mgp=c(5,2,0))

# Loop over different data sets.
   SHW_KU <- NULL

   #for (k in 1:length(vec_datasets) {
   for (k in c(5,6,7,8)) {
   #for (k in 6) {
      #vec_lat_sat1 <- paste("0",50:55,"N",sep="")
      #vec_lon_sat1 <- paste(200:210,"E",sep="")
# 46006
      nc1_filename <- paste(vec_datasets[k],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[k],"_FV02_040N-223E-DM00.nc",sep="")
# 46066
      #nc1_filename <- paste(vec_datasets[k],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[k],"_FV02_052N-205E-DM00.nc",sep="")
# 51004
      #nc1_filename <- paste(vec_datasets[k],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[k],"_FV02_017N-202E-DM00.nc",sep="")

      nc1 = nc_open(paste(data_path,nc1_filename,sep=""))
   
# SWH.
      nc1_SHW_KU <- ncvar_get(nc1,"SWH_KU")
      #SHW_KU <- c(SHW_KU,nc1_SHW_KU)
      nc1_SHW_C <- ncvar_get(nc1,"SWH_C")

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

#-----------------------------------------------------------------------#
# Compute temporal residual distribution.
#-----------------------------------------------------------------------#
# Loop over satellite obs.
      vec_min_time <- numeric(length(nc1_date))
      for (nc1_idx in 1:length(nc1_date)) {
         vec_search <- vec_buoy_dates[ mat_b_date_search[[ which(year_list == format(nc1_date[nc1_idx],"%Y")) , as.numeric(format(nc1_date[nc1_idx],"%m")) ]] ]
         vec_min_time[nc1_idx] <- tryCatch( min( abs( sapply(X=vec_search,FUN=difftime,time2=nc1_date[nc1_idx],units = "mins") ) ), error = function(e) { AA <- NA; return(AA) } )
      }
# Histogram.
      X11()
      hist(vec_min_time[vec_min_time < 201], xlim=c(0,200),breaks=25)

#-----------------------------------------------------------------------#
# Analysis:
#-----------------------------------------------------------------------#
# Wave height in each block. Matrix dim should not be the size of the entire record!
      mat_nc1_block_first <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      mat_nc1_block_max <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      mat_nc1_block_mean <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
      nc1_block_cor <- 0

      for (i in 1:dim(mat_nc1_breaks)[1]) {
         CC <- mat_nc1_breaks[i,1]:mat_nc1_breaks[i,2]
         mat_nc1_block_first[i,] <- c( nc1_SHW_KU[CC[1]],nc1_SHW_C[CC[1]] )
         mat_nc1_block_max[i,] <- c( max(nc1_SHW_KU[CC]),max(nc1_SHW_C[CC]) )
         mat_nc1_block_mean[i,] <- c( mean(nc1_SHW_KU[CC]),mean(nc1_SHW_C[CC]) )
         nc1_block_cor[i] <- cor(nc1_SHW_KU[CC],nc1_SHW_C[CC])
      }
# Seasonal info.
      nc1_month <- sapply(X=nc1_date[vec_nc1_break_mid],FUN=substr,start=6,stop=7)
      nc1_ONDJFM <- which(nc1_month == "01" | nc1_month == "02" | nc1_month == "03" | nc1_month == "10" | nc1_month == "11" | nc1_month == "12")

      #cor(mat_nc1_block_mean)

# Close files.
      nc_close(nc1)
   #}

   #X11()
   #qqplot(hist_buoy_hs,mat_nc1_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main="Comparison at NDBC 46066")
   func_qqplot(hist_buoy_hs[which(!is.na(hist_buoy_hs))],nc1_SHW_KU,
   #func_qqplot(hist_buoy_hs[which(!is.na(hist_buoy_hs))],mat_nc1_block_mean[,1],
          #xlim=c(0,12),ylim=c(0,12),
          xlim=c(0,4),ylim=c(0,4),regress = FALSE,
          xlab="Buoy",ylab="Satellite",
          main=paste(vec_datasets[k],": QQ at NDBC ",buoy_name,sep=""),cex.main=3.0,cex.lab=3.0,cex.axis=3.0,lwd.reg=5.0)
# Quantiles.
   for (kk in 1:2) {
      segments(x0=q_plot[kk], y0=0, x1=q_plot[kk], y1=q_plot[kk], lty=3, lwd=3, col="blue")
      segments(x0=0, y0=q_plot[kk], x1=q_plot[kk], y1=q_plot[kk], lty=3, lwd=3, col="blue")
   }
   mtext(text = paste("Q50 =",q_plot[1],"m"), side = 3, line = -5, adj = 0.04, cex = 2.5, col = "blue")
   mtext(text = paste("Q90 =",q_plot[2],"m"), side = 3, line = -8, adj = 0.04, cex = 2.5, col = "blue")
# Seasonal.
   #buoy_month <- mat_buoy_obs[,2]
   #buoy_ONDJFM <- which(buoy_month == "1" | buoy_month == "2" | buoy_month == "3" | buoy_month == "10" | buoy_month == "11" | buoy_month == "12")
   #qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc1_block_mean[nc1_ONDJFM,],xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 1",main="Seasonal (ONDJFM) at NDBC 46066")

   }

   #X11()
   #qqplot(hist_buoy_hs,mat_nc2_block_mean,xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main="Comparison at NDBC 46066")
# Seasonal.
   #qqplot(hist_buoy_hs[buoy_ONDJFM],mat_nc2_block_mean[nc2_ONDJFM,],xlim=c(0,12),ylim=c(0,12),xlab="Hs at NDBC 46066",ylab="Hs Jason 2",main="Seasonal (ONDJFM) at NDBC 46066")
   #abline(a=0,b=1)

   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))


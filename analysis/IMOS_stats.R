# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/IMOS_stats.R")

# Script to load a single dataset from IMOS and obtain summary statistics.
# BT 04/2019

# Libraries.
   library(ncdf4)
   library(extRemes)

# ================================================================= #
# Specify domain.
   lon_range <- 140:240
   lon_mid <- lon_range+0.5
   lat_range <- 30:60
   lat_mid <- lat_range+0.5

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")

# Data set selection.
   data_idx <- 5

# File base name.
   nc_sat1 <- paste(vec_datasets[data_idx],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[data_idx],"_FV02_",sep="")
   vec_lat_sat <- paste("0",rev(lat_range),"N",sep="")
   vec_lon_sat <- paste(lon_range,"E",sep="")
   vec_lonlat_sat <- paste(rep(vec_lat_sat,times=length(vec_lon_sat)),"-",rep(vec_lon_sat,each=length(vec_lat_sat)),sep="")
   vec_filenames_sat1 <- paste(nc_sat1,vec_lonlat_sat,"-DM00.nc",sep="")
   mat_filenames_sat1 <- matrix(vec_filenames_sat1,nrow=length(vec_lat_sat))

# ================================================================= #
# Set up variables.
   array_stats <- array(NA,dim=c(length(vec_lat_sat),ncol=length(vec_lon_sat),10,2))

# Loop over different locations.
   for (lon_idx in 1:length(lon_range)) {
   #for (lon_idx in 16) {
      for (lat_idx in 1:length(lat_range)) {
      #for (lat_idx in 4) {

#   X11()
#   #pdf(file = "./figures/J1_J2_46066_annual.pdf",width=10,height=3)
#   if ( length(file_list) == 9) {
#      par(mfrow=c(3,3))
#   } else {
#      par(mfrow=c(1,3))
#   }

# Test for filename existance.
      nc1_f <- paste(data_path,mat_filenames_sat1[lat_idx,lon_idx],sep="")

         if(! (file.exists(nc1_f) )) {
            print("File does not exist!")
            print(paste("LAT:",lat_range[lat_idx],"; LON:",lon_range[lon_idx]))
            array_stats[lat_idx,lon_idx,,] <- NA
         } else {
# Open file.
            nc1 = nc_open(paste(data_path,mat_filenames_sat1[lat_idx,lon_idx],sep=""))
   
# SWH.
            nc1_SHW_KU <- ncvar_get(nc1,"SWH_KU")
            nc1_SHW_C <- ncvar_get(nc1,"SWH_C")
            mat_nc1 <- cbind(nc1_SHW_KU,nc1_SHW_C)

# Summary data.
# Raw counts [1].
            array_stats[lat_idx,lon_idx,1,1:2] <- apply(X=mat_nc1,MAR=2,length)
# Pass counts [2].
# See below.
# Mean [3].
            array_stats[lat_idx,lon_idx,3,1:2] <- apply(X=mat_nc1,MAR=2,mean)
            #array_sum_out[lat_idx,lon_idx,1,3] <- mean(c(nc1_SHW_KU,nc1_SHW_C))
# Variance [4].
            array_stats[lat_idx,lon_idx,4,1:2] <- apply(X=mat_nc1,MAR=2,var)
# Quantiles [5:9].
            array_stats[lat_idx,lon_idx,5:9,1:2] <- apply(X=mat_nc1,MAR=2,FUN=quantile,probs=c(0.5,0.75,0.90,0.95,0.99),na.rm=T)
# Max [10].
            array_stats[lat_idx,lon_idx,10,1:2] <- apply(X=mat_nc1,MAR=2,max)

# Time.
            nc1_time_idx <- ncvar_get(nc1,"TIME")
            nc1_date <- as.POSIXct(nc1_time_idx*3600*24, origin = '1985-01-01', tz='GMT')
            nc1_start <- nc1_date[1]

            nc1_lon <- ncvar_get(nc1,"LONGITUDE")
            nc1_lat <- ncvar_get(nc1,"LATITUDE")

# Look for breaks.
            nc1_breaks <- NULL
            if ( length(nc1_time_idx) > 1 ) {
               for (i in 2:length(nc1_time_idx)) {
                  AA <- nc1_time_idx[i] - nc1_time_idx[i-1]
                     if ( AA > 1e-4 ) {
                     nc1_breaks <- c(nc1_breaks,i)
                  #print(paste(" Break before:",i))
                  }
               }
            }
   
# Matrix of 'block' indices for contiguous readings.
            mat_nc1_breaks <- cbind(c(1,nc1_breaks),c(nc1_breaks-1,length(nc1_time_idx)))
            vec_nc1_break_mid <- apply(X=mat_nc1_breaks,MAR=1,FUN=function(x) { floor(mean(x[1]:x[2])) })

# Close file.
            nc_close(nc1)

# Pass counts [2].
            array_stats[lat_idx,lon_idx,2,1:2] <- dim(mat_nc1_breaks)[1]

# ================================================================= #
# Analysis:
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
## Seasonal info.
#   nc1_month <- sapply(X=nc1_date[vec_nc1_break_mid],FUN=substr,start=6,stop=7)
#   nc1_ONDJFM <- which(nc1_month == "01" | nc1_month == "02" | nc1_month == "03" | nc1_month == "10" | nc1_month == "11" | nc1_month == "12")

         }
      }
   }


# Create data structure including metadata.
   array_meta <- list(data_idx=data_idx,data_name=vec_datasets[data_idx],lat=lat_mid,lon=lon_mid,lab_stat=c("raw_counts","pass_counts","mean","variance","Q50","Q75","Q90","Q95","Q99","maxium"))   
   list_IMOS_stats <- list(array_meta,array_stats)
# Write out data.
   stats_file <- paste("./test_output/list_",vec_datasets[data_idx],"_stats.Robj",sep="")
   save(list_IMOS_stats,file = stats_file)


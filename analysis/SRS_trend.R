# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/SRS_trend.R")

# Script to load multiple datasets from IMOS and obtain summary statistics and annual trend.
# BT 04/2019

# Libraries.
   library(ncdf4)
   library(abind)
   library(extRemes)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")

# ================================================================= #
# Specify domain.
   lon_range <- 190:200
   #lon_range <- 210:220
   lon_mid <- lon_range+0.5
   lat_range <- 30:40
   #lat_range <- 40:50
   lat_mid <- lat_range+0.5

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   mat_valid_dates <- cbind(c(1986,1990),c(1992,1995),c(1993,2004),c(1996,2008),c(1999,2007),c(2002,2012),c(2009,2017),c(2011,2017),c(2012,2017),c(2014,2017),c(2016,2017),c(2017,2017))

# Data set selection.
   mission_idx <- 2:6

# File base name.
   #array_filenames <- array(NA,dim=c(length(lat_range),length(lon_range),length(mission_idx)))
   array_filenames <- NULL
   for (m_idx in mission_idx) {
      nc_sat1 <- paste(vec_datasets[m_idx],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[m_idx],"_FV02_",sep="")
      vec_lat_sat <- paste("0",rev(lat_range),"N",sep="")
      vec_lon_sat <- paste(lon_range,"E",sep="")
      vec_lonlat_sat <- paste(rep(vec_lat_sat,times=length(vec_lon_sat)),"-",rep(vec_lon_sat,each=length(vec_lat_sat)),sep="")
      vec_filenames_sat1 <- paste(nc_sat1,vec_lonlat_sat,"-DM00.nc",sep="")
      mat_filenames_sat1 <- matrix(vec_filenames_sat1,nrow=length(vec_lat_sat))
      array_filenames <- abind(array_filenames,mat_filenames_sat1,along=3)
   }

# Analysis duration (typically the coverage of the chosen missions).
   anal_years <- min(mat_valid_dates[,mission_idx]):max(mat_valid_dates[,mission_idx])

# ================================================================= #
# Set up variables.
# Matrix to store "block" Q50s and time index per mission.
   array_block <- array(list(),dim=c(length(vec_lat_sat),length(vec_lon_sat),2,length(mission_idx)))

# Loop over different locations.
   for (lon_idx in 1:length(lon_range)) {
   #for (lon_idx in 16:17) {
      for (lat_idx in 1:length(lat_range)) {
      #for (lat_idx in 10:11) {

# Loop over missions.
# Data structures.
         list_nc1_block_Q50 <- list(length(mission_idx))
# Matrix to store "block" Q50s and time index per mission.
         #mat_nc1_block <- matrix(list(),nrow=2,ncol=length(mission_idx))
# Loop.
         for (m_idx in 1:length(mission_idx)) {

# Test for filename existance.
            nc1_f <- paste(data_path,array_filenames[lat_idx,lon_idx,m_idx],sep="")

            if(! (file.exists(nc1_f) )) {
               print(paste(nc1_f,": File does not exist!"))
               print(paste("LAT:",lat_range[lat_idx],"; LON:",lon_range[lon_idx]))
               array_block[lat_idx,lon_idx,,m_idx] <- NA
            } else {
# Open file.
               nc1 = nc_open(nc1_f)
   
# Time.
               nc1_time_idx <- ncvar_get(nc1,"TIME")
               nc1_date <- as.POSIXct(nc1_time_idx*3600*24, origin = '1985-01-01', tz='GMT')
               nc1_start <- nc1_date[1]

               nc1_lon <- ncvar_get(nc1,"LONGITUDE")
               nc1_lat <- ncvar_get(nc1,"LATITUDE")

# SWH.
               nc1_SHW_KU <- ncvar_get(nc1,"SWH_KU")
               #nc1_SHW_C <- ncvar_get(nc1,"SWH_C")
               #mat_nc1 <- cbind(nc1_SHW_KU,nc1_SHW_C)

# Close file.
               nc_close(nc1)

# ================================================================= #
# Find passes through grid cell (blocks) and take an average (Q50).
# ================================================================= #
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

# Time index per block (mean value, although the spacing is only ~1s).
               #mat_nc1_block[[1,m_idx]] <- nc1_date[vec_nc1_break_mid]
               array_block[[lat_idx,lon_idx,1,m_idx]] <- nc1_date[vec_nc1_break_mid]

# Wave height in each block. Matrix dim should not be the size of the entire record!
               #mat_nc1_block_mean <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
               #mat_nc1_block_Q50 <- matrix(0,nrow=dim(mat_nc1_breaks)[1],ncol=2)
               temp_block_Q50 <- numeric(dim(mat_nc1_breaks)[1])

               for (i in 1:dim(mat_nc1_breaks)[1]) {
                  CC <- mat_nc1_breaks[i,1]:mat_nc1_breaks[i,2]
                  #mat_nc1_block_mean[i,] <- c( mean(nc1_SHW_KU[CC]),mean(nc1_SHW_C[CC]) )
                  #mat_nc1_block_Q50[i,] <- c( median(nc1_SHW_KU[CC]),median(nc1_SHW_C[CC]) )
                  temp_block_Q50[i] <- median(nc1_SHW_KU[CC])
               }
               #mat_nc1_block[[2,m_idx]] <- temp_block_Q50
               array_block[[lat_idx,lon_idx,2,m_idx]] <- temp_block_Q50
            }
         }
      }
   }

# Create data structure including metadata.
   array_meta <- list(mission_idx=mission_idx,data_name=vec_datasets[mission_idx],lat=lat_mid,lon=lon_mid,lab_stat=c("date","block_Q50"))
   list_SRS_blocks <- list(array_meta,array_block)
# Write out data.
   data_file <- paste("./output/test_block/list_",paste(vec_datasets[mission_idx],collapse='_'),"_blocks.Robj",sep="")
   save(list_SRS_blocks,file = data_file)

# ================================================================= #
# Summary data and temporal processing.
# Loop over years to get obs per year across all missions.

# Data structures.
   mat_list_annual_KU <- matrix(list(),nrow=length(vec_lat_sat),ncol=length(vec_lon_sat))
   mat_list_annual_trend <- matrix(list(),nrow=length(vec_lat_sat),ncol=length(vec_lon_sat))
   vec_q <- c(0.5,0.9,0.95)

# Loop.
   for (lon_idx in 1:length(lon_range)) {
      for (lat_idx in 1:length(lat_range)) {
# Array for annual data (stats, CIs).
         array_annual_stats <- array(NA,dim=c(length(anal_years),length(vec_q),3))
         if (!all(is.na(unlist(array_block[lat_idx,lon_idx,,])))) {
            for (y_idx in 1:length(anal_years)) {
               temp_annual_hs <- NULL
               for (m_idx in 1:length(mission_idx)) {
                  if (!is.na(array_block[lat_idx,lon_idx,1,m_idx])) {
                     temp_annual_hs <- c(temp_annual_hs,array_block[[lat_idx,lon_idx,2,m_idx]][which(format(array_block[[lat_idx,lon_idx,1,m_idx]],"%Y") == anal_years[y_idx])])
                  }
               }
# Do stats.
               array_annual_stats[y_idx,,2] <- quantile(temp_annual_hs,probs=vec_q)
# Do CIs.
               array_annual_stats[y_idx,,c(1,3)] <- t(sapply(X=vec_q,FUN=function(x) { sort(temp_annual_hs)[quantile.CI(length(temp_annual_hs),q=x)$Interval] }))
            }
            mat_list_annual_KU[[lat_idx,lon_idx]] <- array_annual_stats
# Trend (linear regression).
            mat_b_q <- mat_list_annual_KU[[lat_idx,lon_idx]][,,2]
            colnames(mat_b_q) <- paste("Q",100*vec_q,sep="")
            df_Q <- data.frame(cbind(year=anal_years,mat_b_q))

            mat_trend <- matrix(NA,nrow=length(vec_q),ncol=2)
            for (qq in 1:length(vec_q)) {
# Catch if lack of data for regression.
               if (!all(is.na(df_Q[,(qq+1)]))) {
                  lm_Q <- lm(df_Q[,(qq+1)] ~ year,data=df_Q)
                  sum_lm_Q <- summary(lm_Q)
                  mat_trend[qq,] <- c(sum_lm_Q$coefficients[2],sum_lm_Q$coefficients[8])
               }
            }
            colnames(mat_trend) <- c("year_slope","Pr(>|t|)")
            mat_list_annual_trend[[lat_idx,lon_idx]] <- mat_trend
         } else {
            mat_list_annual_KU[[lat_idx,lon_idx]] <- NA
            mat_list_annual_trend[[lat_idx,lon_idx]] <- NA
         }
      }
   }

# Create data structure for saving, including metadata.
   array_meta <- c(list_SRS_blocks[[1]],list(trend_stats=paste("Q",100*c(0.5,0.9,0.95),sep=""),trend=c("slope","P-val")))
   list_SRS_trend <- list(array_meta,mat_list_annual_trend)
# Write out data.
   data_file <- paste("./output/test_block/list_",paste(vec_datasets[mission_idx],collapse='_'),"_trend.Robj",sep="")
   save(list_SRS_trend,file = data_file)

## ================================================================= #
## Plotting the trend for a single grid cell.
#   mat_b_q <- mat_list_annual_KU[[lat_idx,lon_idx]][,,2]
#   colnames(mat_b_q) <- paste("Q",100*c(0.5,0.9,0.95),sep="")
#   df_Q <- data.frame(cbind(year=anal_years,mat_b_q))
#
## Weights for weighted least squares (crude).
#   array_q_CI <- abind(mat_list_annual_KU[[lat_idx,lon_idx]][,,1],mat_list_annual_KU[[lat_idx,lon_idx]][,,3],along=3)
#   mat_CI_w <- matrix(NA,nrow=length(anal_years),ncol=length(vec_q))
#   for (qq in 1:length(vec_q)) {
#      mat_CI_w[,qq] <- max(mat_list_annual_KU[[lat_idx,lon_idx]][,qq,3] - mat_list_annual_KU[[lat_idx,lon_idx]][,qq,1],na.rm=T) /
#                  (mat_list_annual_KU[[lat_idx,lon_idx]][,qq,3] - mat_list_annual_KU[[lat_idx,lon_idx]][,qq,1])
#   }
#   #vec_CI_w[is.na(vec_CI_w)] <- 1
#
## Plot.
#   fig_file_name <- paste("./figures/",paste(c(lat_range[lat_idx],lon_range[lon_idx],vec_datasets[mission_idx]),collapse='_'),"_annual_trend.png",sep="")
#   #X11()
#   png(filename = fig_file_name, width = 2200, height = 2200)
#   par(mfrow=c(2,2),oma=c(1.5,1.5,3,1),mar=c(7.0,7.0,5.0,3),mgp=c(5,2,0))
#   for (qq in 1:length(vec_q)) {
## Weights.
#      #lm_Q <- lm(df_Q[,(qq+1)] ~ year,data=df_Q,weights=vec_CI_w)
## No weighting.
#      lm_Q <- lm(df_Q[,(qq+1)] ~ year,data=df_Q)
#      sum_lm_Q <- summary(lm_Q)
#      print(summary(lm_Q))
#
#      plot(df_Q[,c(1,(qq+1))],ylim=c(0,7),cex.main=3.0,cex.lab=3.0,cex.axis=3.0,lwd=3.0)
#      arrows(df_Q$year, array_q_CI[,qq,1], df_Q$year, array_q_CI[,qq,2], length=0.05, angle=90, code=3, lwd=2)
#      abline(lm_Q,lwd=3)
#
#      if ( sum_lm_Q$coefficients[8] < 0.1 ) {
#         trend_sig <- paste("Trend: ",format(sum_lm_Q$coefficients[2],digits=2),"m per year. Significant at 10%, Pr(>|t|) = ",format(sum_lm_Q$coefficients[8],digits=2),sep="")
#      } else {
#         trend_sig <- paste("Trend: ",format(sum_lm_Q$coefficients[2],digits=2),"m per year. Not significant, Pr(>|t|) = ",format(sum_lm_Q$coefficients[8],digits=2),sep="")
#      }
#      mtext(text = trend_sig, side = 3, line = -5, cex = 3)
#   }
#   mtext(text = paste("Trend in Q50, Q90, Q95 at LON:",lon_mid[lon_idx],"LAT:",lat_mid[lat_idx],"from",paste(vec_datasets[mission_idx],collapse=', ')), outer = TRUE, side = 3, line = -2, cex = 4)
#   dev.off()
#   system(paste("okular",fig_file_name,"&> /dev/null &"))
#
### Seasonal info.
##   nc1_month <- sapply(X=nc1_date[vec_nc1_break_mid],FUN=substr,start=6,stop=7)
##   nc1_ONDJFM <- which(nc1_month == "01" | nc1_month == "02" | nc1_month == "03" | nc1_month == "10" | nc1_month == "11" | nc1_month == "12")
##

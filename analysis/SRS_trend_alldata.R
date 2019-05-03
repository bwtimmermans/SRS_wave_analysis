# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/SRS_trend_alldata.R")

# Script to load multiple datasets from IMOS and obtain summary statistics and annual trend.
# BT 04/2019

# Libraries.
   library(ncdf4)
   library(abind)
   library(extRemes)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")

# ================================================================= #
# Edit here
# ----------------------------------------------------------------- #

# Specify domain.
   lon_range <- 140:239
   #lon_range <- 210:219
   lon_mid <- lon_range+0.5
   lat_range <- 0:59
   #lat_range <- 0:11
   lat_mid <- lat_range+0.5

# Resolution (1 = 1 degree, 2 = 2 degree, etc).
   res <- 2
   mat_lat_grid_idx <- matrix(1:length(lat_range),nrow=res)
   mat_lon_grid_idx <- matrix(1:length(lon_range),nrow=res)

# Data set selection.
   mission_idx <- 2:6

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   mat_valid_dates <- cbind(c(1986,1990),c(1992,1995),c(1993,2004),c(1996,2008),c(1999,2007),c(2002,2012),c(2009,2017),c(2011,2017),c(2012,2017),c(2014,2017),c(2016,2017),c(2017,2017))

# File base name.
   #array_filenames <- array(NA,dim=c(length(lat_range),length(lon_range),length(mission_idx)))
   array_filenames <- NULL
   for (m_idx in mission_idx) {
      nc_sat1 <- paste(vec_datasets[m_idx],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[m_idx],"_FV02_",sep="")
      vec_lat_sat <- unlist( lapply(X=rev(lat_range),FUN=function(x) { paste(tail(strsplit(paste("00",x,"N",sep=""),split='')[[1]],n=4),collapse='') } ) )
      vec_lon_sat <- paste(lon_range,"E",sep="")
      vec_lonlat_sat <- paste(rep(vec_lat_sat,times=length(vec_lon_sat)),"-",rep(vec_lon_sat,each=length(vec_lat_sat)),sep="")
      vec_filenames_sat1 <- paste(nc_sat1,vec_lonlat_sat,"-DM00.nc",sep="")
      mat_filenames_sat1 <- matrix(vec_filenames_sat1,nrow=length(vec_lat_sat))
      array_filenames <- abind(array_filenames,mat_filenames_sat1,along=3)
   }

# Analysis duration (typically the coverage of the chosen missions).
   anal_years <- min(mat_valid_dates[,mission_idx]):max(mat_valid_dates[,mission_idx])

# ================================================================= #
# Data structures.
   mat_list_annual_KU <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
   mat_list_annual_trend <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
   vec_q <- c(0.5,0.9,0.95)

# Loop over cells by resolution.
   for (lon_res_idx in 1:dim(mat_lon_grid_idx)[2]) {
      for (lat_res_idx in 1:dim(mat_lat_grid_idx)[2]) {
# Loop over missions.
# Data structures.
# Matrix to store time index and SWH obs per mission.
         #mat_block_temp <- matrix(list(),nrow=2,ncol=length(mission_idx))
         mat_annual_hs <- matrix(list(),nrow=length(anal_years),ncol=length(mission_idx))
# Loop.
         for (m_idx in 1:length(mission_idx)) {
# Data structures.
            nc1_time_idx <- NULL
            nc1_SHW_KU <- NULL

# Loop over res "sub" indices.
            for (sub_lon_idx in 1:dim(mat_lon_grid_idx)[1]) {
               lon_idx <- mat_lon_grid_idx[sub_lon_idx,lon_res_idx]
               for (sub_lat_idx in 1:dim(mat_lat_grid_idx)[1]) {
                  lat_idx <- mat_lat_grid_idx[sub_lat_idx,lat_res_idx]

# File loading.
# Test for filename existance.
                  nc1_f <- paste(data_path,array_filenames[lat_idx,lon_idx,m_idx],sep="")

                  if(! (file.exists(nc1_f) )) {
                     print(paste(nc1_f,": File does not exist!"))
                     print(paste("LAT:",lat_range[lat_idx],"; LON:",lon_range[lon_idx]))
                     #mat_block_temp[,m_idx] <- NA
                  } else {
# Open file.
                     nc1 = nc_open(nc1_f)
                     #nc1_lon <- ncvar_get(nc1,"LONGITUDE")
                     #nc1_lat <- ncvar_get(nc1,"LATITUDE")
# Time.
                     nc1_time_idx <- c(nc1_time_idx,ncvar_get(nc1,"TIME"))
# SWH.
                     nc1_SHW_KU <- c(nc1_SHW_KU,ncvar_get(nc1,"SWH_KU"))
                     #nc1_SHW_C <- ncvar_get(nc1,"SWH_C")

# Close file.
                     nc_close(nc1)

                  }
               }
            }
# Time index and wave height obs.
            #mat_block_temp[[1,m_idx]] <- nc1_date_temp
            #mat_block_temp[[2,m_idx]] <- nc1_SHW_KU_temp
# Temporal processing.
# Loop over years to get obs per year across all missions.
            if (!is.null(nc1_time_idx)) {
               nc1_date_temp <- as.POSIXct(nc1_time_idx*3600*24, origin = '1985-01-01', tz='GMT')
               for (y_idx in 1:length(anal_years)) {
                  mat_annual_hs[[y_idx,m_idx]] <- nc1_SHW_KU[which(format(nc1_date_temp,"%Y") == anal_years[y_idx])]
               }
            }
         }
         if (!is.null(unlist(mat_annual_hs))) {
# Array for annual data (stats, CIs) and temps to hold sub-grid data.
            array_annual_stats <- array(NA,dim=c(length(anal_years),length(vec_q),3))
# Loop over years to get total obs from each mission.
            for (y_idx in 1:length(anal_years)) {
               temp_annual_hs <- unlist(mat_annual_hs[y_idx,])
# Do stats.
               array_annual_stats[y_idx,,2] <- quantile(temp_annual_hs,probs=vec_q)
# Do CIs.
               array_annual_stats[y_idx,,c(1,3)] <- t(sapply(X=vec_q,FUN=function(x) { sort(temp_annual_hs)[quantile.CI(length(temp_annual_hs),q=x)$Interval] }))
            }
# Trend (linear regression).
            mat_b_q <- array_annual_stats[,,2]
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
# Store for writing.
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- array_annual_stats
            mat_list_annual_trend[[lat_res_idx,lon_res_idx]] <- mat_trend
         } else{
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- NA
            mat_list_annual_trend[[lat_res_idx,lon_res_idx]] <- NA
         }
      }
   }

# Create data structure including metadata.
   array_meta <- list(mission_idx=mission_idx,data_name=vec_datasets[mission_idx],
                      orig_lat_cell=lat_range, orig_lon_cell=lon_range,
                      lat_cell=(matrix(lat_range,nrow=res)[1,] + (res/2)),lon_cell=(matrix(lon_range,nrow=res)[1,] + (res/2)),
                      lat_mid=(matrix(lat_mid,nrow=res)[1,] + (res/2)),lon_mid=(matrix(lon_mid,nrow=res)[1,] + (res/2)),
                      trend_stats=paste("Q",100*c(0.5,0.9,0.95),sep=""),trend=c("slope","P-val"))
   list_SRS_trend <- list(array_meta,mat_list_annual_trend)
# Write out data.
   data_file <- paste("./output/test_block/list_",paste(vec_datasets[mission_idx],collapse='_'),"_all_data_trend.Robj",sep="")
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

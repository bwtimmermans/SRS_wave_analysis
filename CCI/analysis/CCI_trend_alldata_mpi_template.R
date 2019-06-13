# source("/home/ben/research/NOC/SRS_wave_analysis/CCI/analysis/CCI_trend_alldata_mpi.R")

# Script to load CCI dataset and obtain long term trends.
# BT 04/2019

# Libraries.
   #.comm.size <- 4; .comm.rank <- 3
   library(ncdf4)
   library(abind)
# MPI.
   library(pbdMPI, quietly = TRUE)
   init()
   .comm.size <- comm.size()
   .comm.rank <- comm.rank()

# ================================================================= #
# Edit here
# ----------------------------------------------------------------- #

# Data tpye (KU-all, KU-pass, CU etc).
   data_type <- "KU"

# Specify domain.
   lon_range <- 0:359
   #lon_range <- 210:219
   lon_mid <- lon_range + 0.5
   lat_range <- -72:71
   #lat_range <- 0:11
   lat_mid <- lat_range + 0.5

# Actual domain.
   nc1 = nc_open("/backup/datasets/CCI/ESACCI-SEASTATE-L4-SWH-MULTI_1M-201001-fv01.nc")
   vec_lon <- ncvar_get(nc1,"lon")
   vec_lat <- ncvar_get(nc1,"lat")
   nc_close(nc1)

# Resolution (1 = 1 degree, 2 = 2 degree, etc).
   res <- 4

# Analysis duration.
# Note for 'flag_winter', first year must be >= 1993.
   anal_years <- <anal_years>

# Flag for complete winter season.
   flag_winter <- FALSE
   if (flag_winter) { lab_y_centre = "n_winter" } else { lab_y_centre = "n_summer" }

   lab_years <- paste(anal_years[c(1,length(anal_years))],collapse='-')
   all_months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

# Months for analysis.
   flag_annual <- FALSE
   #anal_months <- c("01","02","03")
   #anal_months <- c("04","05","06")
   #anal_months <- c("07","08","09")
   #anal_months <- c("10","11","12")
   #anal_months <- c("01","02","11","12")
   anal_months <- c("01","02","03","10","11","12")
   #anal_months <- c("04","05","06","07","08","09")
   #anal_months <- c("01","02","03","04","05","06")
   #anal_months <- c("07","08","09","10","11","12")
   if (!flag_annual) {
      lab_months <- paste(c("J","F","M","A","M","J","J","A","S","O","N","D")[as.numeric(anal_months)],collapse='')
   } else {
      lab_months <- "annual"
   }

# Flag for regression.
   #flag_reg <- "ONI"
   #flag_reg <- "PDO"
   #flag_reg <- "NAO"
   #flag_reg <- "AO"
   #flag_reg <- "none"
   flag_reg <- "<flag_reg>"

# Parallelise over the geographic range, by longitude.
# Parallelise over longitude, divide the range by number of processing cores.
   dim.lon <- length(lon_range)
   lon.base <- floor(dim.lon / .comm.size)
   #vec.lon_dim <- rep( lon.base, .comm.size ) + c( rep(1,dim.lon %% .comm.size), rep(0,(.comm.size - dim.lon %% .comm.size)) )
   vec.lon_dim <- res * ( floor( ( dim.lon / res ) / .comm.size) + c( rep(1,( dim.lon / res ) %% .comm.size), rep(0,(.comm.size - ( dim.lon / res ) %% .comm.size)) ) )

   lon.start.idx <- sum( vec.lon_dim[c(0:.comm.rank)] ) + 1
   lon.dim.node <- vec.lon_dim[.comm.rank+1]

   lon_range_node <- lon_range[lon.start.idx:(lon.start.idx + lon.dim.node - 1)]
   lon.start.idx.nc <- which(vec_lon == (lon_range_node - 179.5)[1])
   print(paste("lon_range_node:",lon_range_node))

# Latitude.
   lat.start.idx <- which(vec_lat == lat_mid[1])
   lat.stop.idx <- which(vec_lat == lat_mid[length(lat_mid)])
   dim.lat <- length(lat_range)

# Grid indices for aggregation.
   mat_lat_grid_idx <- matrix(1:length(lat_range),nrow=res)
   mat_lon_grid_idx <- matrix(1:length(lon_range_node),nrow=res)
   #mat_lon_grid_idx <- matrix(lon.start.idx:(lon.start.idx+length(lon_range_node)-1),nrow=res)

# ================================================================= #
# Specify data files and load all data (only summary statistics, so do all in one).

# Data path.
   data_path <- "/backup/datasets/CCI/"

# File base name.
   mat_filenames <- matrix(NA,nrow=length(anal_years),ncol=12)
   for (y_idx in 1:length(anal_years)) {
      mat_filenames[y_idx,] <- paste(data_path,"ESACCI-SEASTATE-L4-SWH-MULTI_1M-",anal_years[y_idx],all_months,"-fv01.nc",sep='')
   }

# Load data.
   array_CCI <- array(NA,dim=c(length(lon_range_node),length(lat_range),length(anal_years),12,3))
   mat_time <- matrix(NA,nrow=length(anal_years),ncol=12)
   for (y_idx in 1:length(anal_years)) {
      for (m_idx in 1:12) {
         nc1 = nc_open(mat_filenames[y_idx,m_idx])
         mat_time[y_idx,m_idx] <- ncvar_get(nc1,"time")
# Statistics.
         array_CCI[,,y_idx,m_idx,1] <- ncvar_get(nc1,"swh_count",start=c(lon.start.idx.nc,lat.start.idx,1), count=c((lon.dim.node),dim.lat,1))
         array_CCI[,,y_idx,m_idx,2] <- ncvar_get(nc1,"swh_mean",start=c(lon.start.idx.nc,lat.start.idx,1), count=c((lon.dim.node),dim.lat,1))
# Squared sum appears to be average (i.e. sum of squares / n).
         array_CCI[,,y_idx,m_idx,3] <- ncvar_get(nc1,"swh_squared_sum",start=c(lon.start.idx.nc,lat.start.idx,1), count=c(lon.dim.node,dim.lat,1))
         #array_CCI[,,y_idx,m_idx,4] <- ncvar_get(nc1,"swh_sum")
         #array_CCI[,,y_idx,m_idx,5] <- ncvar_get(nc1,"swh_rms")

         nc_close(nc1)
      }
   }
# Clean up missing values.
   array_CCI[array_CCI > 50] <- NA

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
   ONI_years <- which( as.numeric(rownames(mat_ONI)) == anal_years[1] ):which( as.numeric(rownames(mat_ONI)) == anal_years[length(anal_years)] )

# NAO.
   df_NAO <- read.csv("/home/ben/research/NOC/SRS_wave_analysis/datasets/indices/norm.nao.monthly.b5001.current.ascii.table",header=FALSE,sep=',')
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
   NAO_years <- which( as.numeric(rownames(mat_NAO)) == anal_years[1] ):which( as.numeric(rownames(mat_NAO)) == anal_years[length(anal_years)] )

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
   PDO_years <- which( as.numeric(rownames(mat_PDO)) == anal_years[1] ):which( as.numeric(rownames(mat_PDO)) == anal_years[length(anal_years)] )

# AO.
   df_AO <- read.csv("/home/ben/research/NOC/SRS_wave_analysis/datasets/indices/AO_NOAA.csv",header=TRUE, sep=',')
# Matrix containing mean, var and "min or max".
   mat_AO_temp <- t(matrix(df_AO$Value[-c(829:833)],nrow=12))

   mat_AO <- matrix(NA,nrow=dim(mat_AO_temp)[1],ncol=3)
   mat_AO[,1] <- apply(X=mat_AO_temp,MAR=1,FUN=mean)
   mat_AO[,2] <- apply(X=mat_AO_temp,MAR=1,FUN=var)
# Min or max based upon whether the mean is positive or negative.
   for (i in 1:dim(mat_AO_temp)[1]) {
      if (mat_AO[i,1] < 0) {
         mat_AO[i,3] <- min(mat_AO_temp[i,])
      } else {
         mat_AO[i,3] <- max(mat_AO_temp[i,])
      }
   }
# Column and row names.
   colnames(mat_AO) <- c("AO_mean","AO_var","AO_minmax")
   rownames(mat_AO) <- 1950:2018
# AO range based upon 2018.
   AO_years <- which( as.numeric(rownames(mat_AO)) == anal_years[1] ):which( as.numeric(rownames(mat_AO)) == anal_years[length(anal_years)] )

# ================================================================= #
# Data structures.
   mat_list_annual_KU <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
   mat_list_annual_trend_node <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
   vec_q <- c(0.5,0.9,0.95)

# Months for continuous winter.
   vec_anal_months <- as.numeric(anal_months)
   list_split_months <- list(2)
   list_split_months[[1]] <- vec_anal_months[vec_anal_months > 6]
   list_split_months[[2]] <- vec_anal_months[vec_anal_months <= 6]

# Loop over cells by resolution.
   for (lon_res_idx in 1:dim(mat_lon_grid_idx)[2]) {
      for (lat_res_idx in 1:dim(mat_lat_grid_idx)[2]) {
# Loop over years.
# Data structures.
         mat_stats_temp <- matrix(NA,nrow=length(anal_years),4)
# Data aggregation.
# Test for at least 75% observations.
         if ( (sum(!is.na(as.vector(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],1,,2]))) /
              length(as.vector(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],1,,2]))) > 0.50 ) {
# Centre year on winter season.
         if (flag_winter) {
            if (flag_annual) {
# Loop over statistics for aggregation.
               for (y_idx in 2:length(anal_years)) {
                  mat_stats_temp[y_idx,1] <- sum(
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),7:12,1],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,1:6,1],na.rm=T
                                                ) 
                  mat_stats_temp[y_idx,2] <- mean(
                                                  c(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),7:12,2],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,1:6,2]),na.rm=T
                                                 ) 
                  mat_stats_temp[y_idx,3] <- mean(
                                                  c(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),7:12,3],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,1:6,3]),na.rm=T
                                                 ) 
               }
               mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
            } else {
               for (y_idx in 2:length(anal_years)) {
                  mat_stats_temp[y_idx,1] <- sum(
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),list_split_months[[1]],1],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,list_split_months[[2]],1],na.rm=T
                                                ) 
                  mat_stats_temp[y_idx,2] <- mean(
                                                  c(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),list_split_months[[1]],2],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,list_split_months[[2]],2]),na.rm=T
                                                 ) 
                  mat_stats_temp[y_idx,3] <- mean(
                                                  c(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),list_split_months[[1]],3],
                                                  array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,list_split_months[[2]],3]),na.rm=T
                                                 ) 
               }
               mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
            }
         } else {
# Centre year on summer (normal).
            if (flag_annual) {
# Loop over statistics for aggregation.
               mat_stats_temp[,1] <- sapply(X=1:length(anal_years),FUN=function(x) { sum(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,,1],na.rm=T) } )
               mat_stats_temp[,2] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,,2],na.rm=T) } )
               mat_stats_temp[,3] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,,3],na.rm=T) } )
               mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
            } else {
               mat_stats_temp[,1] <- sapply(X=1:length(anal_years),FUN=function(x) { sum(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,as.numeric(anal_months),1],na.rm=T) } )
               mat_stats_temp[,2] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,as.numeric(anal_months),2],na.rm=T) } )
               mat_stats_temp[,3] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_CCI[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,as.numeric(anal_months),3],na.rm=T) } )
               mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
            }
         }
         df_stats <- as.data.frame(mat_stats_temp)
         colnames(df_stats) <- c("swh_counts","swh_mean","swh_squared_sum","swh_var")
         #rownames(df_stats) <- anal_years

# Trend (linear regression).
         df_Q <- data.frame(cbind(year=anal_years,df_stats,mat_ONI[ONI_years,],mat_NAO[NAO_years,],mat_PDO[PDO_years,],mat_AO[AO_years,]))

         mat_trend <- matrix(NA,nrow=length(vec_q),ncol=2)
         for (qq in 1:1) {
# Catch if lack of data for regression.
            if (all(!is.nan(df_Q[,(qq+2)]))) {
               if (flag_reg == "NAO") {
                  lm_Q <- lm(df_Q[,(qq+2)] ~ year + NAO_mean,data=df_Q)
               } else if (flag_reg == "ONI") {
                  lm_Q <- lm(df_Q[,(qq+2)] ~ year + ONI_mean,data=df_Q)
               } else if (flag_reg == "PDO") {
                  lm_Q <- lm(df_Q[,(qq+2)] ~ year + PDO_mean,data=df_Q)
               } else if (flag_reg == "AO") {
                  lm_Q <- lm(df_Q[,(qq+2)] ~ year + AO_mean,data=df_Q)
               } else {
                  lm_Q <- lm(df_Q[,(qq+2)] ~ year,data=df_Q)
                  #lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q)
               }
               sum_lm_Q <- summary(lm_Q)
               mat_trend[qq,] <- c(sum_lm_Q$coefficients[2,1],sum_lm_Q$coefficients[2,4])
            }
         }
         colnames(mat_trend) <- c("year_slope","Pr(>|t|)")
         if (!all(is.na(mat_trend))) {
# Store for writing.
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- df_stats
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- mat_trend
         } else{
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- NA
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- NA
         }
         } else{
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- NA
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- NA
         }
      }
   }

# Row and column names (lon and lat).
   colnames(mat_list_annual_trend_node) <- (matrix(lon_range_node,nrow=res)[1,] + (res/2))
   rownames(mat_list_annual_trend_node) <- (matrix(lat_range,nrow=res)[1,] + (res/2))

# Gather all the output into a single array, along longitude.
   mat_list_annual_trend <- do.call( cbind, allgather( mat_list_annual_trend_node ) )

# Create data structure including metadata.
   array_meta <- list(dataset_name="CCI L4",band="unknown",
                      years=lab_years,months=lab_months,year_centre=lab_y_centre,resolution=res,
                      orig_lat_cell=(vec_lat-0.5), orig_lon_cell=(vec_lon-0.5),
                      orig_lat_mid=vec_lat, orig_lon_mid=vec_lon,
                      lat_cell=matrix(lat_range,nrow=res)[1,],lon_cell=matrix((lon_range-180),nrow=res)[1,],
                      lat_mid=(matrix(lat_range,nrow=res)[1,] + (res/2)),lon_mid=(matrix((lon_range-180),nrow=res)[1,] + (res/2)),
                      trend_stats="swh_mean",trend=c("slope","P-val"))
   list_CCI_trend <- list(array_meta,mat_list_annual_trend)

# Write out data.
   data_file <- paste("./output/",res,"deg/list_trend_",lab_years,"_",lab_months,"_",lab_y_centre,"_",flag_reg,".Robj",sep="")
   save(list_CCI_trend,file = data_file)

   finalize()


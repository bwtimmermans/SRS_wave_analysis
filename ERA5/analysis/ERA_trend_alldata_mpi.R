# source("/home/ben/research/NOC/SRS_wave_analysis/ERA5/analysis/ERA_trend_alldata_mpi.R")

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
   lon_range <- c(0,359.5)
   lat_range <- c(-71.5,72)

# Actual domain.
   nc1 = nc_open("/backup/datasets/ERA5/era5_swh_1979-2018.nc")
   vec_lon <- ncvar_get(nc1,"longitude")
   vec_lat <- ncvar_get(nc1,"latitude")
   vec_time <- ncvar_get(nc1,"time")
   nc_close(nc1)

# Resolution (1 = 1 degree, 2 = 2 degree, etc).
   res <- 4

# Mid points for meta data.
   lon_mid <- matrix((lon_range[1]:lon_range[2]),nrow=res)[1,] + (res/2) - 1/4
   lat_mid <- rev( matrix(lat_range[1]:lat_range[2],nrow=res)[1,] + (res/2) - 1/4 )

# Find indices for ERA longitude, latitude and time.
   lon_idx_range <- which(lon_range[1] == vec_lon):which(lon_range[2] == vec_lon)
   lat_idx_range <- rev(which(lat_range[1] == vec_lat):which(lat_range[2] == vec_lat))

# Analysis duration.
# Note for 'flag_winter', first year must be >= 1993.
   anal_years <- 1992:2000
   anal_years <- 2001:2009
   anal_years <- 2010:2018
   anal_years <- 1992:2018
   anal_years <- 1980:2018

   time_idx_range <- c(which(as.numeric(format(as.POSIXct(vec_time*3600, origin = '1900-01-01', tz='GMT'),'%Y')) == anal_years[1])[1]:
                       which(as.numeric(format(as.POSIXct(vec_time*3600, origin = '1900-01-01', tz='GMT'),'%Y')) == anal_years[length(anal_years)])[12])

# Flag for complete winter season.
   flag_winter <- FALSE
   if (flag_winter) { lab_y_centre = "n_winter" } else { lab_y_centre = "n_summer" }

   lab_years <- paste(anal_years[c(1,length(anal_years))],collapse='-')
   all_months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

# Months for analysis.
   flag_annual <- TRUE
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
   flag_reg <- "ONI"
   flag_reg <- "NAO"
   flag_reg <- "none"

# Parallelise over the geographic range, by longitude.
# Parallelise over longitude, divide the range by number of processing cores.
   dim.lon <- length(lon_idx_range)
   lon.base <- floor(dim.lon / .comm.size)
   #vec.lon_dim <- rep( lon.base, .comm.size ) + c( rep(1,dim.lon %% .comm.size), rep(0,(.comm.size - dim.lon %% .comm.size)) )
   vec.lon_dim <- (2*res) * ( floor( ( dim.lon / (2*res) ) / .comm.size) + c( rep(1,( dim.lon / (2*res) ) %% .comm.size), rep(0,(.comm.size - ( dim.lon / (2*res) ) %% .comm.size)) ) )

   lon.start.idx <- sum( vec.lon_dim[c(0:.comm.rank)] ) + 1
   lon.dim.node <- vec.lon_dim[.comm.rank+1]

   lon_range_node <- lon_idx_range[lon.start.idx:(lon.start.idx + lon.dim.node - 1)]
   lon.start.idx.nc <- lon_range_node[1]
   print(paste("lon_range_node:",lon_range_node))

# Latitude.
   lat.start.idx <- which(vec_lat == lat_range[2])
   lat.stop.idx <- which(vec_lat == lat_range[1])
   dim.lat <- length(lat_idx_range)

# Grid indices for aggregation.
   mat_lat_grid_idx <- matrix(1:length(lat_idx_range),nrow=2*res)
   mat_lon_grid_idx <- matrix(1:length(lon_range_node),nrow=2*res)
   #mat_lon_grid_idx <- matrix(lon.start.idx:(lon.start.idx+length(lon_range_node)-1),nrow=res)

# ================================================================= #
# Specify data files and load all data (only summary statistics, so do all in one).

# Data path.
   data_path <- "/backup/datasets/ERA5/era5_swh_1979-2018.nc"

# Load data.
   nc1 = nc_open(data_path)
   array_ERA_temp <- ncvar_get(nc1,"swh",start=c(lon.start.idx.nc,lat.start.idx,time_idx_range[1]), count=c((lon.dim.node),dim.lat,(12*length(anal_years))))
   nc_close(nc1)
# Reformat array to separate years.
   array_ERA <- array(array_ERA_temp,dim=c(dim(array_ERA_temp)[1],dim(array_ERA_temp)[2],12,length(anal_years)))

# Clean up missing values.
#   array_ERA[array_ERA > 50] <- NA

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

# ================================================================= #
# Data structures.
   mat_list_annual_stats_node <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
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
# Test for at least 15% observations.
         if ( (sum(!is.na(as.vector(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],1,1]))) /
              length(as.vector(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],1,1]))) > 0.15 ) {
# Centre year on winter season.
            if (flag_winter) {
               if (flag_annual) {
# Loop over statistics for aggregation.
                  for (y_idx in 2:length(anal_years)) {
                     mat_stats_temp[y_idx,1] <- sum(
                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],7:12,(y_idx-1)],
                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],1:6,y_idx],na.rm=T
                                                   ) 
#                     mat_stats_temp[y_idx,2] <- mean(
#                                                     c(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),7:12,2],
#                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,1:6,2]),na.rm=T
#                                                    ) 
#                     mat_stats_temp[y_idx,3] <- mean(
#                                                     c(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),7:12,3],
#                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,1:6,3]),na.rm=T
#                                                    ) 
                  }
#                  mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
               } else {
                  for (y_idx in 2:length(anal_years)) {
                     mat_stats_temp[y_idx,1] <- sum(
                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],list_split_months[[1]],(y_idx-1)],
                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],list_split_months[[2]],y_idx],na.rm=T
                                                   ) 
#                     mat_stats_temp[y_idx,2] <- mean(
#                                                     c(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),list_split_months[[1]],2],
#                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,list_split_months[[2]],2]),na.rm=T
#                                                    ) 
#                     mat_stats_temp[y_idx,3] <- mean(
#                                                     c(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],(y_idx-1),list_split_months[[1]],3],
#                                                     array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],y_idx,list_split_months[[2]],3]),na.rm=T
#                                                    ) 
                  }
#                  mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
               }
            } else {
# Centre year on summer (normal).
               if (flag_annual) {
# Loop over statistics for aggregation.
                  mat_stats_temp[,1] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],,x],na.rm=T) } )
#                  mat_stats_temp[,2] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,,2],na.rm=T) } )
#                  mat_stats_temp[,3] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,,3],na.rm=T) } )
#                  mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
               } else {
                  mat_stats_temp[,1] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],as.numeric(anal_months),x],na.rm=T) } )
#                  mat_stats_temp[,2] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,as.numeric(anal_months),2],na.rm=T) } )
#                  mat_stats_temp[,3] <- sapply(X=1:length(anal_years),FUN=function(x) { mean(array_ERA[mat_lon_grid_idx[,lon_res_idx],mat_lat_grid_idx[,lat_res_idx],x,as.numeric(anal_months),3],na.rm=T) } )
#                  mat_stats_temp[,4] <- mat_stats_temp[,3] - mat_stats_temp[,2]^2
               }
            }
            df_stats <- as.data.frame(mat_stats_temp)
            colnames(df_stats) <- c("swh_mean","empty","empty","empty")
            rownames(df_stats) <- anal_years

# Trend (linear regression).
            df_Q <- data.frame(cbind(year=anal_years,df_stats,mat_ONI[ONI_years,],mat_NAO[NAO_years,]))

            mat_trend <- matrix(NA,nrow=length(vec_q),ncol=6)
            for (qq in 1:1) {
# Catch if lack of data for regression.
               if (all(!is.nan(df_Q[,(qq+1)]))) {
# Mean and var.
                  mat_trend[qq,5] <- mean(df_Q[,(qq+1)])
                  mat_trend[qq,6] <- var(df_Q[,(qq+1)])
# Trend.
                  if (flag_reg == "NAO") {
                     lm_Q <- lm(df_Q[,(qq+1)] ~ year + NAO_mean,data=df_Q)
                  } else if (flag_reg == "ONI") {
                     lm_Q <- lm(df_Q[,(qq+1)] ~ year + ONI_mean,data=df_Q)
                  } else {
                     lm_Q <- lm(df_Q[,(qq+1)] ~ year,data=df_Q)
                     #lm_Q <- lm(eval(parse(text=paste(colnames(mat_b_q)[qq]))) ~ year,data=df_Q)
                  }
                  sum_lm_Q <- summary(lm_Q)
                  mat_trend[qq,1:4] <- c(sum_lm_Q$coefficients[2,1],sum_lm_Q$coefficients[2,4],sum_lm_Q$coefficients[2,2],sum_lm_Q$sigma)
               }
            }
            colnames(mat_trend) <- c("year_slope","Pr(>|t|)","trend_SE","reg_SE","mean","var")
            if (!all(is.na(mat_trend))) {
# Store for writing.
               mat_list_annual_stats_node[[lat_res_idx,lon_res_idx]] <- df_stats
               mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- mat_trend
            } else{
               mat_list_annual_stats_node[[lat_res_idx,lon_res_idx]] <- NA
               mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- NA
            }
         } else{
            mat_list_annual_stats_node[[lat_res_idx,lon_res_idx]] <- NA
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- NA
         }
      }
   }

# Gather all the output into a single array, along longitude.
   mat_list_annual_trend <- do.call( cbind, allgather( mat_list_annual_trend_node ) )
   mat_list_annual_stats <- do.call( cbind, allgather( mat_list_annual_stats_node ) )

# Create data structure including metadata.
   array_meta <- list(dataset_name="ERA5",band="NA",
                      years=lab_years,months=lab_months,year_centre=lab_y_centre,resolution=res,
                      orig_lat_cell=vec_lat, orig_lon_cell=vec_lon,
                      orig_lat_mid=vec_lat, orig_lon_mid=vec_lon,
                      lat_mid=lat_mid,lon_mid=lon_mid,
                      trend_stats="swh_mean",trend=c("slope","P-val","trend_SE","reg_SE","mean","var"))
# Trends.
   list_ERA_trend <- list(array_meta,mat_list_annual_trend)
# Stats.
   list_ERA_stats <- list(array_meta,mat_list_annual_stats)

# Write out data.
# Trend.
   data_file <- paste("./output/",res,"deg/list_trend_",lab_years,"_",lab_months,"_",lab_y_centre,"_",flag_reg,"_mean.Robj",sep="")
   save(list_ERA_trend,file = data_file)
# Stats.
   data_file <- paste("./output/",res,"deg/list_stats_",lab_years,"_",lab_months,"_",lab_y_centre,"_",flag_reg,"_mean.Robj",sep="")
   save(list_ERA_stats,file = data_file)

   finalize()


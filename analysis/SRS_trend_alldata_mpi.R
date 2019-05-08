# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/SRS_trend_alldata.R")

# Script to load multiple datasets from IMOS and obtain summary statistics and annual trend.
# BT 04/2019

# Libraries.
   library(ncdf4)
   library(abind)
   library(extRemes)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")
# MPI.
   library(pbdMPI, quietly = TRUE)
   .comm.size <- comm.size()
   .comm.rank <- comm.rank()

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

# Parallelise over the geographic range, by longitude.
# Parallelise over longitude, divide the range by number of processing cores.
   dim.lon <- length(lon_range)
   lon.base <- floor(dim.lon / .comm.size)
   #vec.lon_dim <- rep( lon.base, .comm.size ) + c( rep(1,dim.lon %% .comm.size), rep(0,(.comm.size - dim.lon %% .comm.size)) )
   vec.lon_dim <- res * ( floor( ( dim.lon / res ) / .comm.size) + c( rep(1,( dim.lon / res ) %% .comm.size), rep(0,(.comm.size - ( dim.lon / res ) %% .comm.size)) ) )

   lon.start.idx <- sum( vec.lon_dim[c(0:.comm.rank)] ) + 1
   lon.dim.node <- vec.lon_dim[.comm.rank+1]

   lon_range_node <- lon_range[lon.start.idx:(lon.start.idx + lon.dim.node - 1)]
   print(paste("lon_range_node:",lon_range_node))

# Grid indices for aggregation.
   mat_lat_grid_idx <- matrix(1:length(lat_range),nrow=res)
   mat_lon_grid_idx <- matrix(1:length(lon_range_node),nrow=res)

# Data set selection.
   mission_idx <- 2:10

# Months for analysis.
   flag_annual <- TRUE
   #anal_months <- c("01","02","03")
   #anal_months <- c("04","05","06")
   #anal_months <- c("07","08","09")
   #anal_months <- c("10","11","12")
   anal_months <- c("01","02","11","12")
   #anal_months <- c("01","02","03","10","11","12")
   #anal_months <- c("04","05","06","07","08","09")
   #anal_months <- c("01","02","03","04","05","06")
   #anal_months <- c("07","08","09","10","11","12")
   if (!flag_annual) {
      lab_months <- paste(c("J","F","M","A","M","J","J","A","S","O","N","D")[as.numeric(anal_months)],collapse='')
   } else {
      lab_months <- "annual"
   }

   #anal_months <- c("04","05","06","07","08","09")

# ================================================================= #
# Data path.
   data_path <- "/home/ben/research/NOC/SRS_wave_analysis/datasets/IMOS/"
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   mat_valid_dates <- cbind(c(1986,1990),c(1992,1995),c(1993,2004),c(1996,2008),c(1999,2007),c(2002,2012),c(2009,2017),c(2011,2017),c(2012,2017),c(2014,2017),c(2016,2017),c(2017,2017))

# File base name.
   #array_filenames <- array(NA,dim=c(length(lat_range),length(lon_range_node),length(mission_idx)))
   array_filenames <- NULL
   for (m_idx in mission_idx) {
      nc_sat1 <- paste(vec_datasets[m_idx],"/IMOS_SRS-Surface-Waves_MW_",vec_datasets[m_idx],"_FV02_",sep="")
      vec_lat_sat <- unlist( lapply(X=rev(lat_range),FUN=function(x) { paste(tail(strsplit(paste("00",x,"N",sep=""),split='')[[1]],n=4),collapse='') } ) )
      vec_lon_sat <- paste(lon_range_node,"E",sep="")
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
   mat_list_annual_trend_node <- matrix(list(),nrow=dim(mat_lat_grid_idx)[2],ncol=dim(mat_lon_grid_idx)[2])
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
                     print(paste("LAT:",lat_range[lat_idx],"; LON:",lon_range_node[lon_idx]))
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
               if (!flag_annual) {
                  nc1_date_temp <- nc1_date_temp[ format(nc1_date_temp,"%m") %in% anal_months ]
               }
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
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- mat_trend
         } else{
            mat_list_annual_KU[[lat_res_idx,lon_res_idx]] <- NA
            mat_list_annual_trend_node[[lat_res_idx,lon_res_idx]] <- NA
         }
      }
   }

# Gather all the output into a single array, along longitude.
   source("/home/ben/research/code/R/func_abind/func_abind.R")
   #array.chi <- do.call( func_abind, list( allgather( mat_list_annual_trend_node ), along=1 ) )
   mat_list_annual_trend <- do.call( cbind, allgather( mat_list_annual_trend_node ) )

# Create data structure including metadata.
   array_meta <- list(mission_idx=mission_idx,data_name=vec_datasets[mission_idx],time_period=lab_months,
                      orig_lat_cell=lat_range, orig_lon_cell=lon_range,
                      lat_cell=matrix(lat_range,nrow=res)[1,],lon_cell=matrix(lon_range,nrow=res)[1,],
                      lat_mid=(matrix(lat_range,nrow=res)[1,] + (res/2)),lon_mid=(matrix(lon_range,nrow=res)[1,] + (res/2)),
                      trend_stats=paste("Q",100*c(0.5,0.9,0.95),sep=""),trend=c("slope","P-val"))
   list_SRS_trend <- list(array_meta,mat_list_annual_trend)

# Write out data.
   data_file <- paste("./output/mpi_test/list_",paste(vec_datasets[mission_idx],collapse='_'),"_all_data_trend_",lab_months,".Robj",sep="")
   save(list_SRS_trend,file = data_file)

   finalize()


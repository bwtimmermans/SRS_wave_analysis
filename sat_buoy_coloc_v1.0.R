# source("/home/ben/research/NOC/waves/sat_buoy_coloc_v1.0.R")

# Librarys.
   library(ncdf4)

# Ensemble members.
   ens_size <- 150

# Experiment name.
   exp_home <- paste("/scratch2/scratchdirs/timmer/hurricanes/ike/4.5km/X_ike_ens4/",1:ens_size,"/output/run1",sep="")
   exp_home_hist <- paste("/scratch2/scratchdirs/timmer/hurricanes/ike/4.5km/X_ike_nl1/output/run1",sep="")

# Simualtion length (days).
   sim_len <- 8

# Time steps per day in NetCDF.
   ncdf_step <- 96

# Buoys.
   #buoy_list <- c("42001","42002","42036","42039","42040","42055")
   buoy_list <- c("46066")
   vec_t_time <- numeric(length(buoy_list))
   buoy_46066_loc <- c((360-155.047),52.785)

# Observation per hour.
   buoy_obs_freq <- c(1,1,1,1,1,1,1,1)
   vec_buoy_minute <- c(50,50,50,50,50,50)
   mfrow_form <- cbind(c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),c(3,3),c(3,3))

#-----------------------------------------------------------------------#
# Observed data.
#-----------------------------------------------------------------------#

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

#-----------------------------------------------------------------------#
# Find colocations.
#-----------------------------------------------------------------------#

# Loop over SAT observations (time mid point).
  nc1_breaks_buoy_time <- 0
  for (i in 1:length(vec_nc1_break_mid)) {
  #for (i in 1:100) {
      print(paste("i:",i))

# Find the required date string (from SAT data.)
      SAT_year <- substr(nc1_date[vec_nc1_break_mid[i]],start=1,stop=4)
      SAT_month <- substr(nc1_date[vec_nc1_break_mid[i]],start=6,stop=7)
      SAT_day <- substr(nc1_date[vec_nc1_break_mid[i]],start=9,stop=10)
      SAT_hour <- substr(nc1_date[vec_nc1_break_mid[i]],start=12,stop=13)
      #SAT_minute <- vec_buoy_minute[b.idx]
      SAT_minute <- "00"

# Start a counter.
      test_SAT_date <- strptime(paste(SAT_year,"-",SAT_month,"-",SAT_day," ",SAT_hour,":",SAT_minute,sep=""),format="%Y-%m-%d %H:%M",tz="GMT")
      print(paste("SAT date:",test_SAT_date))
# Set start time to the first entry of relevant year (to save time).
      t_time <- which(array_buoy_obs[,1,1] == SAT_year)[1]
# If there are no matches per year, skip.
      if (!is.na(t_time)) {

      if (array_buoy_obs[t_time,2,1] <= as.numeric(SAT_month)) {

         #test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":",array_buoy_obs[t_time,5,b.idx],":00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
         test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
         print(paste("Obs date1:",test_obs_date))

         while ( test_obs_date != test_SAT_date & t_time != 94934 ) {
            t_time <- t_time + 1
            #test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":",array_buoy_obs[t_time,5,b.idx],":00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
            test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
#print(paste("Obs date1:",test_obs_date))
         }
         print(paste("Obs date2:",test_obs_date))

         #if ( test_obs_date == "NA" ) {
         if ( t_time == 94934 ) {
            nc1_breaks_buoy_time[i] <- NA
         } else {
            nc1_breaks_buoy_time[i] <- t_time
         }
      } else {
         nc1_breaks_buoy_time[i] <- NA
      }
   } else {
      nc1_breaks_buoy_time[i] <- NA
   }
   }

#=======================================================================#
# Plotting.
#-----------------------------------------------------------------------#
# 2.a Plot wave stats (Hs, Tp).
#-----------------------------------------------------------------------#
   cex_lab <- 1.2
   X11()
   #pdf("/global/homes/t/timmer/waves/hurricanes/hurricanes_waves_OC/figures/ike_hs_tp_1buoy_samples_OSM_1.pdf", width=10, height=13)
   #pdf("./test.pdf", width=10, height=12)
   par(mfrow=mfrow_form[,length(buoy_list)],oma=c(0,0,0,0),mgp=c(3,1,0),mar=c(4,4,4,5))

#   for (b.idx in 1)
   for (b.idx in 1:length(buoy_list))
   {
      plot_title <- paste("NDBC", buoy_list[b.idx])
      plot(NULL,xlim=wind_time[c(time_start,(time_start+sim_len*ncdf_step))],ylim=c(-9,15),xaxt='n',yaxt='n',xlab=NA,ylab=NA,main=plot_title,cex.main=1.5*cex_lab)
      axis(side=1,at=time_labels,labels=as.POSIXct(time_labels*3600, origin = '1900-01-01', tz='GMT'))
      #abline(v=time_labels,lty=4)
      axis(side=2,at=seq(0,12,2),labels=seq(0,12,2),cex.axis=cex_lab)
      mtext(expression(H[s]~(m)),side=2,line=2,cex=cex_lab)
      abline(h=0,lty=1)
      abline(h=c(5,10),lty=2)

# Find the required date string (from SAT data.)
      ww3_year <- substr(nc1_date[800],start=1,stop=4)
      ww3_month <- substr(nc1_date[800],start=6,stop=7)
      ww3_day <- substr(nc1_date[800],start=9,stop=10)
      ww3_hour <- substr(nc1_date[800],start=12,stop=13)
      #ww3_minute <- array_ww[4,3,b.idx,1]
      #ww3_minute <- vec_buoy_minute[b.idx]
      ww3_minute <- "00"

# Start a counter.
      test_ww3_date <- strptime(paste(ww3_year,"-",ww3_month,"-",ww3_day," ",ww3_hour,":",ww3_minute,sep=""),format="%Y-%m-%d %H:%M",tz="GMT")
print(paste("WW3 date:",test_ww3_date))
      t_time <- 4
      #test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":",array_buoy_obs[t_time,5,b.idx],":00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
      test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
print(paste("Obs date1:",test_obs_date))

      while ( test_obs_date != test_ww3_date ) {
         t_time <- t_time + 1
         #test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":",array_buoy_obs[t_time,5,b.idx],":00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
         test_obs_date <- strptime(paste(array_buoy_obs[t_time,1,b.idx],"-",array_buoy_obs[t_time,2,b.idx],"-",array_buoy_obs[t_time,3,b.idx]," ",array_buoy_obs[t_time,4,b.idx],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT")
print(paste("Obs date1:",test_obs_date))
            if ( b.idx == 8 ) {
	          t_time <- 3201
                  test_obs_date <- 1
                  test_ww3_date <- 1
               }
          }
print(paste("Obs date2:",test_obs_date))

# Capture t_time.
      vec_t_time[b.idx] <- t_time

# Observed data.
      obs_plot_seq <- seq(wind_time[time_start],,(1/buoy_obs_freq[b.idx]),(sim_len*24*buoy_obs_freq[b.idx]))
      obs_plot_data <- array_buoy_obs[( t_time:(t_time - 1 + sim_len*24*buoy_obs_freq[b.idx]) ),9,b.idx]
      print(paste("T_Time:",t_time))

# Print diagnostics.
#      print(paste("Buoy:", buoy_list[b.idx],"Start time:",as.POSIXct(min(obs_plot_seq)*3600, origin = '1900-01-01', tz='GMT')))
#      print(paste("Buoy:", buoy_list[b.idx],"Stop time:",as.POSIXct(max(obs_plot_seq)*3600, origin = '1900-01-01', tz='GMT')))
#      print(paste("Data:", array_buoy_obs[1:2,,b.idx]))
#      print(paste("Data:", array_buoy_obs[1:2,,b.idx]))
# WW3 Hs data.
   #for (m.idx in 1:ens_size) {
   for (m.idx in seq(1,ens_size,4)) {
      lines(seq(wind_time[time_start],wind_time[sim_len*ncdf_step],0.25), array_ww[time_start:(sim_len*ncdf_step),5,b.idx,m.idx],col="darkgrey",lwd=1.0)
   }
# Case hist.
      lines(seq(wind_time[time_start],wind_time[sim_len*ncdf_step],0.25), array_ww_hist[time_start:(sim_len*ncdf_step),5,b.idx],col="darkblue",lwd=2.5)
# Plot obs. points after lines.
      points(obs_plot_seq, obs_plot_data, pch=1, cex=0.8, col="black")
#   }
#dev.off()
# Highlight points.
#   points(wind_time[770], mat_ww_41010[32,5],col="darkred",pch=4,cex=3.0)
#   points(wind_time[770], mat_ww_41010[32,5],col="darkred",pch=1,cex=2.0)
# Boxplot of MC data.
#   Hs_MC_output <- read.table("../R_Hs_max/Hs_MC_output")
#   boxplot(Hs_MC_output,add=TRUE,at=wind_time[770],pars = list(boxwex = 12000),yaxt='n',outline=TRUE)

# Plot Tp data.
   par(new=TRUE)
   plot(NULL,xlim=wind_time[mat_xlim],ylim=c(2.0,40),xlab=NA,ylab=NA,xaxt='n',yaxt='n')
   mtext(expression(T[p]~(s)),side=4,line=3,cex=cex_lab)
   axis(side=4,at=seq(5,20,5))
   #abline(h=5,lty=2)
   #points(seq(wind_time[1],,(1/buoy_obs_freq[b.idx]),(14*24*buoy_obs_freq[b.idx])),array_buoy_obs[1:(14*24*buoy_obs_freq[b.idx]),10,b.idx],pch=3,cex=0.8,col="darkblue")

   obs_plot_data_tp <- array_buoy_obs[( t_time:(t_time - 1 + sim_len*24*buoy_obs_freq[b.idx]) ),10,b.idx]

# WW3 Tp data.
   for (m.idx in seq(1,ens_size,4)) {
      lines(seq(wind_time[time_start],wind_time[sim_len*ncdf_step],0.25), 1/array_ww[time_start:(sim_len*ncdf_step),10,b.idx,m.idx],col="darkgrey",lwd=1.0)
#      lines(wind_time[sim_start_idx:sim_stop_idx],1/mat_ww_41010_sample[,10,i],col="blue",lwd=0.5)
   }
#   lines(wind_time[sim_start_idx:sim_stop_idx],1/array_ww[seq(1,,2,53),10,b.idx],col="darkblue",lwd=2.5)
   lines(seq(wind_time[time_start],wind_time[sim_len*ncdf_step],0.25), 1/array_ww_hist[time_start:(sim_len*ncdf_step),10,b.idx],col="darkblue",lwd=2.5,lty=2)
# Plot obs. points after lines.
   points(obs_plot_seq, obs_plot_data_tp, pch=3, cex=0.8, col="black")

## Boxplot of MC data.
##   Tp_MC_output <- read.table("../R_Tp_max/Tp_MC_output")
##   boxplot(Tp_MC_output,add=TRUE,at=wind_time[770],pars = list(boxwex = 12000),yaxt='n')
#
# Legend.
      legend(
             x=wind_time[floor(mat_xlim[2]/40)],y=40,
             c(expression(H[s]~(Obs.)),expression(H[s]~(Hist.~4~km)),expression(H[s]~samples~(4~km)),expression(T[p]~(Obs.)),expression(T[p]~(Hist.~4~km)),expression(T[p]~samples~(4~km))),
             lty=c(-1,1,1,-1,2,1),lwd=c(-1,2.5,1.5,-1,2.5,1.5),pch=c(1,-1,-1,3,-1,-1),col=c("black","darkblue","darkgrey","black","darkblue","darkgrey"),bty="o",bg="white",
             x.intersp=0.7,y.intersp=1.0,cex=1.1,xjust=0.1
            )

   }
   #dev.off()

#-----------------------------------------------------------------------#
# 2.b Scatter plot of obs. vs simulated data.
#-----------------------------------------------------------------------#
#   X11()
#   #pdf("/home/ben/research/code/R/ww3/global/exp_6_ecmwf_2011_JAN_025_loSpec/scatter_ecmwf_025_loSpec_buoys.pdf", paper="USr",width=12, height=12)
#   par(mfrow=c(2,2),oma=c(0,0,0,0),mgp=c(3,1,0),mar=c(4,4,4,5)+0.1)
#
#   for (b.idx in 1:length(buoy_list)) {
#         plot(array_buoy_obs[seq(vec_t_time[b.idx],,buoy_obs_freq[b.idx],(mat_ww_dims[1,1]/4)),9,b.idx],array_ww[seq(1,,4,(mat_ww_dims[1,1]/4)),5,b.idx],xlim=c(0,12),ylim=c(0,12),xlab=paste("Observation at NDBC",buoy_list[b.idx]),ylab="WW3 data")
#         abline(a=0,b=1)
#      }
#   #dev.off()


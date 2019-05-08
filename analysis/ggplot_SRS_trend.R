# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/ggplot_SRS_trend.R")
   library(ggplot2)
   library(grid)
   library(gridExtra)
   library(RColorBrewer)
   library(colorspace)
   library(fBasics)
   library(abind)

#-------------------------------------------------------#
# Missions.
   mission_idx <- 2:10
   lab_months <- "JFMAMJ"
   lab_months <- "JASOND"
   lab_months <- "JFND"
   lab_months <- "annual"

# Read analysis data (amtrix of lists).
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   #data_path <- paste("./output/test_block/list_",paste(vec_datasets[mission_idx],collapse='_'),"_trend.Robj",sep="")
   #data_path <- paste("./output/test_block/list_",paste(vec_datasets[mission_idx],collapse='_'),"_all_data_trend.Robj",sep="")
   #data_path <- paste("./output/all_KU/list_",paste(vec_datasets[mission_idx],collapse='_'),"_all_data_trend_",lab_months,".Robj",sep="")
# MPI.
   data_path <- paste("./output/mpi_test/list_",paste(vec_datasets[mission_idx],collapse='_'),"_all_data_trend_",lab_months,".Robj",sep="")

# Matrix to hold datasets.
   attach(data_path[1])
   attached_data <- list_SRS_trend
   detach(pos=2)

# Set up data structures.
   mat_plot_data <- matrix(NA,nrow=length(attached_data[[1]]$lat_cell),ncol=length(attached_data[[1]]$lon_cell))
   mat_plot_CI <- matrix(NA,nrow=length(attached_data[[1]]$lat_cell),ncol=length(attached_data[[1]]$lon_cell))
   lab_dataset <- attached_data[[1]]$data_name
   lat_mid <- attached_data[[1]]$lat_mid
   lon_mid <- attached_data[[1]]$lon_mid
   lab_stats <- attached_data[[1]]$trend_stats
   lab_period <- attached_data[[1]]$time_period

# Analysis data.
# 1: Q50
# 2: Q90
# 3: Q95
   stat_idx <- 2
   list_data <- attached_data[[2]]

# list_SRS_trend dimensions:
# 1,2: lat, lon
# Matrix.
# 3: Q50,Q90,Q95
# 4: year_slope, P-val
   for (lat_idx in 1:dim(mat_plot_data)[1]) {
     for (lon_idx in 1:dim(mat_plot_data)[2]) {
# Temporal trend.
         if ( !is.na(list_data[lat_idx,lon_idx]) ) { mat_plot_data[lat_idx,lon_idx] <- list_data[[lat_idx,lon_idx]][[stat_idx,1]] }
         #if ( ! all(is.na(array.chi[[j,i,1,1]])) ) { array_qq[j,i,k] <- mat_qq_out[[j,i]]$regression[[1]][2]
# Trend CI.
         if ( !is.na(list_data[lat_idx,lon_idx]) ) { mat_plot_CI[lat_idx,lon_idx] <- list_data[[lat_idx,lon_idx]][[stat_idx,2]] }
      }
   }
# Fix for CI info.
   mat_plot_CI[mat_plot_CI > 0.1] <- NA
   mat_plot_CI[!is.na(mat_plot_CI)] <- 1
   #mat_plot_CI[is.na(mat_plot_CI)] <- 0

#-------------------------------------------------------#
# Set up data frame for ggplot.
   df_plot <- NULL
   plot_labels <- paste(paste(lab_dataset,collapse=', '),": Temporal trend in ",lab_stats[stat_idx]," (",lab_period," KU-band)",sep="")
   for (k in 1:1) {
      df_plot <- rbind( df_plot,
                            cbind( expand.grid( lat=rev(lat_mid), lon=lon_mid ), plot_stat=as.vector(mat_plot_data), plot_CI=as.vector(mat_plot_CI), anal=plot_labels[k] )
                          )
   }

# Custom 3.
   #colfunc <- colorRampPalette(c("white","khaki1","orangered","purple","black"))
   #colfunc <- colorRampPalette(rev(c("red","white","blue")))
   colfunc <- colorRampPalette(rev(c("firebrick2","white","royalblue")))
   seq.plot_cols <- colfunc(20)

# Colour ranges.
# Colour ranges for counts.
   if ( stat_idx == 1) {
      #plot_breaks_lo <- 0.9*min(df_plot$plot_stat,na.rm=T)
      #plot_breaks_hi <- 1.1*max(df_plot$plot_stat,na.rm=T)
      #col_breaks <- seq(plot_breaks_lo,plot_breaks_hi,500)
      col_breaks <- seq(-0.04,0.04,0.01)
      #plot_lims <- c(0.95*min(df_plot$plot_stat,na.rm=T),1.05*max(df_plot$plot_stat,na.rm=T))
      plot_lims <- c(-0.04,0.04)
   } else if ( stat_idx == 2 | stat_idx == 3 ) {
      col_breaks <- seq(-0.07,0.07,0.01)
      plot_lims <- c(-0.075,0.075)
   }

# Chi plot.
   p1 <- ggplot(data = df_plot) + #scale_x_continuous(limits = c(-100, -50), expand = c(0,0)) + scale_y_continuous(limits = c(10, 50), expand = c(0,0)) +
         ylab("Latitude\n") + xlab("\nLongitude") +
         #coord_map(xlim = range(lon_mid), ylim = range(lat_mid), projection = "lambert", lat0 = 33, lat1 = 45) +
         coord_map(xlim = c(105,275), ylim = c(0,55), projection = "lambert", lat0 = 10, lat1 = 40) +
         #coord_map(projection = "rectangular", lat0 = 0, xlim = c(0,360), ylim = c(-75,75)) +
         #coord_map(projection = "rectangular", lat0 = 0, xlim = c(120,310), ylim = c(-50,65)) +
#  scale_x_continuous(limits = c(-180, 180)) +
#  scale_y_continuous(limits = c(-75, 75)) +
# Fill.
         geom_tile(aes(x = lon, y = lat, fill = plot_stat)) +
	 #geom_point(size = 1, shape = 15 ) +
         scale_fill_gradientn( colours = seq.plot_cols,
                                #name = expression("SD(" * hat(chi) * "*)"),
                                name = paste(lab_stats[stat_idx],"trend (m/year)"),
                                limits = plot_lims,
                                breaks=col_breaks,
                                na.value = "transparent",
                                guide = guide_colorbar(title.position = "bottom",label.position = "left",ticks.colour = "black",ticks.linewidth = 8.0) ) +
# CI information.
         geom_point(aes(x = lon, y = lat, colour = plot_CI), size = 3, shape = 16) +
         scale_colour_gradientn(colors = c("black","black"), values = c("0.0","1.0"), guide = FALSE, na.value = "transparent") +
# Map.
         geom_polygon(data = map_data("world2"), aes(x = long, y = lat, group = group), color = "#000000", fill = NA, size = 0.55) +
# Multiple plots.
         facet_wrap(~anal, ncol=2) +
         #ggtitle(plot_labels) +
# Theme stuff.
         theme(axis.title.x=element_blank(),
               axis.title.y=element_blank(),

               strip.text = element_text(size = 50, margin = margin(25,0,25,0)),
               strip.background = element_rect(fill = "white"),
               panel.spacing.x = unit(1, "lines"),
               panel.spacing.y = unit(2, "lines"),
               axis.text = element_blank(),
               axis.ticks = element_blank(),

               legend.position = "left",
               legend.margin = margin(0,50,0,0),
               legend.key.width = unit(1, "inch"),
               legend.key.height = unit(2, "inch"),
               legend.title = element_text(size = 30, margin = margin(25,0,0,0)),
               legend.title.align = 0.5,
               legend.text = element_text(size = 25, margin = margin(0,0,0,25))
         )

# Include blank panels.
   #fig_file_name <- paste("./figures/SRS_trends/trend_CI_block_",lab_stats[stat_idx],"_",paste(lab_dataset,collapse='_'),".png",sep="")
   fig_file_name <- paste("./figures/SRS_trends/trend_CI_all_KU_",lab_months,"_",lab_stats[stat_idx],"_",paste(lab_dataset,collapse='_'),"_mpi.png",sep="")
   png(filename = fig_file_name, width = 2400, height = 1300)
   plot(p1)
   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))


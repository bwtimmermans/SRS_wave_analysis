# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/ggplot_stats.R")
   library(ggplot2)
   library(grid)
   library(gridExtra)
   library(RColorBrewer)
   library(colorspace)
   library(fBasics)
   library(abind)

#-------------------------------------------------------#
# Read analysis data (amtrix of lists).
   vec_datasets <- c("GEOSAT","ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","HY-2","SARAL","JASON-3","SENTINEL-3A")
   #data_path <- paste("./test_output/list_IMOS_stats.Robj",sep="")
   data_path <- paste("./test_output/list_",vec_datasets[4],"_stats.Robj",sep="")

# Matrix to hold datasets.
   attach(data_path[1])
   mat_plot_data <- matrix(NA,nrow=length(list_IMOS_stats[[1]]$lat),ncol=length(list_IMOS_stats[[1]]$lon))
   lab_dataset <- list_IMOS_stats[[1]]$data_name
   lat_mid <- list_IMOS_stats[[1]]$lat
   lon_mid <- list_IMOS_stats[[1]]$lon
   #detach(pos=2)

# array_stats dimensions:
# 1,2: lat, lon
# 3: data set index
# 4:
#      for (j in 1:dim(mat_plot_data)[1]) {
#         for (i in 1:dim(mat_plot_data)[2]) {
#            if ( ! all(is.na(mat_qq_out[[j,i]])) ) { mat_plot_data[j,i] <- mat_qq_out[[j,i]]$regression[[1]][2] }
#            #if ( ! all(is.na(array.chi[[j,i,1,1]])) ) { array_qq[j,i,k] <- mat_qq_out[[j,i]]$regression[[1]][2]
#         }
#      }
#   #   detach(pos=2)
#   #}

# Direct write from array.

# Available statistics.
#   "raw_counts","pass_counts","mean","variance","Q50","Q75","Q90","Q95","Q99","maxium","Q50_mean-pass","Q50_Q50-pass"
   stat_idx <- 1
   band_idx <- 1

   lab_stat <- list_IMOS_stats[[1]]$lab_stat[stat_idx]
# Raw counts.
   mat_plot_data <- list_IMOS_stats[[2]][,,stat_idx,band_idx]
# Pass counts.
#   mat_plot_data <- array_stats[,,5,1]
# Difference in mean based on raw data from KU band.
#   mat_plot_data <- 100 * ( array_stats[,,1,1,1] - array_stats[,,2,1,1] ) / array_stats[,,1,1,1]
# Normalised mean difference.
#   mat_plot_data <- ( array_stats[,,1,1,1] - array_stats[,,2,1,1] ) / array_stats[,,1,1,2]

#-------------------------------------------------------#
# Set up data frame for CPC CONUS grid.
   df_plot <- NULL
   plot_labels <- paste(lab_dataset,": ",lab_stat," (KU-band only)",sep="")
   for (k in 1:1) {
      df_plot <- rbind( df_plot,
                            cbind( expand.grid( lat=rev(lat_mid), lon=lon_mid ), plot_stat=as.vector(mat_plot_data), anal=plot_labels[k] )
                          )
   }

# Custom 3.
   colfunc <- colorRampPalette(c("white","khaki1","orangered","purple","black"))
   #colfunc <- colorRampPalette(c("red","white","blue"))
   seq.plot_cols <- colfunc(20)

# Colour ranges.
# Colour ranges for counts.
   if ( stat_idx == 1) {
      plot_breaks_lo <- 0.9*min(df_plot$plot_stat,na.rm=T) - ((0.9*min(df_plot$plot_stat,na.rm=T)) %% 100)
      plot_breaks_hi <- 100 + (1.1*max(df_plot$plot_stat,na.rm=T) - ((1.1*max(df_plot$plot_stat,na.rm=T)) %% 100))
      col_breaks <- seq(plot_breaks_lo,plot_breaks_hi,500)
      plot_lims <- c(0.95*min(df_plot$plot_stat,na.rm=T),1.05*max(df_plot$plot_stat,na.rm=T))
   } else if ( stat_idx == 3 | stat_idx == 5 ) {
      col_breaks <- seq(0,3.5,0.50)
      plot_lims <- c(0,3.5)
   } else if ( stat_idx == 6 ) {
      col_breaks <- seq(2,5,0.5)
      plot_lims <- c(1,5)
   } else if ( stat_idx == 7 ) {
      col_breaks <- seq(2,6,0.5)
      plot_lims <- c(2,6)
   } else if ( stat_idx == 8 ) {
      col_breaks <- seq(2,8,0.5)
      plot_lims <- c(2,8)
   } else if ( stat_idx == 11 | stat_idx == 12 ) {
      col_breaks <- seq(0,3.5,0.5)
      plot_lims <- c(0,3.5)
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
         #geom_point(aes(x = lon, y = lat, color = plot_stat), size = 64, shape = 15) +
         scale_fill_gradientn( colours = seq.plot_cols,
                                #name = expression("SD(" * hat(chi) * "*)"),
                                name = lab_stat,
                                limits = plot_lims,
                                breaks=col_breaks,
                                na.value = "transparent",
                                guide = guide_colorbar(title.position = "bottom",label.position = "left",ticks.colour = "black",ticks.linewidth = 8.0) ) +
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
   fig_file_name <- paste("./figures/",lab_dataset,"_",lab_stat,".png",sep="")
   png(filename = fig_file_name, width = 2400, height = 1300)
   plot(p1)
   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))


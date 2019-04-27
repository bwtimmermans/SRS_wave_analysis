# source("/home/ben/research/NOC/SRS_wave_analysis/datasets/altimeter_tracks/plot_cryosat.R")
# Credit: http://www.mazamascience.com/WorkingWithData/?p=1494

# Libraries.
   library(ggplot2)
   library(sp)
   library(rgdal)
   library(rgeos)

# Load data.
   dataProjected <- readOGR(dsn="./CRYOSAT/SUBCYCLE_KML_FILES_FOR_CYCLE_01_2010-06-14_to_2011-06-17/CYCLE_01_SUBCYCLE_01.kml")
   dataProjected@data$id <- rownames(dataProjected@data)
   watershedPoints <- fortify(dataProjected, region = "id")
   watershedDF <- merge(watershedPoints, dataProjected@data, by = "id")

# Create plot object.
   ggWatershed <- ggplot(data = watershedDF, aes(x=long, y=lat, group = group)) +
                         xlim(0,50) + ylim(0,50) + coord_equal() +
                         geom_polygon() + geom_path(color = "blue") + scale_fill_hue(l = 40) +
                         theme(legend.position = "none", title = element_blank(),axis.text = element_blank())

# Plot.
   X11()
   print(ggWatershed)


---
title: "Uncovering urban circadian pulses based on an animated cartogram: the example of Bogotá"
author: "[Hugo Thomas](https://perso.univ-rennes2.fr/hugo.thomas) and [Florent Demoraes](https://perso.univ-rennes2.fr/florent.demoraes)"
date: "06/02/2023"
output: 
  html_document: 
    toc: yes
    
---

## Preliminary work
### Defining the working environment
```{r setup, include=FALSE}
    knitr::opts_chunk$set(warning = FALSE, message = FALSE, verbose = FALSE)
```

### Loading packages

```{r }
    library(downloader) # to download files from a server
    library(tidyverse) # to simplify the syntax
    library(sf) # to use SimpleFeature objects
    library(sp) # to use Spatial*DataFrame objects
    library(spatstat) # to compute the spatial smoothing
    library(maptools) # to coerce sp object into ppp object (spatstat)
    library(cartography) # to map the results
    library(cartogramR) # to create cartograms
    library(raster)  # to use rasters
    library(terra) # to use rasters
    library(stars) # to vectorize  rasters
    library(rmapshaper) # to generalize edges (simplifying geometries)
    library(dplyr) # to use data tables
    library(openxlsx) # to read EXCEL files
    library(ggplot2) # to make attractive graphs and charts
    library(tidyr) # to use matrix
    library(rgdal) # to use files with a spatial component
    library(rgeos) # to compute centroids
    library(spdep) # to compute spatial autocorrelation
    library(geoR) # to compute the empirical semivariogram
    library(magick) # to make an animated gif
    library(av) # to make an mp4 video
    library(knitr) # to display the exported pictures
    library(filesstrings) # to move a file from a folder to another one
```

### Downloading and unzipping the dataset
```{r }
    download("https://github.com/ESO-Rennes/Animated-Cartograms/raw/main/Data_Pulsation.zip", destfile=paste0(getwd(),"/Data_Pulsation.zip"), mode="wb", overwrite = TRUE)
    unzip(zipfile = "Data_Pulsation.zip", exdir = ".")
```

### Loading the GIS layer (map background)
```{r results='hide'}
    EMU2019 <- st_read(dsn = "Data_Pulsation", layer = "EMU2019", stringsAsFactors = FALSE)
```

### Loading the table with the survey data
```{r}
    viajes <- read.xlsx("Data_Pulsation/ViajesEODH2019.xlsx") #134497 observations
```

### Defining a function to render transparent color (useful for the crossfade)

```{r}
    t_col <- function(color, percent = 50, name = NULL) {
      rgb.val <- col2rgb(color)
      t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                   max = 255,
                   alpha = (100 - percent) * 255 / 100,
                   names = name)
    }  
    t_col_palette <- function(palette, percent = 50) {
      palette_trans <- palette
      for (j in 1:length(palette_trans)){
        palette_trans[j] <- t_col(palette_trans[j], percent = trans)
      }
      return(palette_trans)
    }
```

## Creating the table with the population stock per UTAM with a 15-minute timestep

Computing the trip duration (we take into account the trips that end the following day)
    
```{r}
    viajes$duracion <- viajes$p31_hora_llegada - viajes$hora_inicio_viaje

    for(i in 1:nrow(viajes)){
        if(viajes$p31_hora_llegada[i]<viajes$hora_inicio_viaje[i]){
        viajes$duracion[i] <- viajes$duracion[i]+1
        }
    }
```
    
Computing the average trip duration in hours. The average trip duration (one-way) is 50 minutes.

```{r}
    avg = 24*mean(viajes$duracion) 
```

Removing the 0-minute trips
    
```{r}
    viajes2 <- viajes %>% filter(duracion>0)
```

Processing the trips whose departure or arrival points are outside Bogotá DC. We create a virtual UTAM which pools all the zones outside the DC.

```{r}
    # tracking the trips which leave or enter the DC 
    excl <- c("N/A", "UPR1", "UPR2", "UPR3")

    viajes2$utam_destino[is.na(viajes2$utam_destino)] <- "N/A"
    viajes2$utam_origen[is.na(viajes2$utam_origen)] <- "N/A"


    # creating a virtual UTAM "UTAM99"
    for (i in 1:nrow(viajes2)){
      if (viajes2$utam_destino[i] %in% excl){
        viajes2$utam_destino[i] <- "UTAM999"
      }
      if (viajes2$utam_origen[i] %in% excl){
        viajes2$utam_origen[i] <- "UTAM999"
      }
    }
```

People leave the beginning UTAM at the departure time and enter the destination UTAM at the arrival time.

```{r}
    # Hours and duration are expressed as a fraction of the day. We multiply them by a factor 24 to convert them in hours of the day. 
    viajes2$h_arrivee <- 24*viajes2$p31_hora_llegada
    viajes2$h_depart <- 24*viajes2$hora_inicio_viaje

    # Animation timestep: default 15 minutes
    pas_de_temps <- 0.25

    # We round the hours down using 15-minute slots. E.g. 8.27 am --> 8.15 am
    viajes2$h_arrivee_round <- pas_de_temps*(viajes2$h_arrivee-viajes2$h_arrivee%%pas_de_temps)%/%pas_de_temps
    viajes2$h_depart_round <- pas_de_temps*(viajes2$h_depart-viajes2$h_depart%%pas_de_temps)%/%pas_de_temps
```

We filtered the trips dataset to select only those which are part of a daily circular chain to produce a loop animation without any “jitter” caused by different nighttime populations on day D and day D+1. For every trip of a given individual, we compute the sum *Arrival UTAM - Departure UTAM*, then we add up these values for every individual. The people who make circular chain show a "0" value. We can remove the other ones, who show a strictly positive or negative value. *This processing removes 5.7% of the initial information.*
   
```{r}
    viajes2$utam_origen_id <- as.numeric(gsub("UTAM", "", viajes2$utam_origen))
    viajes2$utam_destino_id <- as.numeric(gsub("UTAM", "", viajes2$utam_destino))
    viajes2$utam_test <- viajes2$utam_destino_id - viajes2$utam_origen_id

    viajes_circulares_2 <- viajes2 %>%
      group_by(id_hogar, id_persona) %>%
      summarize(circular_2 = sum(utam_test))

    viajes2 <- viajes2 %>%
      left_join(viajes_circulares_2, by = c("id_hogar" = "id_hogar", "id_persona" = "id_persona"))

    viajes2 <- viajes2 %>%
      filter(circular_2 == 0)
```

We remove local trips (whithin an UTAM) which are equivalent to no trip at our geographic level of analysis, because they do not trigger a change in an UTAM population.
   
```{r}
    viajes2 <- viajes2 %>% filter(viajes2$utam_origen != viajes2$utam_destino)
```

We create a table with the population balance per UTAM.

```{r}
    # Population balance per UTAM for a given slot, weighted by the expansion factor
    sorties_utam <- viajes2 %>%
      group_by(h_depart_round,utam_origen) %>%
      summarize(sorties = sum(f_exp))

    entrees_utam <- viajes2 %>%
      group_by(h_arrivee_round,utam_destino) %>%
      summarize(entrees = sum(f_exp))


    # Creating a table with the population balance per UTAM.
    flux <- data.frame(EMU2019$UTAM)
    names(flux)[names(flux) == 'EMU2019.UTAM'] <- 'UTAM'

    flux <- flux %>%
      full_join(sorties_utam, by = c("UTAM" = "utam_origen"))

    flux <- flux %>%  
      full_join(entrees_utam, by = c("UTAM" = "utam_destino",  "h_depart_round" = "h_arrivee_round"))

    flux[is.na(flux)] <- 0

    flux$solde <- flux$entrees - flux$sorties
    names(flux)[names(flux) == 'h_depart_round'] <- 'h_round'


    # Sorting the values by hour (slot) and then by UTAM
    flux <- arrange(flux, h_round, UTAM)

    # Converting the table as a matrix
    f <- pivot_wider(flux, id_cols = UTAM, names_from = h_round, values_from = solde)
    f[is.na(f)] <- 0
```

We chose the base time slot when we may consider that everybody is at home. To this end, we compute the overall number of departures for each time slot and chose the slot with minimum departures. We plot a histogram. We assign the variable *NUMPERSTOT*, which gives the night population according to the 2018 census, to this slot.

```{r}
    #arrivals
    ggplot(flux, aes(x = h_round, y = entrees)) +
      stat_summary(geom = "bar", position = "dodge", fun = sum)

    #departures
    ggplot(flux, aes(x = h_round, y = sorties)) +
      stat_summary(geom = "bar", position = "dodge", fun = sum)
```



```{r}
    # Creating a table with the population stock per UTAM as an iterative sum on the time slots.

    stock <- f %>%
      inner_join(EMU2019, by = 'UTAM', copy = TRUE) %>% # joining the EMU_2019
      relocate('NUMPERSTOT', .after = '2.75') # Assigning NUMPERSTOT to the 2.45 am time slot

    # Reordering columns starting at 3.00 am 
    stock <- stock[,1:98] 
    stock1 <- stock[,1]
    stock2 <- stock[,2:13]
    stock3 <- stock[,14:98]

    stock <-cbind(stock1, stock3, stock2)

    # Filling the stock matrix (iterative sum)
    for (i in 3:ncol(stock)){
      stock[,i]<-stock[,i]+stock[,i-1]
    }

    # The possible negative values are replaced by zeros (if everything went right in the previous step, no UTAM is affected)
    for (i in 3:ncol(stock)){
      for (j in 1:nrow(stock)){
        if(stock[j,i]<0){
          stock[j,i] <- 0
        }
      }
    }
```


Interchanging rows and columns in the table to plot graphs of the population change by time of the day in a couple of UTAM.
    
```{r}
    stock2 <- stock
    rownames(stock2)<-stock2$UTAM # Renaming rows with UTAM codes
    stock2 <- stock2[,3:98] # Removing the first two rows
    tstock <- data.frame(t(stock2))

    tstock$hours <- rownames(tstock)
    tstock$rank <- c(1:96)

    # Plotting the graph for two UTAM as an example: EL LUCERO and CHAPINERO
    ggplot(tstock, aes(x = rank)) +                    
      geom_line(aes(y=UTAM67, col="EL LUCERO")) +  
      geom_line(aes(y=UTAM99, col="CHAPINERO")) + 
      scale_colour_manual("", breaks = c("EL LUCERO", "CHAPINERO"), values = c("red", "blue")) +
      scale_x_continuous(breaks=c(1+4*(0:23)), labels=c(3:23,0:2)) +
      scale_y_continuous(breaks=c(50000,100000), labels=c("50,000","100,000")) +
      labs(x="Time of the day", y="Population") +
      theme(legend.position = "left")
```
      

## First cartogram

```{r}
    # Joining the stock table to the data table with map background
    EMU2019_2 <- EMU2019[,c(2,55)] %>%
      left_join(stock, by = "UTAM")

    # Creating a list with well-ordered hours
    heures <- c(0.25*(12:95),0.25*(0:11))
    # Creating the labels to be inserted in the title of the maps from the list above (English time notation)
    heures2 <- heures
    for(i in 41:84){
      heures2[i] <- heures2[i]-12
    }
    for(i in 85:88){
      heures2[i] <- heures2[i]+12
    }
    
    heures_lab <- as.character(heures2)
    
    for(i in 1:length(heures_lab)){
      heures_lab[i] <- paste(as.character(heures2[i]%/%1), 
                         ".",as.character(60*heures2[i]%%1), sep = "")
      if(i%%4 == 1){heures_lab[i] <- str_replace(heures_lab[i], fixed(".0"), fixed(".00"))}
    }
    
    for(i in 1:length(heures_lab)){
      if(i %in% 37:84){
        heures_lab[i] <- paste(heures_lab[i], "pm", sep = " ")
      }
      else{
        heures_lab[i] <- paste(heures_lab[i], "am", sep = " ")
      }
    }
    
    # Renaming the columns with the hours
    for (i in 1:ncol(EMU2019_2)){ 
      names(EMU2019_2)[names(EMU2019_2) == heures[i]] <- paste("stock",heures[i], sep = "_")
    }
```

First trial to export cartograms, individual hours (e.g. 3.00 am). We define an absolute error of 10 hectares, or 100,000 m². Computing limited to 100 iterations.
    
```{r }
    EMU2019_carto <- cartogramR(EMU2019_2, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 100, absrel = FALSE, abserror = 100000, L = 256))

    titre <- paste("Population per UTAM in Bogotá", heures_lab[1], sep = " - ")
    png(paste(titre,".png"), width = 960, height = 960)
    plot(EMU2019_carto, color = "black")
    title(main = titre, cex.main = 2)
    dev.off()
    
    include_graphics(paste(titre,".png"))
        
    #conversion en fichier sf
    EMU2019_carto_sf <- as.sf(EMU2019_carto)
```

## Rendering the second variable: the population multiplication factor with respect to night population, per 15-minute slot. 

We compute multiplication factors from the stock table. We fill the table one column after another dividing the stock for each slot by the night stock *NUMPERSTOT*.

```{r}
    variation <- stock

    for (i in 3:ncol(variation)){ 
      variation[,i]<-variation[,i]/variation[,2]
      # Renaming the multiplication factor columns to distinguish them from the stock columns when we join both tables 
      names(variation)[names(variation) == heures[i]] <- paste("var",heures[i], sep = "_") 
    }
    names(variation)[names(variation) == heures[1]] <- paste("var",heures[1], sep = "_")
    names(variation)[names(variation) == heures[2]] <- paste("var",heures[2], sep = "_")

    # Joining the multiplication factors table to the data table with map background
    EMU2019_3 <- EMU2019_2 %>%
      left_join(variation, by = "UTAM")
```

We define breaks for 6 classes
    
```{r}
    bks <- c(0,0.7,0.9,1.1,2,5,50)
```

We create a color palette (diverging color range)
    
```{r}
    cols <- c("#303b98","#00beed","#fafaaa","#f2bf26","#e96831","#bb2b30")
    # Transparency factor
    trans <- 50
    # Applying the transparency factor
    for (j in 1:6){
      cols[j] <- t_col(cols[j], percent = trans)
    }
```

### Displaying the choropleth map at 12.00 pm as an example

```{r}
      titre <- paste("Multiplication factor with respect to night population", heures_lab[37], sep = " - ")
      png(paste(titre, ".png"), width = 960, height = 960)
      plot(st_geometry(EMU2019_3))
      choroLayer(x = EMU2019_3, var = paste("var", heures[37], sep = "_"), breaks = bks, col = cols, legend.title.txt = "Multiplication factor", legend.values.rnd = 2, legend.values.cex = 1.5, legend.title.cex = 2)
      title(main = paste("Multiplication factor with respect to night population", heures_lab[37], sep = "\n"), cex.main = 2)
      dev.off()
      
      include_graphics(paste(titre,".png"))
```

## Spatial smoothing
### Preliminary processing (computing the spatial autocorrelation and variogram)

```{r}
    # Coercing the UTAM layer into a spatial object (required by geoR)

    UTAM <- as(EMU2019_3, "Spatial")

    # Computing UTAM centroids (necessary to compute the nearest neighbors)
    
    UTAMCentroids <- gCentroid(UTAM,byid=TRUE)

    # Recovering UTAM initial data on centroids

    UTAMCentroids <- SpatialPointsDataFrame(UTAMCentroids, UTAM@data)

    # Computing the nearest neighbors

    listPPV <- knearneigh(UTAMCentroids@coords, k = 1) # Finding each UTAM's nearest neighbor
    PPV <- knn2nb(listPPV, row.names = UTAM$UTAM) # Coercing knn objects into nb objects
    distPPV <- nbdists(PPV, UTAMCentroids@coords) # Computing the distance between nearest neighbors
    print(as.data.frame(t(as.matrix(summary(unlist(distPPV))))))
    hist(unlist(distPPV), breaks = 20,
         main = "Distance to the closest neighbor",
         col = "black", border = "white", xlab = "Distance", ylab = "Frequence")
```
         
Most of the UTAM have at least one neighbor in a 1500m bandwidth.

```{r}
    # Coercing UTAMS into nb objects
    nbUTAM <- poly2nb(pl = UTAM,
                     row.names = UTAM$UTAM,
                     snap = 50,
                     queen = TRUE)

    # Identifying the UTAM without topological neighbor 
    summary(nbUTAM)
```
    
    

```{r}
    # Creating a list of the 16 UTAM without neighbors found at the previous step
    UTAM_isolees <- c("UTAM89",
              "UTAM563",
              "UTAM540",
              "UTAM580",
              "UTAM640",
              "UTAM650",
              "UTAM660",
              "UTAM670",
              "UTAM680",
              "UTAM690",
              "UTAM700",
              "UTAM600",
              "UTAM630",
              "UTAM620",
              "UTAM610",
              "UTAM590")

    # Deleting "UTAM_isolees" otherwise we cannot compute the Moran's I
    UTAM <- UTAM[which(!UTAM$UTAM %in% UTAM_isolees),]
    EMU2019_3 <- EMU2019_3[which(!EMU2019_3$UTAM %in% UTAM_isolees),]
    nbUTAM <- poly2nb(pl = UTAM,
                     row.names = UTAM$UTAM,
                     snap = 50,
                     queen = TRUE)

    # Computing the Moran's I for each time slot
    moran <- c(0*(1:96))
    for(i in 1:96){
      m <- moran.test(UTAM@data[,100+i], listw = nb2listw(nbUTAM)) # In the UTAM table, the columns with the multiplication factors are numbered 101 to 196
      moran[i] <- as.numeric(m$estimate[1])
    }
    # Displaying the Moran's I
    df <- data.frame(cbind(moran,heures))
    ggplot(df) +  
      theme_classic() +
      geom_point(aes(x = heures, y=moran), shape = 21, fill = "white", color = "black", size = 2, stroke = 0.5)  +
      labs(x="Time of the day", y="Moran") +
      theme(legend.position = "left", 
            panel.border = element_rect(colour = "black", fill=NA, size=0.5))
```

### Variogram

```{r}
    # Computing the UTAM pseudo-centroids without the isolated UTAM (necessary to comupute the semi-variogram)
    UTAMCentroids <- gPointOnSurface(UTAM,byid=TRUE)

    # Recovering UTAM initial data on centroids
    UTAMCentroids <- SpatialPointsDataFrame(UTAMCentroids, UTAM@data)

    # Coercing the SpatialPointsDataFrame into geodata object
    UTAMCentroids.geodata <- as.geodata(UTAMCentroids, data.col = "var_12")

    # Computing the empirical semivariogram
    vario.ex<- variog(UTAMCentroids.geodata, bin.cloud=TRUE, option = "bin")
    plot(vario.ex, main = "Semivariogram of the multiplication factor at 12.00 pm in function of the distance", cex.main = 1)
    lines(vario.ex, type ="l", lty = 2, col="red")
```

### Smoothed map at 12.00 pm

```{r}
    # Defining the edge of Bogotá as the extent of the map
    Emprise <- as.owin(gBuffer(UTAM, width=0))

    # Creating a ppp object (spatstat) and including the extent and values to be smoothed (i.e. the multiplication factor at 12.00 pm) 
    UTAM.ppp <- ppp(UTAMCentroids@coords[,1], UTAMCentroids@coords[,2], window = Emprise, marks = UTAM$var_12)

    # Computing the smoothed surface (smoothing bandwidth: 1 km and picture spatial resolution: 1 ha) --> takes some time
    cartelissee <- Smooth(UTAM.ppp, kernel = "gaussian", sigma = 1000, weights = UTAM.ppp$marks, eps=c(100,100))

    # Coercing the smoothed surface into a raster
    cartelissee.raster <- raster(cartelissee)
    crs(cartelissee.raster) <- st_crs(UTAM)$srid # to specify a CRS to the raster object

    # Configuring the window margins to maximize the extent of the map
    par(mar = c(0, 0, 0, 0))

    # Displaying all the UTAM in the background and centering the map on Bogotá.
    plot(st_geometry(EMU2019_3), col = NA, border = "black", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))

    # Displaying the smoothed map
    plot(cartelissee.raster, breaks = bks, col=cols, add = T, legend=F)

    # Displaying the UTAM boundaries
    plot(st_geometry(EMU2019_3), border = "black", lwd = 0.05, lty=3, add = T)

    legendChoro(
      pos = "bottomleft",
      title.txt = "Multiplication factor",
      breaks = bks, 
      nodata = FALSE,
      values.rnd = 2,
      col = cols,
      cex = 1.2,
      values.cex = 1,
      title.cex = 1
    )

    title("Population multiplication factor at 12.00 pm \n Smoothing radius 1000m", cex.main = 1, line = -2)
```

We performed various trials before selecting a smoothing bandwidth of 1,000 meters. We also tested various kernels before choosing a Gaussian function. This function offers a good compromise between the weight given to close neighbors and that given to remote ones. 

## Isolating the UTAM "El Rincon de Suba" to display it as a legend on the side of the map and show the scale of deformation

```{r}
  EMU2019_carto <- cartogramR(EMU2019_2, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 100, absrel = FALSE, abserror = 100000, L = 256))
  EMU2019_carto_sf <- as.sf(EMU2019_carto)
  RINCON <- EMU2019_carto_sf[which(EMU2019_carto_sf$UTAM == "UTAM28"),]
  st_geometry(RINCON)[[1]] = st_geometry(RINCON)[[1]] + st_point(c(-25000,-10000))
  RINCON$legend <- paste(RINCON$NUMPERSTOT, "people", sep = "\n")
```


## Initial distortion of the map
### Cartogram only
To do so, we carry out an interpolation between the undistorted map of Bogotá and the map distorted according the the night population at 3.00 am.

```{r results='hide'}
    
    # Distorting the initial map - cartogramm only
    dir.create("Animated_Cartogram")

    # Computing areas
    EMU2019_area <- EMU2019[which(!EMU2019$UTAM %in% UTAM_isolees),]

    EMU2019_area$Area <- as.numeric(st_area(st_geometry(EMU2019_area)))

    # Computing the total population
    Population <- sum(EMU2019_area$NUMPERSTOT)

    # Computing the total area
    Superficie <- sum(EMU2019_area$Area)

    # Number of intermediate pictures during the distortion
    imgs <- 8

    # Distributing the total population proportionally to the geographic area of each UTAM - Picture 0
    EMU2019_area$NUMPERSTOT_0 <- EMU2019_area$Area*Population/Superficie

    # The field coding the fictitious area for each picture
    EMU2019_area$NUMPERSTOT_i <- EMU2019_area$NUMPERSTOT

    # Creating intermediate pictures with a linear interpolation (produces a smoother change than the annual average growth rate would)
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    par(mar = c(1, 1, 1, 1))

    # Displaying the undistorted map
    png(paste0(getwd(),"/Animated_Cartogram/", "deformation 0.png"), width = 960, height = 960)
    plot(st_geometry(EMU2019_area), border = "black", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
    dev.off()

    for (i in 1:imgs){
      EMU2019_area$NUMPERSTOT_i <- EMU2019_area$NUMPERSTOT*i/imgs + EMU2019_area$NUMPERSTOT_0*(imgs-i)/imgs
      carto <- cartogramR(EMU2019_area, count = "NUMPERSTOT_i", options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      png(paste0(getwd(),"/Animated_Cartogram/", paste0("deformation ", i,".png")), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = NA, xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      carto_sf <- as.sf(carto)
      plot(st_geometry(carto_sf), border = "black", lwd = 1, add = TRUE)
      dev.off()
    }
        
```

### Choropleth cartograms

```{r results='hide'} 
  #Déformation du fond de carte initial pour le cartogramme seconde variable. 
    dir.create("Choropleth_animated_cartogram")

    # Computing areas
    EMU2019_area <- EMU2019[which(!EMU2019$UTAM %in% UTAM_isolees),]

    EMU2019_area$Area <- as.numeric(st_area(st_geometry(EMU2019_area)))

    # Computing the total population
    Population <- sum(EMU2019_area$NUMPERSTOT)

    # Computing the total area
    Superficie <- sum(EMU2019_area$Area)

    # Number of intermediate pictures during the distortion
    imgs <- 8

    # Distributing the total population proportionally to the geographic area of each UTAM - Picture 0
    EMU2019_area$NUMPERSTOT_0 <- EMU2019_area$Area*Population/Superficie

    # The field coding the fictitious area for each picture
    EMU2019_area$NUMPERSTOT_i <- EMU2019_area$NUMPERSTOT

    # Creating intermediate pictures with a linear interpolation (produces a smoother change than the annual average growth rate would do)
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    par(mar = c(1, 1, 1, 1))

    # Displaying the undistorted map
    png(paste0(getwd(),"/Choropleth_animated_cartogram/", "deformation 0.png"), width = 960, height = 960)
    plot(st_geometry(EMU2019_area), col = NA, border = "black", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))

    dev.off()

    for (i in 1:imgs){
      EMU2019_area$NUMPERSTOT_i <- EMU2019_area$NUMPERSTOT*i/imgs + EMU2019_area$NUMPERSTOT_0*(imgs-i)/imgs
      carto <- cartogramR(EMU2019_area, count = "NUMPERSTOT_i", options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      png(paste0(getwd(),"/Choropleth_animated_cartogram/", paste0("deformation ", i,".png")), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      carto_sf <- as.sf(carto)
      plot(st_geometry(carto_sf), border = "black", lwd = 1, add = TRUE)
      dev.off()
    }
    
    # Duplicate the first initial pictures
    dir.create("Smoothed_binned_colored_animated_cartogram")
    initial_files <- list.files(paste0(getwd(),"/Choropleth_animated_cartogram/"))
    file.copy(from = paste0(getwd(),"/Choropleth_animated_cartogram/", initial_files),
          to = paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram/", initial_files))

    dir.create("Smoothed_continuous_colored_animated_cartogram")
    file.copy(from = paste0(getwd(),"/Choropleth_animated_cartogram/", initial_files),
          to = paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram/", initial_files))
    
```

## Crossfade arrival of the legend
### Cartogram only

```{r}
    setwd(paste0(getwd(),"/Animated_Cartogram"))
    carto <- cartogramR(EMU2019_3, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
  
    # Level of transparency (%)
  
  for (i in c(90, 80, 70, 60, 50)) {
    trans <- i
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    par(mar = c(1, 1, 1, 1))
    
    # Displaying all the UTAM in the background
    titre <- paste("transition", 100 - trans, "percent", sep = " - ")
    png(paste(titre,".png"), width = 960, height = 960)
    plot(st_geometry(EMU2019_3), col = NA, border = "white", bg = "white", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
    
    par(col.main = paste("grey", 1.8*i-90))
    title(main = "Population per UTAM in Bogotá - 3.00 am", cex.main = 2, line = -2) 
    
    # Coercing the cartogram into a sf object
    carto_sf <- as.sf(carto)
    plot(st_geometry(carto_sf), border = "black", lwd = 1, add = TRUE)
    
    # Displaying the UTAM EL RINCON in a corner
    plot(st_geometry(RINCON), col = NA, border = paste("grey", 1.8*i-90), add = TRUE)
    labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0, col = paste("grey", 1.8*i-90))
    
    # Displaying the undistorted map
    mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
    split.screen(figs = mat, erase=FALSE)
    plot(st_geometry(EMU2019_3), border = paste("grey", 0.6*i+30), add = TRUE)
    mtext("Initial undistorted map", side = 1, line = 3, col = paste("grey", 1.8*i-90), outer = FALSE)
    close.screen(all = TRUE)
    dev.off()
  }
```

## Crossfade arrival of the legend and the second variable
### Choropleth animated cartogram

```{r}
    setwd(paste0(getwd(),"/Choropleth_animated_cartogram"))
    
    Correspondance_palette <- as.data.frame(cbind(c(1:6), cols))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"
    
    carto <- cartogramR(EMU2019_3, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
    
    EMU2019_carto_sf <- as.sf(carto)
    
    for (i in c(90,80,70,60,50)) { 
      trans <- i
      # Displaying the background
      par(mar = c(1,1,1,1))
      titre <- paste("transition", 100 - trans, "percent", sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      choroLayer(x = EMU2019_carto_sf, var = paste("var", heures[1], sep = "_"), breaks = bks, col = cols, legend.title.txt = "Multiplication factor", legend.values.rnd = 2, legend.values.cex = 1, legend.title.cex = 1.5, add = TRUE)
      
      title(main = "Population per UTAM in Bogotá - 3.00 am", cex.main = 2, line = -2) 
      
      # Displaying the UTAM EL RINCON in a corner
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)
      
      # Displaying the undistorted map
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(n = 1)
      
      background <- t_col("grey90",200-2*i)
      mat <- matrix(c(0,1.0,0,1), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), col = NA, border = background, bg = background, xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      plot(st_geometry(EMU2019_carto_sf), border = "black", lwd = 1, add = TRUE)
      close.screen(all = TRUE)
      dev.off()
    }
```


### Smoothed binned-colored animated cartogram

```{r}
    setwd(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram"))
    
    Correspondance_palette <- as.data.frame(cbind(c(1:6), cols))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"
    
    carto <- cartogramR(EMU2019_3, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))

    # Creating a ppp object (spatstat) and including the extent and values to be smoothed (i.e. the multiplication factor at 12.00 pm) 
    UTAM.ppp <- ppp(UTAMCentroids@coords[,1], UTAMCentroids@coords[,2], window = Emprise, marks = UTAM@data[,196])
    
    # Computing the smoothed surface (smoothing bandwidth: 1 km and picture spatial resolution: 1 ha) --> takes some time
    cartelissee <- Smooth(UTAM.ppp, kernel = "gaussian", sigma = 1000, weights = UTAM.ppp$marks, eps=c(100,100))
    
    # Coercing the smoothed surface into a raster
    cartelissee.raster <- raster(cartelissee)
    crs(cartelissee.raster) <- st_crs(UTAM)$srid # pour spécifier un SCR à l'objet raster
    
    # Reclassifying the smoothed surface
    cartelissee.reclass <- cut(cartelissee.raster, breaks = bks)
    
    # Vectorizing the reclassed surface
    cartelissee.vecteur <- sf::st_as_sf(stars::st_as_stars(cartelissee.reclass), as_points = FALSE, merge = TRUE) # requires the sf, sp, raster and stars packages
    
    # Sorting the layer column information
    cartelissee.vecteur <- cartelissee.vecteur[order(cartelissee.vecteur$layer, decreasing = FALSE), ]
    
    # Creating a list of classes which are present at the time slot
    liste_classes_presentes <- as.character(cartelissee.vecteur$layer)
    
    # Slightly cropping cartelissee.vecteur (-50 m) so that the extent of this object fits into the extent of the cartogram (only do this if the next instruction does not work)
    cartelissee.vecteur <- st_crop(cartelissee.vecteur, ext(UTAM)-50)
    
    # Distorting the smoothed surface to fit closely to the cartogram
    cartelissee.vecteur.drapee <- geom_cartogramR(cartelissee.vecteur, carto)
    
    # Generalizing the edges of cartelissee.vecteur.drapee (we keep 5% of the original vertices)
    cartelissee.vecteur.drapee <- ms_simplify(cartelissee.vecteur.drapee, keep = 0.05, keep_shapes = TRUE)
    
    # Level of transparency (%)
    
    for (i in c(90, 80, 70, 60, 50)) {
      trans <- i
      cols_2 <- cols
      for (j in 1:10){
        cols_2[j] <- t_col(cols[j], percent = trans)
      }
      
      Emprise <- as.owin(gBuffer(UTAM, width=0))
      par(mar = c(1, 1, 1, 1))
      
      # Creating a correlation table between each class and its color code
      Correspondance_palette <- as.data.frame(cbind(c(1:6), cols_2))
      names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"
      
      # Creating a color palette with only the classes which are present
      couleurs_classes_presentes <- Correspondance_palette[which(Correspondance_palette$Id_Classe %in% liste_classes_presentes),]
      
      # Displaying all the UTAM in the background
      titre <- paste("transition", 100 - trans, "percent", sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      
      # Displaying the distorted smoothed surface
      typoLayer(
        x = cartelissee.vecteur.drapee,
        var="layer",
        col = as.character(couleurs_classes_presentes$cols),
        lwd = 0.1,
        border = as.character(couleurs_classes_presentes$cols),
        legend.pos = "n",
        add=T)
      
      legendChoro(
        pos = "bottomleft",
        title.txt = "Multiplication factor",
        breaks = c("","0.72","0.9","1.1","2","5",""),
        nodata = FALSE,
        values.rnd = 2,
        col = cols,
        cex = 1.2,
        values.cex = 1,
        title.cex = 1.5
      )
      
      par(col.main = paste("grey", 1.8*i-90))
      title(main = "Population per UTAM in Bogotá - 3.00 am", cex.main = 2, line = -2) 
      
      # Coercing the cartogram into a sf object
      carto_sf <- as.sf(carto)
      
      # Displaying the cartogram
      plot(st_geometry(carto_sf), border = paste("grey", 200-2*i), lwd = 1, add = TRUE)
      
      # Displaying the edges of the localities 
      Localidades_2 <- EMU2019_3 %>% group_by(LocMuni) %>% summarize()
      Localidades_carto_2 <- geom_cartogramR(Localidades_2, carto)
      plot(st_geometry(Localidades_carto_2), col = NA, border = "black", lwd = 3.25 - i/40, add = TRUE)
      
      
      # Displaying the UTAM EL RINCON as the legend
      plot(st_geometry(RINCON), col = NA, border = paste("grey", 1.8*i-90), add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0, col = paste("grey", 1.8*i-90))
      
      # Displaying the undistorted map in the bottom-right corner
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = paste("grey", 0.6*i+30), add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, col = paste("grey", 1.8*i-90), outer = FALSE)
      close.screen(all = TRUE)
      dev.off()
    }
```


### Smoothed continuous-colored animated cartogram

```{r}
    setwd(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram"))
    # Configuring the discretization (defining breaks and class bounds)
    bks_pseudo_continuous <- c(seq(from=0, to=3, by = 0.025), seq(from=3.5, to=50, by = 0.5))
    
    n0<-length(seq(from=0, to=0.6, by = 0.025))
    n1<-length(seq(from=0.625, to=0.8, by = 0.025))
    n2<-length(seq(from=0.825, to=0.975, by = 0.025))
    n3<-length(seq(from=1, to=1.5, by = 0.025))
    n4<-length(seq(from=1.525, to=3, by = 0.025))
    n5<-length(seq(from=3.5, to=10, by = 0.5))
    n6<-length(seq(from=10.5, to=50, by = 0.5))
    
    
    # Color palette (diverging color range) and transparency factor
    cols_continuous <- c("#1d235c", "#303b98","#00beed","#fafaaa","#f2bf26","#e96831","#bb2b30","#7f0000")
    trans <- 50
    
    pal0 <- colorRampPalette(colors = cols_continuous[1:2], interpolate = "linear")
    palette0 <- pal0(n0)
    palette0 <- t_col_palette(palette0, trans)
    
    pal1 <- colorRampPalette(colors = cols_continuous[2:3], interpolate = "linear")
    palette1 <- pal1(n1)
    palette1 <- t_col_palette(palette1, trans)
    
    pal2 <- colorRampPalette(colors = cols_continuous[3:4], interpolate = "linear")
    palette2 <- pal2(n2)
    palette2 <- t_col_palette(palette2, trans)
    
    pal3 <- colorRampPalette(colors = cols_continuous[4:5], interpolate = "linear")
    palette3 <- pal3(n3)
    palette3 <- t_col_palette(palette3, trans)
    
    pal4 <- colorRampPalette(colors = cols_continuous[5:6], interpolate = "linear")
    palette4 <- pal4(n4)
    palette4 <- t_col_palette(palette4, trans)
    
    pal5 <- colorRampPalette(colors = cols_continuous[6:7], interpolate = "linear")
    palette5 <- pal5(n5)
    palette5 <- t_col_palette(palette5, trans)
    
    pal6 <- colorRampPalette(colors = cols_continuous[7:8], interpolate = "linear")
    palette6 <- pal6(n6)
    palette6 <- t_col_palette(palette6, trans)
    
    palette <- c(palette0, palette1, palette2, palette3, palette4, palette5, palette6)
    
    # Continuous legend apraisal
    leg1 <- pal1(10)
    leg1 <- t_col_palette(leg1, trans)
    leg2 <- pal2(10)
    leg2 <- t_col_palette(leg2, trans)
    leg3 <- pal3(10)
    leg3 <- t_col_palette(leg3, trans)
    leg4 <- pal4(10)
    leg4 <- t_col_palette(leg4, trans)
    leg5 <- pal5(10)
    leg5 <- t_col_palette(leg5, trans)
    leg <- cbind(leg1, leg2, leg3, leg4, leg5)
    
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    
    Correspondance_palette <- as.data.frame(cbind(c(1:(n0+n1+n2+n3+n4+n5+n6)), palette))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V2"] <- "palette"
    
    carto <- cartogramR(EMU2019_3, count = paste("stock",heures[1], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
    
    # Creating a ppp object (spatstat) and including the extent and values to be smoothed (i.e. the multiplication factor at 12.00 pm) 
    UTAM.ppp <- ppp(UTAMCentroids@coords[,1], UTAMCentroids@coords[,2], window = Emprise, marks = UTAM@data[,196])
    
    # Computing the smoothed surface (smoothing bandwidth: 1 km and picture spatial resolution: 1 ha) --> takes some time
    cartelissee <- Smooth(UTAM.ppp, kernel = "gaussian", sigma = 1000, weights = UTAM.ppp$marks, eps=c(100,100))
    
    # Coercing the smoothed surface into a raster
    cartelissee.raster <- raster(cartelissee)
    crs(cartelissee.raster) <- st_crs(UTAM)$srid # pour spécifier un SCR à l'objet raster
    
    # Cutting the raster following the edges of the UTAM
    cartelissee.raster <- mask(cartelissee.raster, UTAM)
    
    # Reclassifying the smoothed surface
    cartelissee.reclass <- cut(cartelissee.raster, breaks = bks_pseudo_continuous)
    
    # Vectorizing the reclassed surface
    cartelissee.vecteur <- sf::st_as_sf(stars::st_as_stars(cartelissee.reclass), as_points = FALSE, merge = TRUE) # requires the sf, sp, raster and stars packages
    
    # Sorting the layer column information
    cartelissee.vecteur <- cartelissee.vecteur[order(cartelissee.vecteur$layer, decreasing = FALSE), ]
    
    # Creating a list of classes which are present at the time slot
    liste_classes_presentes <- as.character(cartelissee.vecteur$layer)
    
    l0 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer <= n0]
    l1 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+1):(n0+n1)]
    l2 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+1):(n0+n1+n2)]
    l3 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+1):(n0+n1+n2+n3)]
    l4 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+n3+1):(n0+n1+n2+n3+n4)]
    l5 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+n3+n4+1):(n0+n1+n2+n3+n4+n5)]
    l6 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer > (n0+n1+n2+n3+n4+n5)]
    
    # Creating a color palette with only the classes which are present
    couleurs_classes_presentes <- Correspondance_palette[which(Correspondance_palette$Id_Classe %in% liste_classes_presentes),]
    
    # Slightly cropping cartelissee.vecteur (-50 m) so that the extent of this object fit into the extent of the cartogram (only do this if the next instruction does not work)
    cartelissee.vecteur <- st_crop(cartelissee.vecteur, ext(UTAM)-50)
    
    # Distorting the smoothed surface to fit closely to the cartogram
    cartelissee.vecteur.drapee <- geom_cartogramR(cartelissee.vecteur, carto)
    
    # Generalizing the edges of cartelissee.vecteur.drapee (we keep 10% of the original vertices)
    cartelissee.vecteur.drapee <- ms_simplify(cartelissee.vecteur.drapee, keep = 0.1, keep_shapes = TRUE)
    
    # Level of transparency (%)
    
    for (i in c(50,60,70,80,90)) { 
      trans <- i
      # Displaying the background
      par(mar = c(1,1,1,1))
      titre <- paste("transition", 100 - trans, "pourcent", sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      
      # Displaying the distorted smoothed surface (values up to 0.6)
      if(length(l0)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l0, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l0]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 0.6 and 0.8)
      if(length(l1)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l1, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l1]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 0.8 and 1)
      if(length(l2)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l2, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l2]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 1 and 1.5)
      if(length(l3)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l3, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l3]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 1.5 and 3)
      if(length(l4)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l4, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l4]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 3 and 10)
      if(length(l5)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l5, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l5]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values greater than 10)
      if(length(l6)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l6, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l6]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the legend levels, without colors
      legendChoro(
        pos = "bottomleft",
        title.txt = "Multiplication factor",
        breaks = c("","0.72","0.9","1.1","2","5",""),
        nodata = FALSE,
        values.rnd = 2,
        col = c("grey90","grey90","grey90","grey90","grey90","grey90"),
        cex = 1.2,
        values.cex = 1,
        title.cex = 1.5,
        border = "grey90"
      )
      
      title(main = "Population per UTAM in Bogotá - 3.00 am", cex.main = 2, line = -2) 
      
      # Coercing the cartogram into a sf object
      carto_sf <- as.sf(carto)
      
      # Displaying the edges of the localities 
      Localidades_2 <- EMU2019_3 %>% group_by(LocMuni) %>% summarize()
      Localidades_carto_2 <- geom_cartogramR(Localidades_2, carto)
      
      # Displaying the UTAM EL RINCON in a corner
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)
      
      # Displaying the continuous legend
      mat <- matrix(c(0.013,0.14,0.013,0.297), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      
      image(1, 1:length(leg), t(as.matrix(1:length(leg))), col = leg, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n",bty = "n")
      close.screen(n = 1)
      
      # Displaying the undistorted map
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(n = 1)
      
      background <- t_col("grey90",200-2*i)
      mat <- matrix(c(0,1.0,0,1), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), col = NA, border = background, bg = background, xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      plot(st_geometry(carto_sf), border = paste("grey", 200-2*i), lwd = 1, add = TRUE)
      plot(st_geometry(Localidades_carto_2), col = NA, border = "black", lwd = 3.25 - i/40, add = TRUE)
      close.screen(all = TRUE)
      dev.off()
    }
```


## Creating the picture sequence with a 15-minute timestep
### Cartogram only

```{r}

    setwd(paste0(getwd(),"/Animated_Cartogram"))
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    
    for (i in 1:96){
      par(mar = c(1, 1, 1, 1))
      carto <- cartogramR(EMU2019_3, count = paste("stock",heures[i], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      titre <- paste("Population per UTAM in Bogotá", heures_lab[i], sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      
      # Displaying the background
      plot(st_geometry(EMU2019_3), col = NA, border = "white", bg = "white", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      
      title(main = titre, cex.main = 2, line = -2) 
      
      # Coercing the cartogram into a sf object
      carto_sf <- as.sf(carto)
      plot(st_geometry(carto_sf), border = "black", lwd = 1, add = TRUE)
      
      # Displaying the UTAM EL RINCON as the legend
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)
      
      # Displaying the undistorted map in the bottom-right corner
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(all = TRUE)
      dev.off()
    }
```

### Choropleth animated cartogram

```{r}
    setwd(paste0(getwd(),"/Choropleth_animated_cartogram"))

    Correspondance_palette <- as.data.frame(cbind(c(1:6), cols))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"

    for (i in 1:96){
      carto <- cartogramR(EMU2019_3, count = paste("stock",heures[i], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      EMU2019_carto_sf <- as.sf(carto)
      par(mar = c(1, 1, 1, 1))
      titre <- paste("Population per UTAM in Bogotá", heures_lab[i], sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      choroLayer(x = EMU2019_carto_sf, var = paste("var", heures[i], sep = "_"), breaks = bks, col = cols, legend.title.txt = "Multiplication factor", legend.values.rnd = 2, legend.values.cex = 1, legend.title.cex = 1.5, add = TRUE)
      title(main = titre, cex.main = 2, line = -2) 
      # Displaying the UTAM EL RINCON in a corner
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)
      # Displaying the undistorted map
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(all = TRUE)
      dev.off()
    }
```


### Smoothed binned-colored animated cartogram

```{r}
    setwd(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram"))
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    Correspondance_palette <- as.data.frame(cbind(c(1:6), cols))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"

    for (i in 1:96){
      par(mar = c(1, 1, 1, 1))
      carto <- cartogramR(EMU2019_3, count = paste("stock",heures[i], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      titre <- paste("Population per UTAM in Bogotá", heures_lab[i], sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)

      # Creating a ppp object (spatstat) and including the extent and values to be smoothed (i.e. the multiplication factor at 12.00 pm) 
      UTAM.ppp <- ppp(UTAMCentroids@coords[,1], UTAMCentroids@coords[,2], window = Emprise, marks = UTAM@data[,100+i])

      # Computing the smoothed surface (smoothing bandwidth: 1 km and picture spatial resolution: 1 ha) --> takes some time
      cartelissee <- Smooth(UTAM.ppp, kernel = "gaussian", sigma = 1000, weights = UTAM.ppp$marks, eps=c(100,100))

      # Coercing the smoothed surface into a raster
      cartelissee.raster <- raster(cartelissee)
      crs(cartelissee.raster) <- st_crs(UTAM)$srid # pour spécifier un SCR à l'objet raster

      # Reclassifying the smoothed surface
      cartelissee.reclass <- cut(cartelissee.raster, breaks = bks)

      # Vectorizing the reclassed surface
      cartelissee.vecteur <- sf::st_as_sf(stars::st_as_stars(cartelissee.reclass), as_points = FALSE, merge = TRUE)

      # Sorting the layer column information
      cartelissee.vecteur <- cartelissee.vecteur[order(cartelissee.vecteur$layer, decreasing = FALSE), ]

      # Creating a list of classes which are present at the time slot
      liste_classes_presentes <- as.character(cartelissee.vecteur$layer)

      # Creating a color palette with only the classes which are present
      couleurs_classes_presentes <- Correspondance_palette[which(Correspondance_palette$Id_Classe %in% liste_classes_presentes),]
      
      # Slightly cropping cartelissee.vecteur (-50 m) so that the extent of this object fit into the extent of the cartogram (only do this if the next instruction does not work)
      cartelissee.vecteur <- st_crop(cartelissee.vecteur, ext(UTAM)-50)

      # Distorting the smoothed surface to fit closely to the cartogram
      cartelissee.vecteur.drapee <- geom_cartogramR(cartelissee.vecteur, carto)

      # Generalizing the edges of cartelissee.vecteur.drapee (we keep 5% of the original vertices)
      cartelissee.vecteur.drapee <- ms_simplify(cartelissee.vecteur.drapee, keep = 0.05, keep_shapes = TRUE)

      # Displaying the background
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))

      # Displaying the distorted smoothed surface
      typoLayer(
        x = cartelissee.vecteur.drapee,
        var="layer",
        col = as.character(couleurs_classes_presentes$cols),
        lwd = 0.1,
        border = NA,
        legend.pos = "n",
        add=T)

      legendChoro(
        pos = "bottomleft",
        title.txt = "Multiplication factor",
        breaks = c("","0.72","0.9","1.1","2","5",""),
        nodata = FALSE,
        values.rnd = 2,
        col = cols,
        cex = 1.2,
        values.cex = 1,
        title.cex = 1.5
      )

      title(main = titre, cex.main = 2, line = -2) 

      # Coercing the cartogram into a sf object
      carto_sf <- as.sf(carto)
      
      # Displaying the cartogram
      plot(st_geometry(carto_sf), border = "white", lwd = 1, add = TRUE)

      # Displaying the edges of the localities 
      Localidades_2 <- EMU2019_3 %>% group_by(LocMuni) %>% summarize()
      Localidades_carto_2 <- geom_cartogramR(Localidades_2, carto)
      plot(st_geometry(Localidades_carto_2), col = NA, lwd = 2, add = TRUE)

      # Displaying the UTAM EL RINCON as the legend
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)

      # Displaying the undistorted map in the bottom-right corner
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(all = TRUE)
      dev.off()
    }
```
    
### Smoothed continuous-colored animated cartogram

```{r}
    setwd(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram"))
    
    Emprise <- as.owin(gBuffer(UTAM, width=0))
    
    Correspondance_palette <- as.data.frame(cbind(c(1:(n0+n1+n2+n3+n4+n5+n6)), palette))
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V1"] <- "Id_Classe"
    names(Correspondance_palette)[names(Correspondance_palette) ==  "V2"] <- "palette"

    for (i in 1:96){
      par(mar = c(1, 1, 1, 1))
      carto <- cartogramR(EMU2019_3, count = paste("stock",heures[i], sep = "_"), options = list(maxit = 20, absrel = FALSE, abserror = 100000, L = 256))
      titre <- paste("Population per UTAM in Bogotá", heures_lab[i], sep = " - ")
      png(paste(titre,".png"), width = 960, height = 960)
      
      # Creating a ppp object (spatstat) and including the extent and values to be smoothed (i.e. the multiplication factor at 12.00 pm) 
      UTAM.ppp <- ppp(UTAMCentroids@coords[,1], UTAMCentroids@coords[,2], window = Emprise, marks = UTAM@data[,100+i])
      
      # Computing the smoothed surface (smoothing bandwidth: 1 km and picture spatial resolution: 1 ha) --> takes some time
      cartelissee <- Smooth(UTAM.ppp, kernel = "gaussian", sigma = 1000, weights = UTAM.ppp$marks, eps=c(100,100))
      
      # Coercing the smoothed surface into a raster
      cartelissee.raster <- raster(cartelissee)
      crs(cartelissee.raster) <- st_crs(UTAM)$srid # pour spécifier un SCR à l'objet raster
      
      # Cutting the raster following the edges of the UTAM
      cartelissee.raster <- mask(cartelissee.raster, UTAM)
      
      # Reclassifying the smoothed surface
      cartelissee.reclass <- cut(cartelissee.raster, breaks = bks_pseudo_continuous)
      
      # Vectorizing the reclassed surface
      cartelissee.vecteur <- sf::st_as_sf(stars::st_as_stars(cartelissee.reclass), as_points = FALSE, merge = TRUE) # requires the sf, sp, raster and stars packages
      
      # Sorting the layer column information
      cartelissee.vecteur <- cartelissee.vecteur[order(cartelissee.vecteur$layer, decreasing = FALSE), ]
      
      # Creating a list of classes which are present at the time slot
      liste_classes_presentes <- as.character(cartelissee.vecteur$layer)
      
      l0 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer <= n0]
      l1 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+1):(n0+n1)]
      l2 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+1):(n0+n1+n2)]
      l3 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+1):(n0+n1+n2+n3)]
      l4 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+n3+1):(n0+n1+n2+n3+n4)]
      l5 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer %in% (n0+n1+n2+n3+n4+1):(n0+n1+n2+n3+n4+n5)]
      l6 <- cartelissee.vecteur$layer[cartelissee.vecteur$layer > (n0+n1+n2+n3+n4+n5)]
      
      # Creating a color palette with only the classes which are present
      couleurs_classes_presentes <- Correspondance_palette[which(Correspondance_palette$Id_Classe %in% liste_classes_presentes),]
      
      # Slightly cropping cartelissee.vecteur (-50 m) so that the extent of this object fit into the extent of the cartogram (only do this if the next instruction does not work)
      cartelissee.vecteur <- st_crop(cartelissee.vecteur, ext(UTAM)-50)
      
      # Distorting the smoothed surface to fit closely to the cartogram
      cartelissee.vecteur.drapee <- geom_cartogramR(cartelissee.vecteur, carto)
      
      # Generalizing the edges of cartelissee.vecteur.drapee (we keep 10% of the original vertices)
      cartelissee.vecteur.drapee <- ms_simplify(cartelissee.vecteur.drapee, keep = 0.1, keep_shapes = TRUE)
      
      # Displaying the background
      plot(st_geometry(EMU2019_3), col = NA, border = "grey90", bg = "grey90", xlim = c(st_bbox(UTAM)[1], st_bbox(UTAM)[3]), ylim = c(st_bbox(UTAM)[2], st_bbox(UTAM)[4]))
      
      # Displaying the distorted smoothed surface (values up to 0.6)
      if(length(l0)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l0, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l0]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 0.6 and 0.8)
      if(length(l1)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l1, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l1]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 0.8 and 1)
      if(length(l2)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l2, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l2]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 1 and 1.5)
      if(length(l3)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l3, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l3]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 1.5 and 3)
      if(length(l4)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l4, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l4]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values between 3 and 10)
      if(length(l5)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l5, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l5]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the distorted smoothed surface (values greater than 10)
      if(length(l6)>0){
        typoLayer(
          x = cartelissee.vecteur.drapee[cartelissee.vecteur.drapee$layer %in% l6, ],
          var="layer",
          col = as.character(couleurs_classes_presentes$palette[couleurs_classes_presentes$Id_Classe %in% l6]),
          lwd = 0.1,
          border = NA,
          legend.pos = "n",
          add=T)
      }
      
      # Displaying the legend levels, without colors
      legendChoro(
        pos = "bottomleft",
        title.txt = "Multiplication factor",
        breaks = c("","0.72","0.9","1.1","2","5",""),
        nodata = FALSE,
        values.rnd = 2,
        col = c("grey90","grey90","grey90","grey90","grey90","grey90"),
        cex = 1.2,
        values.cex = 1,
        title.cex = 1.5,
        border = "grey90"
      )
      
      title(main = titre, cex.main = 2, line = -2) 
      
      # Coercing the cartogram into a sf object
      carto_sf <- as.sf(carto)
      # Displaying the cartogram
      plot(st_geometry(carto_sf), border = "white", lwd = 1, add = TRUE)
      
      # Displaying the edges of the localities 
      Localidades_2 <- EMU2019_3 %>% group_by(LocMuni) %>% summarize()
      Localidades_carto_2 <- geom_cartogramR(Localidades_2, carto)
      plot(st_geometry(Localidades_carto_2), col = NA, lwd = 2, add = TRUE)
      
      # Displaying the UTAM EL RINCON in a corner
      plot(st_geometry(RINCON), col = NA, border = "black", add = TRUE)
      labelLayer(RINCON, txt = "legend", halo = TRUE, bg = "white", r = 0.05, cex = 1.2, pos = 3, font = 2, offset = 0)
      
      # Displaying the continuous legend
      mat <- matrix(c(0.013,0.14,0.013,0.297), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      
      image(1, 1:length(leg), t(as.matrix(1:length(leg))), col = leg, xlab = "", ylab = "", xaxt = "n",
            yaxt = "n",bty = "n")
      close.screen(n = 1)
      
      # Displaying the undistorted map
      mat <- matrix(c(0.7,1.0,0.1,0.4), nrow=1)
      split.screen(figs = mat, erase=FALSE)
      plot(st_geometry(EMU2019_3), border = "grey60", add = TRUE)
      mtext("Initial undistorted map", side = 1, line = 3, outer = FALSE)
      close.screen(all = TRUE)
      dev.off()
    }
```


## Creating an animated GIF from the sequence: Cartogram only

```{r}
    dir.create("Animated_Cartogram_Export")
    
    # Loading the metadata of the pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Animated_Cartogram"), full.names = TRUE))
    
    # Sorting the metadata by creation time
    details = details[with(details, order(as.POSIXct(ctime))), ] 
    
    list <- rep("0",110)
    for (i in 1:10){
      list[i] <- paste("00",i, sep = "")
    }
    
    for (i in 10:110){
      if (i>=10 & i<100){
        list[i] <- paste("0",i, sep = "")
      } else {
        list[i] <- paste(i)
      }
    }
    
    file.rename(rownames(details), paste0(getwd(),"/Animated_Cartogram/anim", list[1:110],".png"))
    
    # Reloading the metadata of the renamed pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Animated_Cartogram"), full.names = TRUE))
    
    details <- details[15:110,]

    imgs = rownames(details)

    img_list <- lapply(imgs, image_read)

    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 20)

    # Creating a GIF
    image_write(image = img_animated, path = "Animated_Cartogram_Export/animation.gif")
    include_graphics(paste0(getwd(),"/Animated_Cartogram_Export/animation.gif"))
```


## Creating a MP4 video from the sequence: cartogram only

```{r results='hide'}
    # Creating a MP4
    rootwd <- getwd()
    setwd(paste0(getwd(),"/Animated_Cartogram"))
    av_encode_video(input = list.files(getwd()),framerate = 15, "animation.mp4")
    file.move(paste0(getwd(),"/animation.mp4"),paste0(rootwd,"/Animated_Cartogram_Export"))
```

## Creating an animated GIF from the sequence: choropleth animated cartogram

```{r}
    dir.create("Choropleth_animated_cartogram_Export")
    
    # Loading the metadata of the pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Choropleth_animated_cartogram"), full.names = TRUE))
    
    # Sorting the metadata by creation time
    details = details[with(details, order(as.POSIXct(ctime))), ] 
    
    list <- rep("0",110)
    for (i in 1:10){
      list[i] <- paste("00",i, sep = "")
    }
    
    for (i in 10:110){
      if (i>=10 & i<100){
        list[i] <- paste("0",i, sep = "")
      } else {
        list[i] <- paste(i)
      }
    }
    
    file.rename(rownames(details), paste0(getwd(),"/Choropleth_animated_cartogram/anim", list[1:110],".png"))
    
    # Reloading the metadata of the renamed pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Choropleth_animated_cartogram"), full.names = TRUE))
    
    details <- details[15:110,]

    imgs = rownames(details)

    img_list <- lapply(imgs, image_read)

    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 20)

    # Creating a GIF
    image_write(image = img_animated, path = "Choropleth_animated_cartogram_Export/animation.gif")
    include_graphics(paste0(getwd(),"/Choropleth_animated_cartogram_Export/animation.gif"))
```


## Creating a MP4 video from the sequence: choropleth animated cartogram

```{r results='hide'}
    # Creating a MP4
    rootwd <- getwd()
    setwd(paste0(getwd(),"/Choropleth_animated_cartogram"))
    av_encode_video(input = list.files(getwd()),framerate = 15, "animation.mp4")
    file.move(paste0(getwd(),"/animation.mp4"),paste0(rootwd,"/Choropleth_animated_cartogram_Export"))
```

## Creating an animated GIF from the sequence: smoothed binned-colored animated cartogram

```{r}
    dir.create("Smoothed_binned_colored_animated_cartogram_Export")
    
    # Loading the metadata of the pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram"), full.names = TRUE))
    
    # Sorting the metadata by creation time
    details = details[with(details, order(as.POSIXct(ctime))), ] 
    
    list <- rep("0",110)
    for (i in 1:10){
      list[i] <- paste("00",i, sep = "")
    }
    
    for (i in 10:110){
      if (i>=10 & i<100){
        list[i] <- paste("0",i, sep = "")
      } else {
        list[i] <- paste(i)
      }
    }
    
    file.rename(rownames(details), paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram/anim", list[1:110],".png"))

    # Reloading the metadata of the renamed pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram"), full.names = TRUE))
        
    details <- details[15:110,]

    imgs = rownames(details)

    img_list <- lapply(imgs, image_read)

    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 20)
        
    # Creating a GIF
    image_write(image = img_animated, path = "Smoothed_binned_colored_animated_cartogram_Export/animation.gif")
    
    include_graphics(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram_Export/animation.gif"))
```

## Creating a MP4 video from the sequence: smoothed binned-colored animated cartogram

```{r results='hide'}
    # Creating a MP4
    rootwd <- getwd()
    setwd(paste0(getwd(),"/Smoothed_binned_colored_animated_cartogram"))
    av_encode_video(input = list.files(getwd()),framerate = 15, "animation.mp4")
    file.move(paste0(getwd(),"/animation.mp4"),paste0(rootwd,"/Smoothed_binned_colored_animated_cartogram_Export"))
```


## Creating an animated GIF from the sequence: smoothed continuous-colored animated cartogram

```{r message=FALSE, warning=FALSE}
    dir.create("Smoothed_continuous_colored_animated_cartogram_Export")
    
    # Loading the metadata of the pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram"), full.names = TRUE))
    
    # Sorting the metadata by creation time
    details = details[with(details, order(as.POSIXct(ctime))), ] 
    
    list <- rep("0",110)
    for (i in 1:10){
      list[i] <- paste("00",i, sep = "")
    }
    
    for (i in 10:110){
      if (i>=10 & i<100){
        list[i] <- paste("0",i, sep = "")
      } else {
        list[i] <- paste(i)
      }
    }
    
    file.rename(rownames(details), paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram/anim", list[1:110],".png"))
  
    # Reloading the metadata of the renamed pictures to be included in the animation
    details = file.info(list.files(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram"), full.names = TRUE))
      
    details <- details[15:110,]

    imgs = rownames(details)

    img_list <- lapply(imgs, image_read)

    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 20)

    # Creating a GIF
    image_write(image = img_animated, path = "Smoothed_continuous_colored_animated_cartogram_Export/animation.gif")
    include_graphics(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram_Export/animation.gif"))
```


## Creating a MP4 video from the sequence: smoothed continuous-colored animated cartogram

```{r results='hide'}
    # Creating a MP4
    rootwd <- getwd()
    setwd(paste0(getwd(),"/Smoothed_continuous_colored_animated_cartogram"))
    av_encode_video(input = list.files(getwd()),framerate = 15, "animation.mp4")
    file.move(paste0(getwd(),"/animation.mp4"),paste0(rootwd,"/Smoothed_continuous_colored_animated_cartogram_Export"))
```




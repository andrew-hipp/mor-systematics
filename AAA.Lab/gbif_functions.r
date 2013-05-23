###GBIF/NICHE MODELING FUNCTIONS (v1.1 5-16-2013 R code)
### To download gbif data for Cariceae sp, clean up data (remove duplicates, remove less precise georef's to create niche maps of different Taxons
### Marlene Hahn May 2013, as part of Carex project,  (based on original R coding by Marcial Escudero and Ian Pearse.)

#v1.1 additions- - now writes out files in download_gbif, and clean_gbif; pdf and jpeg maps creates with axes
#v1.1 issues- download-gbif function is running into some extraction issues that kill the function for some datasets.


## Step 1: install packages (do only once)
#libs <- c("rJava", "rgdal", "sp", "XML", "raster", "dismo", "maps", "maptools","RColorBrewer", "classInt", "mapdata", "MIPHENO")
#lapply(libs, install.packages)

## Step2: Upload species list into R(Genus	species	etc...)
#speciesfile <- read.delim(file.choose(), as.is = TRUE)    ##file must be text, tab delim for this to work; CAREX species names coming from WorldChecklist extraction. (species column needs header as "species")
#note for Carex first run, I only used species marked as accepted species, and excluded intraspecific designations. 5/14/2013- MH

## Step 3: download gbif data for different taxons; you will have to enter in a species list, and a genus name for the species list.
#EX. 							Kobresia_gbif <-download_gbif(Kobresiaspecies,"Kobresia")
#				 		    	Schoenoxiphium_gbif <-download_gbif(Schoenoxiphiumspecies,"Schoenoxiphium")
#				  				Carex_gbifdata <-download_gbif(Carexspecies,"Carex")
#				 				Uncinia_gbif <- download_gbif(Unciniaspecies,"Uncinia")
#				  				Cymophyllus_gbif <- download_gbif(Cymophyllusspecies, "Cymophyllus")

download_gbif = function(specieslist, genus) { 
## We assume that specieslist is either a dataframe or list or matrix with a "species" column, or just a vector of species epithets
  require(rJava) 
  require(rgdal)
  require(sp)
  require(XML)
  require(raster)
  require(dismo)
  if(class(specieslist) %in% c('matrix', 'list', 'data.frame')) specieslist <- specieslist$species
  gbifdata <- lapply(specieslist, function(x) {
    out <- try(gbif(genus, species=x, ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, getAlt=TRUE, ntries=5, nrecs=1000, start=1, end=NULL, feedback=3))
	if(class(out) != "try-error") save(out, file = paste(gsub(" ", "_", out$species[1]), '.Rdata', sep = ''))
	return(out)
	})
  # gbifdata <- vector('list', length(specieslist))  # defines gbifdata as list
  # for (i in 1:length(specieslist)) gbifdata[[i]] <- try(gbif(genus, species=specieslist[i], ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, getAlt=TRUE, ntries=5, nrecs=1000, start=1, end=NULL, feedback=3))
  names(gbifdata) <- specieslist
  return(gbifdata)
 }
				
##Step 4: -Flags specimens with low lat/long precision as false in precise_enough column; flags duplicate specimens (species records with same lat/long coordinates) as false in unique_record column.
			#Ex. Schoenoxiphium_cleaned <- clean_gbif(Schoenoxiphium_gbifdata)
			
 clean_gbif = function(gbifdata, clean.by.locality = FALSE) {
  for (i in 1:length(gbifdata)) {
    if(class(gbifdata[[i]]) == 'try-error') next
        else gbifdata[[i]] <- as.data.frame(gbifdata[[i]]) #Create dataframe of gbif data
        }
  xd <- list() #tempfile to use to compare data that will be flagged as unuseable
  for (i in 1:length(gbifdata)) {
    a = names(gbifdata)[i]
    print(paste("Doing", a))
    # if(a == 'diluta') browser()
    if(class(gbifdata[[i]]) %in% c('try-error', 'NULL')) next
    if(try(nrow(gbifdata[[i]])) == 0) next
# if(dim(gbifdata[[i]][1]) %in% c(NULL, 0)) next
        gbifdata[[i]]$lat <- as.numeric(gbifdata[[i]]$lat)
    gbifdata[[i]]$calc_error <- ifelse(gbifdata[[i]]$lat==as.integer(gbifdata[[i]]$lat), 100, ifelse((10*gbifdata[[i]]$lat)==as.integer(10*gbifdata[[i]]$lat), 10, ifelse((100*gbifdata[[i]]$lat)==as.integer(100*gbifdata[[i]]$lat), 1, ifelse((1000*gbifdata[[i]]$lat)==as.integer(1000*gbifdata[[i]]$lat), 0.1, ifelse((10000*gbifdata[[i]]$lat)==as.integer(10000*gbifdata[[i]]$lat), 0.01, ifelse((100000*gbifdata[[i]]$lat)==as.integer(100000*gbifdata[[i]]$lat), 0.001, 0.0001))))))
    gbifdata[[i]]$precise_enough <- ifelse(gbifdata[[i]]$calc_error < 10, TRUE, FALSE)
        gbifdata[[i]]$unique_record <- ifelse(!duplicated(gbifdata[[i]]$lat) | !duplicated(gbifdata[[i]]$lon), TRUE, FALSE) #cleans by lat and long
    # if(clean.by.locality) gbifdata[[i]]$unique_record <- gbifdata[[i]]$unique_record & ifelse(!duplicated(gbifdata[[i]]$cloc), TRUE, FALSE) -- CLEAN UP NULLS FIRST
        write.table(gbifdata[[i]], file = paste('cleaned_',gbifdata[[i]]$species[[1]],format(Sys.time(),"%Y-%m-%d"),'.txt'), sep = "|")
        xd[[i]]<-subset(gbifdata[[i]], calc_error < 10)  # can be cleaned out
    } # close i
  #browser()
  nrowlistx <- lapply(gbifdata, nrow) #this may fail now
  nrowlistxd <- lapply(xd, nrow)
  number.not.unique <- lapply(gbifdata, function(x) {
    out = 1
    if(class(x) %in% c('try-error', 'NULL')) out <- 0
    if(try(nrow(x) == 0)) out <- 0
    if(out == 1) out <- sum(!x$unique_record)
    }
    ) # end lapply function
  print("Comparison of # of original rows to # of high precision rows for LAT/LONG Coordinates; third column is number of rows not unique based on lat and long")
  print(cbind(nrowlistx,nrowlistxd,number.not.unique))
 
  return(gbifdata)
  }
 
 

##Step 5a: Create pdf maps of data (excludes specimen records flagged as low precision or as duplicate record.)
##note maps are still made if originally there was data data but flagged as false. Map without any points...
				#Ex.  Sch_log_pdf <- map_gbif(Schoenoxiphium_cleaned_dups)
map_gbif = function(gbifdata) {
  require(maps)
  require(maptools)
  require(RColorBrewer)
  require(classInt)
  require(mapdata)
  logbook <- list()
  for (i in 1:length(gbifdata)){  
    if(class(i) == "try-error") {
	  message(paste('Dataset', names(gbifdata[i]), 'is an utter failure'))
	  logbook[i] = (paste('Dataset',names(gbifdata[i]), 'is an utter failure.  Most likely download error'))
	  next
	  } # close if
	pdf(file = paste(names(gbifdata)[i],'_map_',format(Sys.time(),"%Y-%m-%d"),'.pdf',sep =''))
    map.try <- try(map("worldHires", xlim = c(min(gbifdata[[i]]$lon)-10, max(gbifdata[[i]]$lon)+10), ylim = c(min(gbifdata[[i]]$lat)-10, max(gbifdata[[i]]$lat)+10)))
    if(class(map.try) == 'try-error') {
	  message(paste('Dataset', names(gbifdata[i]), 'has some SERIOUS mapping problems. Either Null dataset or perhaps lats and Longs are be switched....Check it out.'))
	  logbook[i] =(paste('Dataset',names(gbifdata[i]), gbifdata[[i]]$species, 'has some SERIOUS mapping problems. Either Null dataset or perhaps lats and Longs are be switched....Check it out.'))
	  # add something here to delete corrupt maps and create a log file for errors
	  dev.off()
	  next
	  } # close if
	points(gbifdata[[i]]$lon[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], gbifdata[[i]]$lat[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], pch = 16, col= 2, cex = 0.5)    
	map.axes()
	title(main = gbifdata[[i]]$species[1], sub = NULL, xlab = NULL, ylab = NULL, line = NA, outer = FALSE)
    dev.off(which = dev.cur())
	logbook[i] = (paste('PDF Map generated for dataset', names(gbifdata[i])))
    } # close i
	write.table(logbook, file = paste('PDF_MAP_log',gbifdata,format(Sys.time(),"%Y-%m-%d"),'.txt'), sep = "|") ##writes out log file
	return(logbook)
  }

 ##Step 5b: Create jpeg maps of data (excludes specimen records flagged as low precision or as duplicate record.)  
 map_gbif_jpeg = function(gbifdata) {
  require(maps)
  require(maptools)
  require(RColorBrewer)
  require(classInt)
  require(mapdata)
   logbook <- list()
  for (i in 1:length(gbifdata)){  
    if(class(i) == "try-error") {
	  message(paste('Dataset', i, names(gbifdata[i]), 'is an utter failure'))
	  logbook[i] = (paste('Dataset',names(gbifdata[i]), 'is an utter failure.  Most likely download error'))
	  next
	  } # close if
	jpeg(filename = paste(names(gbifdata)[i],'_map_',format(Sys.time(),"%Y-%m-%d"),'.jpeg',sep =''), width = 480, height = 480, pointsize = 12, quality = 100, bg = "white")
    map.try <- try(map("worldHires", xlim = c(min(gbifdata[[i]]$lon)-10, max(gbifdata[[i]]$lon)+10), ylim = c(min(gbifdata[[i]]$lat)-10, max(gbifdata[[i]]$lat)+10)))
    if(class(map.try) == 'try-error') {
	  message(paste('Dataset', i, names(gbifdata[i]), 'has some SERIOUS mapping problems. Check it out.'))
	  	  logbook[i] =(paste('Dataset', names(gbifdata[i]), gbifdata[[i]]$species, 'has some SERIOUS mapping problems. Either Null dataset or perhaps lats and Longs are be switched....Check it out.'))
	  # add something here to delete corrupt maps
	  dev.off()
	  next
	  } # close if
	points(gbifdata[[i]]$lon[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], gbifdata[[i]]$lat[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], pch = 16, col= 2, cex = 0.5)    
	map.axes()
	title(main = gbifdata[[i]]$species[1], sub = NULL, xlab = NULL, ylab = NULL, line = NA, outer = FALSE)
    dev.off(which = dev.cur())
	logbook[i] = (paste('Jpeg Map generated for dataset', names(gbifdata[i])))
    } # close i
	return(logbook)
  }

 
##Step6: Download WorldClim Data (http://www.worldclim.org/download) to get bioclim variables
			#EX Sch_bioclim <- world_clim(Sch_clean)
world_clim = function(gbifdata) {
  require(rJava) 
  require(rgdal)
  require(sp)
  require(XML)
  require(raster)
  require(dismo)
  clim <-getData('worldclim', var='bio', res=5)
  bioclim <- list()
  # we also need to exclude flagged data again....  
  #worldclim and gbif have reversed long and lats??  which is why code below is [8:7]- double check with Marcial
  for (i in 1:length(gbifdata)) {
	if (gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record) {
	bioclim[[i]] <- extract(clim, gbifdata[[i]][8:7], method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, fun=NULL, na.rm=TRUE)
	}
		gbifdata[[i]]$bio1 <- ifelse(!duplicated(gbifdata[[i]]$lat) | !duplicated(gbifdata[[i]]$lon), TRUE, FALSE) 

	} #closes i
  return(bioclim)
}

##Step7: Remove OUTLIERS from Bioclim data
rm_outliers = function(bioclim){
 require(MIPHENO)
 bioclimnoout <- list()
 for(i in 1:length(bioclim)) bioclimnoout[[i]] <- rm.outliers(bioclim[[i]], fill = TRUE)
 return(bioclimoout)
}

##Step 8: Calculate MEAN and Standard Deviation (SD), then generates PCA.
mean_SD = function(bioclimoout){
 #bioclimmean <- list()
 #for (i in 1:length(bioclimnoout)) bioclimmean[[i]] <- colMeans(bioclimnoout[[i]], na.rm = TRUE, dims = 1)
 #Means <- do.call(rbind, bioclimmean)
 #row.names(Means) = specieslist
 #bioclimSD <- list()
 #for (i in 1:length(bioclimnoout)) bioclimSD[[i]] <- apply(bioclimnoout[[i]],2, sd, na.rm = TRUE)
 #SDs <- do.call(rbind, bioclimSD)
 #row.names(SDs) = specieslist
 #pca <- prcomp(Means, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL)
 #plot(pca$x[,"PC1"], pca$x[,"PC2"])
}




 
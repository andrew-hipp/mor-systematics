###GBIF/NICHE MODELING FUNCTIONS
### To download gbif data for Cariceae sp, clean up data (remove duplicates, remove less precise georef's to create niche maps of different Taxons
### Marlene Hahn May 2013, as part of Carex project.

## Step 1: install packages (do only once)
#libs <- c("rJava", "rgdal", "sp", "XML", "raster", "dismo", "maps", "maptools","RColorBrewer", "classInt", "mapdata", "MIPHENO")
#lapply(libs, install.packages)

## Step2: Download species list (Genus	species	etc...)
#speciesfile <- read.delim(file.choose(), as.is = TRUE)    ##file must be text, tab delim for this to work; CAREX species names coming from WorldChecklist extraction. (species column needs header as "species")
#note for Carex first run, I only used species marked as accepted species, and excluded intraspecific designations. 5/14/2013- MH

## Step 3: download gbif data for different taxons; you will have to enter in a species list, and a genus name for the species list.
#EX. 							Kobresia_gbifdata <-download_gbif(speciesfile = Kobresiaspecies, genus = "Kobresia")
#				 		    	Schoenoxiphium_gbifdata <-download_gbif(speciesfile = Schoenoxiphiumspecies, genus = "Schoenoxiphium")
#				  				Carex_gbifdata <-download_gbif(speciesfile = Carexspecies, genus = "Carex")
#				 				Uncinia_gbifdata <- download_gbif(speciesfile = Unciniaspecies, genus = "Uncinia")
#				  				Cymophyllus_gbifdata <- download_gbif(speciesfile = Cymophyllus_species, genus = "Cymophyllus")

download_gbif = function(specieslist, genus) { 
## We assume that specieslist is either a dataframe or list or matrix with a "species" column, or just a vector of species epithets
  require(rJava) 
  require(rgdal)
  require(sp)
  require(XML)
  require(raster)
  require(dismo)
  if(class(specieslist) %in% c('matrix', 'list', 'data.frame')) {
    specieslist <- list()
    specieslist <- specieslist$species
	}  
  gbifdata <- lapply(specieslist, function(x) {try(gbif(genus, species=x, ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, getAlt=TRUE, ntries=5, nrecs=1000, start=1, end=NULL, feedback=3))})
  # gbifdata <- vector('list', length(specieslist))  # defines gbifdata as list
  # for (i in 1:length(specieslist)) gbifdata[[i]] <- try(gbif(genus, species=specieslist[i], ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=TRUE, download=TRUE, getAlt=TRUE, ntries=5, nrecs=1000, start=1, end=NULL, feedback=3))
  names(gbifdata) <- specieslist
  return(gbifdata)
 }
		### Done for Schoenoxiphium, Cymophyllus, Uncinia, Kobresia (on May 14, 2013)
			##test success using commands like: Kobresia_gbifdata[[1]] or    Kobresia_gbifdata[[1]][7:8]  ##gives lat and long coordinates for specimens of specific Kobresia spp

			
##Step 4: -Flags specimens with low lat/long precision as LowPrecision in new column called flag_precision
			#Ex. Schoenoxiphium_cleaned <- clean_gbif(gbifdata = Schoenoxiphium_gbifdata)
			
clean_gbif = function(gbifdata, clean.by.locality = FALSE) {
  for (i in 1:length(gbifdata)) gbifdata[[i]] <- as.data.frame(gbifdata[[i]]) #Create dataframe of gbif data
  xd <- list() #tempfile to use to compare data that will be flagged as unuseable
  for (i in 1:length(gbifdata)) {
    gbifdata[[i]]$lat <- as.numeric(gbifdata[[i]]$lat)
    gbifdata[[i]]$calc_error <- ifelse(gbifdata[[i]]$lat==as.integer(gbifdata[[i]]$lat), 100, ifelse((10*gbifdata[[i]]$lat)==as.integer(10*gbifdata[[i]]$lat), 10, ifelse((100*gbifdata[[i]]$lat)==as.integer(100*gbifdata[[i]]$lat), 1, ifelse((1000*gbifdata[[i]]$lat)==as.integer(1000*gbifdata[[i]]$lat), 0.1, ifelse((10000*gbifdata[[i]]$lat)==as.integer(10000*gbifdata[[i]]$lat), 0.01, ifelse((100000*gbifdata[[i]]$lat)==as.integer(100000*gbifdata[[i]]$lat), 0.001, 0.0001))))))
    gbifdata[[i]]$precise_enough <- ifelse(gbifdata[[i]]$calc_error < 10, TRUE, FALSE)
	gbifdata[[i]]$unique_record <- ifelse(!duplicated(gbifdata[[i]]$lat) | !duplicated(gbifdata[[i]]$lon), TRUE, FALSE) #cleans by lat and long
    # if(clean.by.locality) gbifdata[[i]]$unique_record <- gbifdata[[i]]$unique_record & ifelse(!duplicated(gbifdata[[i]]$cloc), TRUE, FALSE) -- CLEAN UP NULLS FIRST
	xd[[i]]<-subset(gbifdata[[i]], calc_error < 10)  # can be cleaned out
    } # close i
  nrowlistx <- lapply(gbifdata, nrow)
  nrowlistxd <- lapply(xd, nrow)
  number.not.unique <- lapply(gbifdata, function(x) sum(!x$unique_record))
  #nrowlistx <- list()
  #nrowlistxd <- list()
  #for (i in 1:length(gbifdata)) nrowlistx[[i]] <- nrow(gbifdata[[i]])
  #for (i in 1:length(gbifdata)) nrowlistxd[[i]] <- nrow(xd[[i]])
  print("Comparison of # of original rows to # of high precision rows for LAT/LONG Coordinates; third column is number of rows not unique based on lat and long")
  print(cbind(nrowlistx,nrowlistxd,number.not.unique))
 
  return(gbifdata)
  }
 
  ##Step 6: Mapping cleaned up data  PROBLEMS #1 ## we now need to make sure we only map specimens that are NOT duplicated and are NOT low precision!!!!
				#Ex.  map_gbif(gbifdata = Schoenoxiphium_cleaned_dups)
map_gbif = function(gbifdata=Schoenoxiphium_cleaned_dups) {
  require(maps)
  require(maptools)
  require(RColorBrewer)
  require(classInt)
  require(mapdata)
  for (i in 1:length(gbifdata)){  ##Problem #2 needs to be able to skip files without any gbif data; right now it kills the function.... or at least the for loop.
    if(class(i) == "try-error") {
	  message(paste('Dataset', i, 'is an utter failure'))
	  next
	  } # close if
	pdf(file = paste(gbifdata[[i]]$species[1],'_map_',format(Sys.time(),"%Y-%m-%d"),'.pdf',sep =''))
    map("worldHires", xlim = c(min(gbifdata[[i]]$lon)-10, max(gbifdata[[i]]$lon)+10), ylim = c(min(gbifdata[[i]]$lat)-10, max(gbifdata[[i]]$lat)+10))
    points(gbifdata[[i]]$lon[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], gbifdata[[i]]$lat[gbifdata[[i]]$precise_enough & gbifdata[[i]]$unique_record], pch = 16, col= 2, cex = 0.5)    
	title(main = gbifdata[[i]]$species[1], sub = NULL, xlab = NULL, ylab = NULL, line = NA, outer = FALSE)
    dev.off(which = dev.cur())
    # gbifdata[[i]] <- gbifdata[[i]][order(gbifdata[[i]][7]),]
    # gbifdata[[i]] <- edit(gbifdata[[i]]) -- this is a throwback to doing the mapping and data cleanup on the command-line
    } # close i
  }

### PRoblem #3- do couldn't get the jpeg part to work, only the PDF

##Step6: Download WorldClim Data
 #clim<-getData(“worldclim”,var=”bio”,res=5)
 
 
 
 
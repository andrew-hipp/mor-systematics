## read data exported from SEINet and format for import into BRAHMS to update georeferencing fields.
##Hipp 2020-02-11,   MH reviewed, commented and updated 2021-03-29
##
##This will create an output file that can be imported into BRAHMS 8 by MATCH-TRANSFER using Collection Events GUID from BRAHMS as the match and transfering georeferenced fields.
##Note this code requires a file that contains the Collection Events table GUID and the matching barcode which can be exported from the Specimen Table of BRAHMS8. (or you will manually need to do a VLOOKUP)
#######This is not the SEINET GUID nor the SPecimen Table GUID, and thus is currently (20210329) not information that is in SEINET
########################################################################

library(openxlsx)
library(magrittr)

doAll <- FALSE  ## if doAll is false it will only pull georeference fields to import, if doAll is TRUE it will grab all the fields listed in fields.noImport and put them into a table.
readData <- TRUE   ## if TRUE will look for new data field to read. (I believe by the last date??) Will want this to be TRUE when downloading new SEINET file
Collect_GUID <- TRUE ##MH ADDED ....if TRUE, you've made a new export from the specimen table containing the collection events table GUID and barcode matchup. This will be used in the match-transfer to get it back into BRAHMS8. If False, it uses the previously exported file for this. Note that this is probably ok, except for records that might have been newly data entered.
editorExclude <- 'Goldberg'
editorCollapse <- list(c(from = 'Bannon|Hipp', to = 'bannon-hipp'))

if(readData) {
  dat <- read.csv(dir('../DATA.FROM.SEINET', patt = 'csv', full = T), as.is = T)
  dat <- dat[order(dat$EditId), ]
  dat$OldValue[dat$OldValue == ''] <- NA
  dat.orig <- dat
  if(length(editorCollapse) > 0) {
    for(i in seq(length(editorCollapse)))
      dat$Editor[grep(editorCollapse[[i]]['from'], dat$Editor)] <-
        editorCollapse[[i]]['to']
  }
}

###MH added. Read in new CollectionEvents to barcode matchup file exported from Specimen table of BRAHMS8. Note this isn't tested yet or implemented.
if(Collect_GUID) {
  CEGUID <- read.xlsx(dir('../BARCODE.Coll.GUID.MATCHUP', patt = 'xlsx', full = T), 1)
}


## exclude everything edited by Goldberg  ##edited fiels are county and stateProvince changes that we're not going to accept back into BRAHMS because we share those county/state formats with Living Collections
if(!identical(editorExclude, NULL)) {
  dat <- dat[grep(editorExclude, dat$Editor, invert = T), ]
}

## shove everything into the table
fields.import <- c(
                    #'catalognumber',
                    'decimallatitude',
                    'decimallongitude',
                    'geodeticdatum',
                    'georeferencedby',
                    'georeferenceremarks',
                    'georeferencesources',
                    'coordinateuncertaintyinmeters'
                  ) # close fields.import

fields.noImport <-c('country',
                    'county',
                    'cultivationstatus',
                    'dateidentified',
                    'eventdate',  #(all BHB records- most likely changed in both SEINET and BRAHMS)
                    'georeferenceprotocol',  #(1 record, by Lee, from blank to I  ??) think we can skip.
                    'habitat',
                    'identificationremarks',
                    'identifiedby',
                    'locality',
                    'localitysecurity',
                    'localitysecurityreason',
                    'municipality',
                    'occurrenceremarks',
                    'processingstatus',
                    'sciname',
                    'startdayofyear',
                    'stateprovince',
                    'stateProvince', #yes, it occurs both ways...
                    'verbatimattributes',
                    'verbatimcoordinates',
                    'verbatimeventdate', #(all BHB records- most likely changed in both SEINET and BRAHMS)
                    'year'
                  ) # close fields.noImport

fields.all <- c(fields.import, fields.noImport)

if(doAll) fields.use <- fields.all
if(!doAll) fields.use <- fields.import
dat <- dat[which(dat$FieldName %in% fields.use), ]
barcodes <- dat$CatalogNumber %>% unique %>% sort

out <- out.lastChange <-
  matrix('', length(barcodes), length(fields.use),
          dimnames = list(barcodes, fields.use)
        )
out.skipped <- numeric(0)
for(i in seq(dim(dat)[1])) {
  if(i %/% 5000 == i / 5000)
    message(paste('doing record', i, 'of', dim(dat)[1]))
  # if(!dat$FieldName[i] %in% fields.use) next # skip if not a field we care about
  dt <- dat[i, ]
  if(is.na(dt$OldValue)) {
    out.lastChange[dt$CatalogNumber, dt$FieldName] <-
      paste(dt$NewValue, dt$Editor, convertToDate(dt$Timestamp), sep = '|')
    out[dt$CatalogNumber, dt$FieldName] <- dt$NewValue
  }
  if(!is.na(dt$OldValue)) {
    lastEditor <-
      strsplit(out.lastChange[dt$CatalogNumber, dt$FieldName], '|', fixed = T)[[1]][2]
    if(identical(dt$Editor, lastEditor)) {
      out[dt$CatalogNumber, dt$FieldName] <- dt$NewValue
      out.lastChange[dt$CatalogNumber, dt$FieldName] <-
        paste(dt$NewValue, dt$Editor, convertToDate(dt$Timestamp), sep = '|')
    } else out.skipped <- c(out.skipped, dt$EditId)
  }
}

stamp <- paste(ifelse(doAll, 'allFields', 'fewFields'),
               format(Sys.time(), '%Y-%m-%d_%H.%M.%S'),
               sep = '_')


##20200329---CODE ABOVE IS NOT WORKING, as sometimes our volunteers have edited their georeference and out table has the first edit not the final edit in the final out table also not true for last edit.

### NEW MH CODE HERE Need to remove row (records) from out where there is no lat and long accepted.
out.import <- as.data.frame(out)
out.import <- out.import[!out.import$decimallatitude == "",]  ##I'm inferring that if there's no lat, there's no long but this could be incorrect assumption....

### ADD THE COLLECTION EVENTS GUID INTO THE TABLE OF IMPORTS, as barcode is not in the collection events table, but our georeference fields are. This will be needed to do the match-transfer back into BRAHMS8
out.import$MatchGUID <- CEGUID$CollectionEventId[match(row.names(out.import), CEGUID$SpecimenBarcode)]  #still need to verify this is done ok
##There are some notice some blanks in the barcode column, whcih might? match up a blank in this other file to the wrong collection events GUID
#####################################  OUTPUT FILES ########################################
write.csv(out, paste('../OUT/FulltableEdits_', stamp, '.csv', sep = ''))  ###This may not be needed. It contains the table version of the edits, and leaves the lat and long blank when we don't accept them, but leaves the row with barcode there.
writeLines(as.character(out.skipped), paste('../OUT/dataEditIdSkipped_', stamp, '.csv', sep = ''))
write.csv(out.lastChange, paste('../OUT/dataLastEdit_', stamp, '.csv', sep = ''))
write.csv(out.import, paste('../OUT/importToBrahmsTEST_', stamp, '.csv', sep = ''))  ##MH added THIS SHOULD BE THE FILE USED TO INPUT INTO BRAHMS8

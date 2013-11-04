### WCM manipulation and display scripts
### ahipp@mortonarb.org
### 2012-06-13: created
### 2013-05-13: updating to update the checklist in a more streamlined manner

do.it.cariceae <- function(datDir = choose.dir(), topNodes = c('Vignea', 'coreCarex', 'Caricoid', 'Siderostictae'), acceptedOnly = T, ...) {
## this wrapper should do everything and then archive the files with date / time stamp
## not done yet
## possible features:
##		(a) DONE: calls traverseChecklist to create a hierarchical checklist (text file, with indents)
##		(b) calls assignSectionsAndClades to assign every taxon to a section and higher clade (spreadsheet)
##		(c) DONE: highlights taxa for which no samples are available
##		(d) aggregates vouchers for synonyms up to accepted names, holding the synonym under 'Deposited as'
## 2013-07-03

  datMat <- voucherBySpecies(datDir)
  if(acceptedOnly) datMat <- datMat[datMat$Usage == 'accepted', ]
  traversedDatFile <- sapply(topNodes, function(x) traverse.checklist(datMat, topParent = x, outfileName = paste('cariceae.', x, format(Sys.time(), ".%Y-%m-%d.txt"), sep = ''), ...))
  #write.csv(datMat, 'cariceae.dataMatrix.', format(Sys.time(), "%Y-%m-%d"), '.csv', sep = '')
  write.csv(datMat, paste('cariceae.dataMatrix.', format(Sys.time(), "%Y-%m-%d"), '.csv', sep = ''))
  # zip up all the files with time stamp
  # upload to web?
  }

read.all.collector.sheets <- function(datDir, special.files = c('CHECKLIST.csv', 'GEO.csv', 'PermittedValues.csv', 'VersionNotes.csv'), rm.csv = T) {
## 13 May 2013 -- assumes all files have been exported using as CSV from standard collecting workbook
  files <- dir(datDir, full = T)
  names(files) <- dir(datDir, full = F)
  out <- files[!names(files) %in% special.files]
  out <- lapply(out, read.csv, as.is = TRUE)
  if(rm.csv) names(out) <- gsub('.csv', '', names(out), fixed = T)
  return(out)
  }
  
taxonByLab <- function(x = 'EXPORT.2013-05-13') {
## This function appends to the WCM matrix a column for each lab, indicating whether that lab has reported any material or sequences
## Arguments:
##  x = data directory to look in for all data collecting sheets and checklist
## FORMERLY CALLED makeCollSheets
## very simplistic function; not super-useful, as it ignores what quality and type of material and sequences are on-hand
  collSheets <- read.all.collector.sheets(x)
  wcm <- read.csv(dir(x, full = T, patt = 'CHECKLIST'))
  headers <- names(collSheets)
  for(i in headers) {
    message('doing', i)
	# browser()
	wcm[[i]] <- (wcm$WCM.ID %in% collSheets[[i]]$ID) & (wcm$WCM.ID != "")
	}
  return(wcm)
  }

voucherBySpecies <- function(datDir = 'EXPORT.2013-05-13', new.rank = "DNA.sample") {
## this function:
##	in a best world, we would: 
##      (a) expand the matrix with a row for each voucher, with Rank = 'VOUCHER', placing the taxon for that voucher as the parent of that voucher; note that this will have the odd side-effect of making the voucher of a species taxonomically equivalent to infraspeces of that species. Perhaps set off with a big character? e.g., -V- or ***; voucher field includes lab ownership / contact for extraction; then
##	    (b) Adds sequences as children of each voucher, with Rank = 'SEQUENCE'; sequence fields include region name and genbank number
##  INSTEAD, we add a line for each voucher with rank = 'DNA.sample', setting the term name to be the lab name, followed by the material.brief field, created uniquely for each spreadsheet.
  
  coll.sheets <- read.all.collector.sheets(datDir)
  spTaxonomy <- read.csv(paste(datDir, '/CHECKLIST.csv', sep = ''), as.is = TRUE)
  spTax.cols <- names(spTaxonomy)
  for(i in names(coll.sheets)) {
    print(paste("Doing sheet", i))
	sheet.working <- coll.sheets[[i]][gsub(" ", "", coll.sheets[[i]]$ID, fixed = TRUE) != "", ] # uses only rows that have a taxon ID matched to CHECKLIST
	sheet.as.mat <- matrix("", dim(sheet.working)[1], length(spTax.cols), dimnames = list(NULL, spTax.cols))
	sheet.as.mat[, c('Term.name', 'Parent.Term.Name', 'Rank', 'Usage')] <- cbind(as.matrix(sheet.working[, c('material.brief', 'sciname.edited')]), rep(new.rank, dim(sheet.working)[1]), rep('accepted', dim(sheet.working)[1]))
	sheet.as.mat[, 'Term.name'] <- paste('DNA VOUCHER --', i, '--', sheet.as.mat[, 'Term.name'])
	spTaxonomy <- rbind(spTaxonomy, sheet.as.mat)
	}
  return(spTaxonomy)
  }

summaryByRegion <- function(wcm = scan(file.choose(), what = character(), sep = "\n", fileEncoding = "UTF-8"), tdwg = read.delim('TDWG.regions.level2', as.is = TRUE), sampleChar = "%", toGetChar = "@", exportByRegion = F, excludeHybrids = TRUE) {
  ## get wcm by reading in a summary file using wcm <- scan(file.choose(), what = character(), sep = "\n", fileEncoding = "UTF-8")
  out <- matrix(NA, nrow = dim(tdwg)[1] + 1, ncol = 5, dimnames = list(c(tdwg$label, "GLOBAL"), c('Total species', 'Species sampled', 'Species sampled or easily obtained', 'Species not accounted for', 'Percent')))
  wcm.spp <- wcm[grep("SPECIES", wcm, fixed = TRUE)]
  if(excludeHybrids) wcm.spp <- wcm.spp[-grep("×", wcm.spp, fixed = TRUE)]
  sppNum <- length(wcm.spp)
  sampled <- seq(sppNum) %in% grep(sampleChar, wcm.spp, fixed = T)
  toSample <- seq(sppNum) %in% grep(toGetChar, wcm.spp, fixed = T)
  sampledOrObtained <- sampled | toSample
  for(i in 1:dim(tdwg)[1]) {
    temp <- seq(sppNum) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm.spp, fixed = F, perl = F)
	out[i, ] <- c(sum(temp), sum(temp & sampled), sum(temp & sampledOrObtained), sum(temp) - sum(temp & sampled), paste(formatC((sum(temp & sampled) / sum(temp)) * 100, 1, format = "f"), "%", sep = ""))
    if(exportByRegion) {
	  outfile = file(tdwg$label[i], open = "w", encoding = "UTF-8")
	  # writeLines(wcm[(seq(length(wcm)) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm, fixed = F, perl = F)) | (!seq(length(wcm)) %in% grep("SPECIES", wcm, fixed = F, perl = F))], con = outfile)
	  writeLines(sort(wcm[seq(length(wcm)) %in% grep(paste(" ", tdwg$code[i], " ", sep = ""), wcm, fixed = F, perl = F)]), con = outfile)
	  close(outfile)
	  }
	}
  out["GLOBAL", ] <- c(sppNum, sum(sampled), sum(sampledOrObtained), sppNum - sum(sampled), paste(formatC((sum(sampled) / sppNum) * 100, 1, format = "f"), "%", sep = ""))
  return(out) 
  }
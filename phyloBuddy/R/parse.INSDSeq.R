### XML GENBANK PARSER
### To Parse XML files downloaded from Geneious search of Genbank in order to get additional info such as voucher info, etc
### 2013-02: Marlene Hahn, as part of Carex project.
### 2015-03-18: Hipp and Hahn, updating to use multicores
### 2015-03-24: use only named references and pull out more needed elements
### 2015-03-25: get the desired voucher features directly from the xml

###Step 1: Save file as XML in Geneious; make sure you know how many files were downloaded
###because you will need to enter that # in function below (usually part of filename in Geneious export)
	### library(XML)
	### xml_file <- xmlTreeParse(file.choose())

 ###Issues see if I can get the program to parse out more detailed information including gene region and voucher info.
		### problem- voucher is in qualifiers in feature_table node, but not at a consistant node.
		##Still has issues parsing to the voucher level- use spliting _metadata_genbank_tables.r function to parse out this info.)



parse.INSDSeq = function(xml_file, do = NA, includeSeqs = F, cores = 1, parse.specimens = T,
                         qualsToUse = c('specimen_voucher', 'country', 'collection_date', 'lat_lon', 'note', 'collected_by', 'isolate', 'pop_variant')) {
  if(cores > 1 & Sys.info()['sysname'] == 'Windows') warning("Multicore is only supported on mac and linux for right now")
  require(ape)
  require(XML)
  require(parallel)
  nRecords <- length(xml_file$doc$children$INSDSet)
  columns <- c('NCBI_accession', 'seq_length','strandedness','moltype','topology','division',
               'update_date','create_date','definition','primary_accession','accession_version',
			   ' otherseq_IDS','seq_source','organism','taxonomy','references','feature_table',
			   'qualifiers1','generegion','Full_sequence', 'authors') ## not needed currently, but might be useful for making the code more flexible
  #open(file('fails.v5.txt', open = 'a'))
  get.a.row <- function(recordNumber) {
    #message(paste('doing record', recordNumber))
    dat <- xml_file$doc$children$INSDSet[[recordNumber]]
    featuresL <- length(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']])
	featuresOut <- character(featuresL)
	for(i in seq(featuresL)) {
	  if(length(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']][[i]]) != 2) next
          featuresOut[i] <- xmlValue(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']][[i]][[2]])
	  names(featuresOut)[i] <- xmlValue(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']][[i]][[1]])
	  }
	readableFeatures = paste(names(featuresOut), featuresOut, collapse = '|', sep = "_:_")
	featuresOutV <- featuresOut[qualsToUse]
	names(featuresOutV) <- qualsToUse
	out <-
	     c(NCBI_accession = try(xmlValue(dat[['INSDSeq_locus']]), silent = T),
             seq_length = try(xmlValue(dat[["INSDSeq_length"]]), silent = T),
             strandedness = try(xmlValue(dat[["INSDSeq_strandedness"]]), silent = T),
             moltype = try(xmlValue(dat[["INSDSeq_moltype"]]), silent = T),
             topology  = try(xmlValue(dat[["INSDSeq_topology"]]), silent = T),
             division = try(xmlValue(dat[["INSDSeq_division"]]), silent = T),
             update_date = try(xmlValue(dat[["INSDSeq_update-date"]]), silent = T),
             create_date = try(xmlValue(dat[["INSDSeq_create-date"]]), silent = T),
             definition = try(xmlValue(dat[["INSDSeq_definition"]]), silent = T),
             primary_accession = try(xmlValue(dat[["INSDSeq_primary-accession"]]), silent = T),
             accession_version = try(xmlValue(dat[["INSDSeq_accession-version" ]]), silent = T),
             otherseq_IDS = try(xmlValue(dat[["INSDSeq_other-seqids" ]]), silent = T),
             seq_source = try(xmlValue(dat[["INSDSeq_source"]]), silent = T),
             organism = try(xmlValue(dat[["INSDSeq_organism"]]), silent = T),
             taxonomy = try(xmlValue(dat[["INSDSeq_taxonomy"]]), silent = T),
             references = try(xmlValue(dat[["INSDSeq_references"]]), silent = T),
             feature_table = try(xmlValue(dat[["INSDSeq_feature-table"]]), silent = T),
             qualifiers1 = readableFeatures,  #part of feature tables, flattened out for readability
             generegion = try(xmlValue(dat[["INSDSeq_feature-table"]][[2]][['INSDFeature_quals']][['INSDQualifier']][['INSDQualifier_value']]), silent = T), ##within feature_table node
             Full_sequence = ifelse(includeSeqs, try(xmlValue(dat[["INSDSeq_sequence"]]), silent = T), ''),
	     authors = try(xmlValue(dat[['INSDSeq_references']][[1]][['INSDReference_authors']]), silent = T)
	     ) # close c
    #if(class(out) == 'try-error') out <- c(xmlValue(dat[['INSDSeq_locus']]), 'failed', rep(0, 19))
    #if(any(substr(out, 1, 5) == 'Error')) writeLines(xmlValue(dat[['INSDSeq_locus']]), con = file('fails.txt', open = 'a'))
    out <- c(out, featuresOutV)
    # write.csv(out, paste(xmlValue(dat[['INSDSeq_locus']]), 'csv', sep = '.'))
    return(out)
    }
  if(!is.na(do[1])) xmlMat <- t(mcmapply(get.a.row, xml_file$doc$children$INSDSet[do]))
  # else xmlMat <- (mcmapply(get.a.row, xml_file$dfor individuals sequencedoc$children$INSDSet, mc.cores = cores))
  # else xmlMat <- t(mclapply(xml_file$doc$children$INSDSet, get.a.row, mc.cores = cores)) # not simplyifying for sake of
  else xmlMat <- do.call(rbind, mclapply(seq(nRecords), get.a.row, mc.cores = cores))
  #debug(get.a.row)
  #xtemp <- lapply(seq(nRecords), get.a.row)
  #xmlMat <- do.call(rbind, xtemp)
  if(parse.specimens) xmlMat <- cbind(xmlMat, parse.specimen(xmlMat))
  #close(file('fails.v5.txt'))
  return(xmlMat)
}

### XML GENBANK PARSER
### To Parse XML files downloaded from Geneious search of Genbank in order to get additional info such as voucher info, etc
### 2013-02: Marlene Hahn, as part of Carex project.
### 2015-03-18: Hipp and Hahn, updating to use multicores
### 2015-03-24: use only named references and pull out more needed elements
### 2015-03-25: get the desired voucher features directly from the xml

###Step 1: Save file as XML in Geneious; make sure you know how many files were downloaded 
###because you will need to enter that # in function below (usually part of filename in Geneious export)
	### libary(XML)
	### xml_file <- xmlTreeParse(file.choose())
 
 ###Issues see if I can get the program to parse out more detailed information including gene region and voucher info.
		### problem- voucher is in qualifiers in feature_table node, but not at a consistant node.
		##Still has issues parsing to the voucher level- use spliting _metadata_genbank_tables.r function to parse out this info.)


parse.INSDSeq = function(xml_file, do = NA, includeSeqs = F, cores = 1, parse.specimens = T,
                         qualsToUse = c('specimen_voucher', 'DNAisolate', 'pop_variant', 'collection_date', 'lat_lon', 'note', 'collected_by')) { 
  if(cores > 1 & Sys.info()['sysname'] == 'Windows') warning("Multicore is only supported on mac and linux for right now")
  require(ape)
  require(XML)
  require(parallel)
  nRecords <- length(xml_file$doc$children$INSDSet)
  columns <- c('NCBI_accession', 'seq_length','strandedness','moltype','topology','division',
               'update_date','create_date','definition','primary_accession','accession_version',
			   ' otherseq_IDS','seq_source','organism','taxonomy','references','feature_table',
			   'qualifiers1','generegion','Full_sequence', 'authors') ## not needed currently, but might be useful for making the code more flexible
  get.a.row <- function(dat) {
    featuresL <- length(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']])
	featuresOut <- character(featuresL)
	for(i in seq(featuresL)) {
	  featuresOut[i] <- xmlValue(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']][[i]][[2]])
	  names(featuresOut)[i] <- xmlValue(dat[['INSDSeq_feature-table']][[1]][['INSDFeature_quals']][[i]][[1]])
	  }
	readableFeatures = paste(names(featuresOut), featuresOut, collapse = '|', sep = "_:_")
	featuresOutV <- featuresOut[qualsToUse]
	names(featuresOutV) <- qualsToUse
	out <- try(
	         c(NCBI_accession = xmlValue(dat[['INSDSeq_locus']]),
             seq_length = xmlValue(dat[["INSDSeq_length"]]),
             strandedness = xmlValue(dat[["INSDSeq_strandedness"]]),
             moltype = xmlValue(dat[["INSDSeq_moltype"]]),
             topology  = xmlValue(dat[["INSDSeq_topology"]]),
             division = xmlValue(dat[["INSDSeq_division"]]),
             update_date = xmlValue(dat[["INSDSeq_update-date"]]),
             create_date = xmlValue(dat[["INSDSeq_create-date"]]),
             definition = xmlValue(dat[["INSDSeq_definition"]]),
             primary_accession = xmlValue(dat[["INSDSeq_primary-accession"]]),
             accession_version = xmlValue(dat[["INSDSeq_accession-version" ]]),
             otherseq_IDS = xmlValue(dat[["INSDSeq_other-seqids" ]]),
             seq_source = xmlValue(dat[["INSDSeq_source"]]),
             organism = xmlValue(dat[["INSDSeq_organism"]]),
             taxonomy = xmlValue(dat[["INSDSeq_taxonomy"]]),
             references = xmlValue(dat[["INSDSeq_references"]]),
             feature_table = xmlValue(dat[["INSDSeq_feature-table"]]),
             qualifiers1 = readableFeatures,  #part of feature tables, flattened out for readability
             generegion = xmlValue(dat[["INSDSeq_feature-table"]][[2]][['INSDFeature_quals']][['INSDQualifier']][['INSDQualifier_value']]), ##within feature_table node
             Full_sequence = ifelse(includeSeqs, xmlValue(dat[["INSDSeq_sequence"]]), ''),
			 authors = xmlValue(dat[['INSDSeq_references']][[1]][['INSDReference_authors']])
	         ), # close c
		    silent = T) # close try
    if(class(out) == 'try-error') out <- c(xmlValue(dat[['INSDSeq_locus']]), 'failed', rep(0, 19))
	out <- c(out, featuresOutV)
	return(out)
	}
  if(!is.na(do[1])) xmlMat <- t(mcmapply(get.a.row, xml_file$doc$children$INSDSet[do]))
  else xmlMat <- t(mcmapply(get.a.row, xml_file$doc$children$INSDSet, mc.cores = cores))
  if(parse.specimens) xmlMat <- cbind(xmlMat, specimen.match(xmlMat))
  return(xmlMat)
}















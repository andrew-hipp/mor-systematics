### XML GENBANK PARSER
### To Parse XML files downloaded from Geneious search of Genbank in order to get additional info such as voucher info, etc
### 2013-02: Marlene Hahn, as part of Carex project.
### 2015-03-18: Hipp and Hahn, updating to use multicores

###Step 1: Save file as XML in Geneious; make sure you know how many files were downloaded 
###because you will need to enter that # in function below (usually part of filename in Geneious export)
	### libary(XML)
	### xml_file <-  xmlTreeParse(file.choose())
 
 ###Issues see if I can get the program to parse out more detailed information including gene region and voucher info.
		### problem- voucher is in qualifiers in feature_table node, but not at a consistant node.
		##Still has issues parsing to the voucher level- use spliting _metadata_genbank_tables.r function to parse out this info.)


parse.INSDSeq = function(xml_file, file_length, do = NA, cores = 1) {  ##filelength = # of specimens in export
  if(cores>1) warning("Multicore is only supported on mac and linux for right now")
  require(ape)
  require(XML)
  require(parallel)
  nRecords <- length(xml_file$doc$children$INSDSet)
  columns <- c('NCBI_accession', 'seq_length','strandedness','moltype','topology','division',
               'update_date','create_date','definition','primary_accession','accession_version',
			   ' otherseq_IDS','seq_source','organism','taxonomy','references','feature_table',
			   'qualifiers1','generegion','Full_sequence') ## not needed currently, but might be useful for making the code more flexible
  get.a.row <- function(dat) {
    out <- c(NCBI_accession = xmlValue(dat[['INSDSeq_locus']]),
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
             qualifiers1 = xmlValue(dat[["INSDSeq_feature-table"]][[1]][[5]]),  #part of feature tables
             generegion = xmlValue(dat[["INSDSeq_feature-table"]][[2]][[5]][['INSDQualifier']][['INSDQualifier_value']]), ##within feature_table node
             Full_sequence = xmlValue(dat[["INSDSeq_sequence"]])
	         )
    return(out)
	}
  if(!is.na(do)) xmlMat <- t(mcmapply(get.a.row, xml_file$doc$children$INSDSet[do]))
  else xmlMat <- t(mcmapply(get.a.row, xml_file$doc$children$INSDSet))
  # write.table(table, file = 'Metadata_tablev2.txt', sep = "|")  
  return(table)
}















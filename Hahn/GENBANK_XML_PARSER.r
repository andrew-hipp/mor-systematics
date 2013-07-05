### XML GENBANK PARSER
### To Parse XML files downloaded from Geneious search of Genbank in order to get additional info such as voucher info, etc
### Marlene Hahn February 2013, as part of Carex project.

###Step 1: Save file as XML in Geneious; make sure you know how many files were downloaded 
###because you will need to enter that # in function below (usually part of filename in Geneious export)
	### libary(XML)
	### xml_file <-  xmlTreeParse(file.choose())
 
 ###Issues see if I can get the program to parse out more detailed information including gene region and voucher info.
		### problem- voucher is in qualifiers in feature_table node, but not at a consistant node.
		##Still has issues parsing to the voucher level- use spliting _metadata_genbank_tables.r function to parse out this info.)


makeDataFrame_XML_voucher = function(xml_file, file_length) {  ##filelength = # of specimens in export
  require(ape)
  require(XML)
  NCBI_accession <- character()
  seq_length <- numeric()
  strandedness <- character()
  moltype <- character()
  topology <- character()
  division <- character()
  update_date <- character()
  create_date  <- character()
  definition  <- character()
  primary_accession  <- character()
  accession_version <- character()
  otherseq_IDS <- character()  
  #keywords <- character()  ###error with different # of rows; is this a problem if there is null information??
  seq_source <- character()
  organism <- character()
  taxonomy <- character() 
  references <- character()
  feature_table <- character() 
  qualifiers1 <- character() 
  generegion <- character()
  Full_sequence <- character()
  for(i in 1:file_length) {
   NCBI_accession[i] <- xmlValue(xml_file$doc$children$INSDSet[[i]][['INSDSeq_locus']])
   seq_length[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_length"]])
   strandedness[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_strandedness"]])
   moltype[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_moltype"]])
   topology [i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_topology"]])
   division[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_division"]])
   update_date[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_update-date"]])
   create_date[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_create-date"]])
   definition[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_definition"]])
   primary_accession[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_primary-accession"]])
   accession_version[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_accession-version" ]])
   otherseq_IDS[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_other-seqids" ]])
   #keywords[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_keywords"]])
   seq_source[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_source"]])
   organism[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_organism"]])
   taxonomy[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_taxonomy"]])
   references[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_references"]])
   feature_table[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_feature-table"]])
   qualifiers1[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_feature-table"]][[1]][[5]])  #part of feature tables
   generegion[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_feature-table"]][[2]][[5]][['INSDQualifier']][['INSDQualifier_value']]) ##within feature_table node
   Full_sequence[i] <-xmlValue(xml_file$doc$children$INSDSet[[i]][["INSDSeq_sequence"]])
	}
	table <- data.frame(NCBI_accession, seq_length,strandedness, moltype, topology, division, update_date, create_date, definition, primary_accession, accession_version, otherseq_IDS, feature_table, qualifiers1, generegion, seq_source, organism, taxonomy, references, Full_sequence)
	write.table(table, file = 'Metadata_tablev2.txt', sep = "|")  
 return(head(table)) 
}















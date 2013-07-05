## functions for working through usearch files
## tailored to cariceae output currently
## Hahn and Hipp, January 2013

## source("usearch.r") to use functions in this document in R

###Step 1: In geneious, query GenBank for all sequences that contain Organism = Carex, Schoenoxiphium, Cymophyllus, Kobresia, Uncinia or Vesicarex. Export sequences as fasta file and tab delimited file.

###Step 2:read in fasta and tab delim files into R using following R commands.
#All_cariceae.fasta <- read.dna(file.choose(), format = 'fasta')    ###Chose fasta file from step 1
#All_carriceae.tabdelim <- read.delim('All_carriceae_2013_1_22.export')


###Step 3: removeMicrosatsSeq function which removes all the microsats and ISSR sequences within fasta file from geneious, creates new file with date appended on.
removeMicrosatsSeq = function(All_cariceae.fasta) {
require(ape)
write.dna(All_cariceae.fasta[-c(grep('microsat', names(All_cariceae.fasta),ignore.case=T), grep('issr', names(All_cariceae.fasta),ignore.case=T))], file = paste('test_all_cariceae.noSSR.',format(Sys.time(),"%Y-%m-%d"),'.fas', sep =''), format = 'fasta')
}

###Step 4: Use new .noSSR file in Usearch to create clusters, using the following command prompts:
#> usearch -sortbylength test_all_cariceae.noSSR.date.fasta  -output  test_all_cariceae.noSSR.date.fasta_seqs_sorted.fasta  ###sorts sequences based upon size
#> usearch -cluster_smallmem test_all_cariceae.noSSR.date.fasta_seqs_sorted.fasta -id 0.3 -centroids test_all_cariceae.noSSR.date.fasta_cluster.id3.fasta -strand both -uc test_all_cariceae.noSSR.date.fasta_cluster.id3.uc  

###Step 5: Put cluster analysis back into R.
#test_cluster_cariceae <- read.delim('test_all_cariceae.noSSR.date.fasta_cluster.id3.uc', header = F)

###Step 6: Add headers to columns in cluster file
AddHeadersToClusters = function(test_cluster_cariceae){               
names(test_cluster_cariceae) <- uc.titles
head(test_cluster_cariceae)
}

## function to give appropriate names (species+accession) to output cluster fasta files- run this prior to makeAllData.accessions
cleanFastaNames <- function(fasta = All_cariceae.fasta, metadata = All_carriceae.tabdelim) {
  names(fasta) <- gsub(" ", "_", paste(metadata$Organism, metadata$Accession))
  fasta
  }
  
###function for creating text files of descriptions of sequences clustered by usearch 
makeAllData.summaries = function(uc.table, metadata = All_carriceae.tabdelim) {
  datOut <- matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, c('clusterNumber', 'description.first', 'clusterLength')))
  for(i in 0:max(as.numeric(uc.table$clusterNumber))) {
    tempOut <- as.character(metadata$Description[which(metadata$Name %in% uc.table$querySeq[which(uc.table$clusterNumber == i)])])
	writeLines(tempOut, con = paste("cluster", i, ".trial_", format(Sys.time(), "%Y-%m-%d"), ".txt", sep = ""))
	datOut <- rbind(datOut, c(i, tempOut[1], length(tempOut)))
	}
  return(datOut)
}


###function for creating fasta files from the cluster analysis. note these final files don't seem to readable by muscle program (should the code be changed to .fa ??)
makeAllData.accession = function(uc.table = clusters.v2, metadata = All_carriceae.tabdelim, fasta = All_cariceae.cleanNames.fasta, basedir = './outdat/') {
  datOut <- matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, c('clusterNumber', 'description.first', 'clusterLength')))
  require(ape)
  for(i in 0:max(as.numeric(uc.table$clusterNumber))) {
    tempOut <- as.character(metadata$Accession[which(metadata$Name %in% uc.table$querySeq[which(uc.table$clusterNumber == i)])])
	writeLines(tempOut, con = paste("accessions_in_cluster", i, ".trial_", format(Sys.time(), "%Y-%m-%d"), ".txt", sep = ""))
	write.dna(fasta[which(metadata$Accession %in% tempOut)], file = paste('cluster_', i, "_", format(Sys.time(),"%Y-%m-%d"),'.fas', sep =''), format = 'fasta')
	datOut <- rbind(datOut, c(i, tempOut[1], length(tempOut)))
	}
  return(datOut)
}

##function to see what strands (if any) need to be reverse-complemented in a particular cluster (enter in the cluster number you're analyzing)
ID_reverse_strands = function(clusternum, metadata = clusters.v2) {
  temp <-metadata[which(metadata$strand == "-"),]
  rev_strand_clust <-temp[which(temp$clusterNumber == clusternum),]
  rev_strand_clust
}


###after doing the makeAllDAta.accession function- take files and muscle them for alignments (using command prompt and muscle program.) (use Muscle_Files.bat file to rerun these)
##muscle -in <inputfile> -out <outputfile>
##muscle -in its.hopefuls.2013-02-28.fa -out its.hopefuls.2013-02-28_muscled.fasta

exportAllFastaFiles <- function(uc.table, fastaIn) {}

makesNewTitles <- function(titleCodes, someFormula, someDataSource) {} 


get.seqs.from.ncbi <- function(coll.dat, dna.title, collector.name = 'noname', concat = c('ITS1', 'ITS2'), dna = cariceae.ncbi, write.files = T, clean.spaces = T, concatChars = paste(rep('N', 100), collapse = ''), silent = TRUE) {
## acc = accession vector to feed to NCBI
## dna.title = gene name (e.g., "ITS", "ETS"), and header to draw data from
## for concat to work, you have to have the string "concat" at the beginning of the gene header (e.g., "concat_1_2" in the ITS column in order to read ITS1 and ITS2 separately and concatenate)
  require(ape)
 
  acc <- coll.dat[[dna.title]]
  which.concat <- grep('concat', acc)
  concat.flag <- length(which.concat) > 0
  if(clean.spaces) acc <- gsub(" ", "", acc)
  sp <- coll.dat$sciname.edited
  work.rows <- acc %in% row.names(dna)
  if(concat.flag) {
    acc2 <- coll.dat[concat]
	if(clean.spaces) acc2 <- apply(acc2, 2, function(x) gsub(' ', '', x))
	concat.rows <- T | work.rows
	for(i in 1:dim(acc2)[2]) concat.rows <- concat.rows & (acc2[, i] %in% row.names(dna)) # throws out any rows in which any concat regions don't work
	work.rows <- work.rows | concat.rows # the row is usable either if the accession number is in DNA title or if both fo the concat regions are
	}
  works <- acc[which(work.rows)] # a label list of which worked
  fails <- acc[which(!work.rows)] # a label list of which failed
  if(concat.flag) {
    out <- structure(character(length(works)), names = works)
    out[grep('concat', works)] <- paste(dna[acc2[concat.rows, 1], 'Full_sequence'], rep(concatChars, length(which.concat)), dna[acc2[concat.rows, 2], 'Full_sequence'], sep = '')
	out[-grep('concat', works)] <- dna[works[-grep('concat', works)], 'Full_sequence']
	}
  else out <- dna[works, 'Full_sequence']
  names(out) <- gsub(" ", "_", as.character(apply(coll.dat[which(work.rows), c('sciname.edited', 'DNA')], 1, paste, collapse = "_"))) 
  if(write.files) {
    write.table(coll.dat[which(work.rows), ], paste(collector.name, dna.title, format(Sys.time(), "ncbiCodes_%Y-%m-%d.tsv"), sep = '_'), sep = '\t')
	write.table(coll.dat[-which(work.rows), ], paste(collector.name, dna.title, format(Sys.time(), "failed_%Y-%m-%d.tsv"), sep = '_'), sep = '\t')
	write.dna(as.DNAbin(lapply(strsplit(out, ''), tolower)), paste(collector.name, dna.title, format(Sys.time(), "dna_%Y-%m-%d.fas"), sep = '_'), format = 'fasta')
	}
  if(silent) return(invisible(out))
  return(out)
  }

do.dna.exports <- function(do = list(Gehrke = c('ITS.original', 'trnLF.original', 'ITS.corrected', 'trnLF.corrected'),
                                     STARR.FORD = c('ITS.original', 'ITS.corrected'),
									 Luceno = 'ITS.original', 'ITS.corrected', 'ETS.original', 'ETS.corrected', 'trnLF.original', 'trnLF.corrected'
									 )) {
  for(i in names(do)) {
    for(j in do[[i]]) {
	  get.seqs.from.ncbi(collectors[[i]], j, i)
	  }}}

orig.dna.exports <- function() {	  
get.seqs.from.ncbi(collectors.dna$escudero, 'ITS', 'escudero')
get.seqs.from.ncbi(collectors.dna$escudero, 'trnK', 'escudero')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'ITS', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'ETS', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'X5IGSf', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'trnT.trnL', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'trnL.trnF', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'matK', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'rpoC1', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'rbcL', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'rpoB', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'atpFH', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'psbKI', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$starr.ford, 'trnH.psbA', 'starr.ford')
get.seqs.from.ncbi(collectors.dna$gehrke, 'ITS', 'gehrke')
get.seqs.from.ncbi(collectors.dna$gehrke, 'trnLF', 'gehrke')
get.seqs.from.ncbi(collectors.dna$luceno, 'trnLF', 'luceno')
get.seqs.from.ncbi(collectors.dna$luceno, 'ITS', 'luceno')
get.seqs.from.ncbi(collectors.dna$luceno, 'ETS', 'luceno')
get.seqs.from.ncbi(collectors.dna$luceno, 'rps16', 'luceno')
}
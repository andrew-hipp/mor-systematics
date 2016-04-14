## functions to make gene matrices based on specimen data parsed out of NCBI using specimen.match.R functions
## Andrew Hipp, June 2013
## 2013-07-05: separated into a separate file

make.shared.gene.matrix <- function(dat = geneMat.2013.06.13, loci = dimnames(dat)[[2]][!dimnames(dat)[[2]] %in% c("orgs","ncbiAcc","numberOfOrgs","numberOfAccessions","numberOfSequences")]) {
## takes output from make.gene.matrix
  out <- matrix(NA, length(loci), length(loci), dimnames = list(loci, loci))
  for(i in 1:length(loci)) {
    for(j in 1:length(loci)) {
	  out[loci[i], loci[j]] <- sum(ifelse(dat[[loci[i]]] == 0, 0, 1) + ifelse(dat[[loci[j]]] == 0, 0, 1) == 2)
	  } # close j
	} # close i
  return(out)
  }

orphan.loci <- function(x)  {
## takes output from make.shared.gene.matrix
  names(which(colSums(x) == 0))
  }

orphan.individuals <- function(dat) {
## takes output from make.gene.matrix
  shared.loci <- make.shared.gene.matrix(dat)
  not.orphan.loci <- names(which(colSums(shared.loci) > 0))
  dat.no.orphan.loci <- dat[, not.orphan.loci]
  not.orphan.individuals <- apply(dat.no.orphan.loci, 1, function(x) sum(x) > 0)
  orphan.individuals <- apply(dat.no.orphan.loci, 1, function(x) sum(x) == 0)
  out <- list(NCBI_accession = dat$ncbiAcc, taxon = dat$orgs, orphans = orphan.individuals, notOrphans = not.orphan.individuals)
  return(out)
  }

make.fasta.files <- function(geneMatrix = geneMat.2013.06.17, seqDat = cariceae.2013.02.28, genes = top12, outdir = paste('fasta.out.', paste(sample(letters,5), collapse = ''), format(Sys.time(), "%Y-%m-%d"), sep = ''), maxtax = 1, treat.multiples = c('discard', 'takefirst'), batchfile = 'muscle') {
## makes fasta files for individuals having the genes indicated
## creates a log file with all individuals used for each gene
## arguments:
##  geneMatrix: output from make.gene.matrix
##  seqDat: sequence data from NCBI
##  genes: a vector of genes to use
##  outdir: directory to write to
##  maxtax: maximum number of taxon names allowable per individual
##  treat.multiples: what to do for individuals that have multiple accessions for a given gene; currently either discards ('discard') them or takes the first NCBI accession ('takefirst')

  if(!outdir %in% dir()) dir.create(outdir)
  geneMatrix <- geneMatrix[geneMatrix$numberOfOrgs <= maxtax, ]
  geneMatrix$seqLabels <- tidyNames(paste(geneMatrix$orgs, row.names(geneMatrix)))
  # make individual gene files
  for(i in genes) {
	# if(i == 'rbcL') browser()
	out.log <- cbind(geneMatrix[geneMatrix[[i]] != '', c('orgs', 'seqLabels', i)], whatToDo = '', comments = '')
	multiples <- grep('|', out.log[, i], fixed = T)
	if(treat.multiples[1] == 'discard') {if(length(multiples) != 0) out.log <- out.log[-multiples, ]}
	if(treat.multiples[1] == 'takefirst') out.log[multiples, i] <- sapply(out.log[multiples, i], function(x) strsplit(x, "|", fixed = T)[[1]][1])
	seqs <- seqDat$Full_sequence[match(out.log[, i], seqDat$NCBI_accession)]
	writeLines(paste(">", out.log$seqLabels, "\n", seqs, sep = ''), paste(outdir, '/', i, format(Sys.time(), ".%Y-%m-%d.fas"), sep = ''))
	write.csv(out.log, paste(outdir, '/', i, '.logfile.', format(Sys.time(), "%Y-%m-%d.csv"), sep = ''))
  }

  ## make batch file
  if(batchfile == 'muscle') {
    dir.create(paste(outdir, '/muscle', sep = ''))
	# files <- dir(outdir, full = TRUE, patt = '.fas')
	littleFiles <- dir(outdir, full = FALSE, patt = '.fas')
	writeLines(paste('muscle3.8.31_i86win32 -in ', paste('../', littleFiles, sep = ''), ' -out ', paste(littleFiles, '.muscled.fas', sep = ''), sep = ''), paste(outdir, '/muscle/muscle.bat', sep = ''))
	}
  return('done!')
  }

refine.fasta <- function(basedir = choose.dir()) {
## goes through the log files in a directory and makes changes according to these codes:
##  D = delete
##  RC = reverse and complement
  require(Biostrings)
  require(ape)
  all.fasta <- sort(dir(basedir, full = T, patt = '.fas'))
  all.logs <- sort(dir(basedir, full = T, patt = 'logfile'))
  files.out <- NULL
  for(i in 1:length(all.logs)) {
    toDo <- read.csv(all.logs[i], as.is = TRUE)
	if(!"whatToDo" %in% names(toDo)) next
	toDo <- toDo[toDo$whatToDo != '', ]
	if(nrow(toDo) == 0) next
	fasta <- read.dna(all.fasta[i], format = 'fasta', as.character = TRUE)
	for(j in 1:nrow(toDo)) {
	  if(toDo[j, 'whatToDo'] == 'RC') fasta[[toDo$seqLabels[j]]] <- strsplit(as.character(reverseComplement(DNAString(paste(fasta[[toDo$seqLabels[j]]], collapse = '')))), '')
	  if(toDo[j, 'whatToDo'] == 'D') fasta <- fasta[-(which(names(fasta) == toDo$seqLabels[j]))]
	  }
	file.out <- paste(all.fasta[i],'.cleaned.fas', sep = '')
	files.out <- c(files.out, file.out)
	write.dna(fasta, file.out, format = 'fasta')
	}
	writeLines(paste('muscle3.8.31_i86win32 -in ', files.out, ' -out ', paste(files.out, '.muscled.fas', sep = ''), sep = ''), paste(basedir, '/muscle/muscle.bat', sep = ''))
  return('done!')
  }

make.unique.vouchers <- function(metadata, voucherFormula = c("Primary.collector.last.name", "Collector.number", "isolate", "CollectionNumber", "Collection")) {
      vouchers.ln.cn <- apply(metadata[voucherFormula], 1, function(x) paste(tidyName(x[!is.na(x)]), collapse = ''))
      names(vouchers.ln.cn) <- metadata$NCBI_voucher # note that NCBI_accession is not actually the accession number used in NCBI! The preferred accession number is primary_accession
      return(vouchers.ln.cn)
      }


make.gene.matrix <- function(metadata, locusCol = 'cleanedGeneRegion', vouchersCol = 'newLabels', ncbiCol = 'NCBI_accession', orgsCol = 'organism', logerrors = TRUE, verbose = FALSE) {
## take vouchers and regions to make a matrix we can use
## Arguments:
##  metadata = parsed data from NCBI genbank
##  loci = translation from verbatim gene regions (NCBI) to a cleaned up gene regions name; relevant columns are "verbatim" and "clean"; every genbank row should correspond with an entry in the verbatim column
##  spmCol = the name of the specimen label column, if it has been creaated; if not, voucherFormula is used to create it

  missingVouchers <- which(gsub(" ", "", metadata[[vouchersCol]], fixed = TRUE) == "")
  uniqueLoci <- unique(sort(as.character(metadata[[locusCol]])))
  uniqueVouchers <- unique(sort(as.character(metadata[[vouchersCol]])))
  out <- matrix('', nrow = length(uniqueVouchers), ncol = length(uniqueLoci), dimnames = list(uniqueVouchers, uniqueLoci))

  # vector of organisms
  orgs <- sapply(uniqueVouchers, function(x) paste(as.character(unique(metadata[[orgsCol]][metadata[[vouchersCol]] %in% x])), collapse = "|"))

  # vector of NCBI accessions
  ncbiAcc <- sapply(uniqueVouchers, function(x) paste(as.character(unique(metadata[[ncbiCol]][metadata[[vouchersCol]] %in% x])), collapse = "|"))

  # populate the matrix
  meta.orig <- metadata
  if(length(missingVouchers) > 0 ) metadata <- metadata[-missingVouchers,]
  # browser()
  for (i in 1:dim(metadata)[1]) {
    if(!any(is.na(metadata[i, c('cleanedGeneRegion', 'cleanedVoucher')]))) {
	  if(verbose) message(paste('doing', i))
	  out[metadata[i, vouchersCol], metadata[i, locusCol]] <- ifelse(out[metadata[i, vouchersCol], metadata[i, locusCol]] == '',
	                                                                 as.character(metadata[i, ncbiCol]),
																	 paste(out[metadata[i, vouchersCol], metadata[i, locusCol]], metadata[i, ncbiCol], sep = '|')
																					 )
	  }
	else(message(paste("Row", i, "of your metadata table seems to have a problem")))
	} # close i
  numberOfOrgs <- sapply(strsplit(as.character(orgs), "|", fixed = T), length)
  numberOfAccessions <- sapply(strsplit(as.character(ncbiAcc), "|", fixed = T), length)
  # numberOfSequences <- apply(out, 1, sum)
  out <- cbind(orgs, ncbiAcc, numberOfOrgs, numberOfAccessions, as.data.frame(out))

  if(logerrors) write.csv(metadata[missingVouchers, ], paste('missingVouchers.log.', paste(sample(letters,5), collapse = ''), '.csv', sep = ''))
  
  # class(out) <- c('geneMat', 'matrix')
  return(out)
  }

gene.matrix.stats <- function(mat = geneMat.2013.06.13, outfile = paste('geneStats.', paste(sample(letters, 5), collapse = ''), '.txt', sep = '')) {
## describes the gene matrix
  loci <- dimnames(mat)[[2]][!dimnames(mat)[[2]] %in% c("orgs","ncbiAcc","numberOfOrgs","numberOfAccessions","numberOfSequences")]
  out <- file(outfile, open = 'a')
  writeLines(paste('Unique taxa:', length(unique(mat$orgs))), con = out)
  writeLines(paste('Unique taxa, only one organism:', length(unique(mat$orgs[mat$numberOfOrgs == 1]))), con = out)
  writeLines('Unique taxa for each gene:\n----------------------', con = out)
  temp <- sapply(loci, function(x) length(unique(mat$orgs[mat[[x]] >0])))
  writeLines(paste(loci[order(temp, decreasing = TRUE)], ": ", sort(temp, decreasing = TRUE), sep = ""), con = out)
  close(out)
  }

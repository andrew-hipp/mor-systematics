match.dna <- function(dna.list = dna.renamed, seq.tables = DNA[c('MOR', 'SHUREN', 'WSU')], source.tables = collectors, taxon.table = gcg, 
                      specimen.table = specimens, dna.list.names = names(dna.list), source.tables.names = substr(names(source.tables), 1, 3),
                      add.parents = T, ...) {
## dna.tables = a list of dna matrices
## seq.tables = a list of tables indexing the dna.tables
## source.tables = list of tables of folks who said they'd provide material
## taxon.table = an SP2 taxonomic table
## specimen.table = a table linking DNA to specimens
## ... = arguments passed on to subset.checklist
  
  rows <- sort(unique(taxon.table[subset.checklist(taxon.table, ...), 'Term.name']))
  if(length(dna.list.names) != length(dna.list)) dna.list.names = paste('dna.dataset', seq(length(dna.list)), sep = '.')
  cols <- c(dna.list.names, 'sequenced', 'who.has.it')
  out <- matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
  dna.taxa <- lapply(dna.list, function(x) get.details.from.DNA(row.names(x), seq.tables, specimen.table))
  names(dna.taxa) <- dna.list.names
  for(i in dna.list.names) {
    missing.data <- row.names(dna.list[[i]])[gsub("[NA ]", "", dna.taxa[[i]]) == '']
	write.csv(cbind(DNA.CODE = missing.data, specimen.table[get.details.from.DNA(missing.data, seq.tables, specimen.table, give = 'specimen'), ]), paste(i, '.missingData.csv', sep = ''))
	}
  for(i in rows) {
    for(j in dna.list.names) {
	  out[i, j] <- sum(i == dna.taxa[[j]], na.rm = T)
	  }
	out[i, 'who.has.it'] <- paste(sort(unique(names(which(sapply(collectors, function(x) i %in% x$sciname.edited))))), collapse = ', ')
	out[i, 'sequenced'] <- any(as.numeric(out[i, dna.list.names]) > 0)
	}
  if(add.parents) out <- cbind(parent.taxon = sapply(row.names(out), get.parent, USE.NAMES = FALSE), out)
  write.csv(out, 'out.temp.csv')
  out <- rbind(out, taxa.sequenced = apply(out, 2, function(x) sum(! x %in% c('FALSE', '0', '')))) 
  out
  }

get.details.from.DNA <- function(x, seq.tables, specimen.table, spm.patt = 'Specimen', extract.patt = 'TUBE', give = 'taxon', clip.subspp = TRUE) {
## returns for a character vector x the taxa corresponding to the extraction codes
  if(class(seq.tables) == 'list') seq.tables <- do.call(rbind, seq.tables)
  spm.col <- grep(spm.patt, names(seq.tables))
  ext.col <- grep(extract.patt, names(seq.tables))
  if(any(duplicated(seq.tables[[spm.col]]))) warning('There is a duplicate specimen number in the tables provided; at least one specimen is reached by two extractions')
  if(any(duplicated(seq.tables[[ext.col]]))) warning('There is a duplicate extraction code in the tables provided')
  x.spm <- seq.tables[match(tidyName(x), tidyName(seq.tables[[ext.col]])), spm.col]
  if(give == 'taxon') {
    x.taxon <- specimen.table[x.spm, "TAXA.Current_determination"]
    if(clip.subspp) x.taxon <- sapply(x.taxon, function(x) paste(strsplit(x, " ")[[1]][1:2], collapse = ' '), USE.NAMES = FALSE)
    return(x.taxon)
	}
  if(give == 'specimen') return(x.spm)
  }

get.parent <- function(x, spm = gcg) spm$Term.name[spm$GUID == spm$Parent.GUID[spm$Term.name == x]][1]
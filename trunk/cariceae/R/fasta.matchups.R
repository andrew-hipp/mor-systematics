match.dna <- function(dna.list, seq.tables, source.tables, taxon.table, specimen.table, 
                      dna.list.names = names(seq.tables), source.tables.names = substr(names(source.tables), 1, 3), ...)
## dna.tables = a list of dna matrices
## seq.tables = a list of tables indexing the dna.tables
## source.tables = list of tables of folks who said they'd provide material
## taxon.table = an SP2 taxonomic table
## specimen.table = a table linking DNA to specimens
## ... = arguments passed on to subset.checklist
  names.to.use <- taxon.table[subset.checklist(taxon.table, ...), 'Term.name']
  
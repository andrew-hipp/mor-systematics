  if(software[1] == 'paup') {
    write.nexus(treeset, file = paste(filebase, '.rtrees.tre', sep = ''))
    write.nexus(x, file = paste(filebase, '.optimal.tre', sep = ''))
	}


  ## write paup block
  if(software[1] == 'paup') {
    paup.out <- file(description = paste(filebase, '.paupBlock.nex', sep = ''), open = "w")
    writeLines("#NEXUS", paup.out)
    writeLines(paste("[PAUP block written using genTrees function in R,", date(), "]"), paup.out)
    writeLines("\nBEGIN PAUP;\n", paup.out)
    writeLines("  exclude constant /only;", paup.out)
    writeLines("  increase = auto autoclose = yes;", paup.out)
    writeLines("  cleartrees;", paup.out)
    writeLines(paste("  gettrees file =", paste(filebase, '.optimal.tre', sep = ''), ";"), paup.out)
    writeLines(paste("  lscores all / sitelikes = yes scorefile =", paste(filebase, '.optimal.scores', sep = ''), ";"), paup.out)
    writeLines("\n  cleartrees;", paup.out)
    writeLines(paste("  gettrees file =", paste(filebase, '.rtrees.tre', sep = ''), ";"), paup.out)
    writeLines(paste("  lscores all / sitelikes = yes scorefile =", paste(filebase, '.rtrees.scores', sep = ''), ";"), paup.out)
    writeLines("\nend;", paup.out)
    close(paup.out)
	}

args <- commandArgs(trailingOnly=TRUE)

suppressMessages(library(rtracklayer))
setwd(args[1])

gtf <- rtracklayer::import.gff(args[2])
cat('gtf file succesfully loaded\n')
saveRDS(gtf, 'main_annotation.rds')
cat('main_annotation.rds saved\n')  

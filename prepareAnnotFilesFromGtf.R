args <- commandArgs(trailingOnly=TRUE)

if (args[1] == '--help' | args[1] == '-h') {
  cat('USAGE: Rscript prepareAnnotFilesFromGtf /path/to/workdir /path/to/gtf annot_type\n')
  cat('/path/to/workdir - absolute path to dir, where PAR-PIPE2 analysis will be performed\n')  
  cat('annot_files dir containing rds annotation files will be created there\n')
  cat('/path/to/annot.gtf - absolute path to your gtf file\n')
  cat('annot_type - one of the main|trna|rep\n')
} else {
  suppressMessages(library(rtracklayer))
  cat(paste0('Your working dir is ', args[1], '\n'))
  setwd(args[1])
  output_dir <- 'annot_files'
  
  if (!dir.exists(output_dir)){
    cat('Creating annot_files dir...\n')
    dir.create(output_dir)
  } else {
    cat('annot_files dir alreasy exists\n')
  }
  
  cat('Reading gtf file...\n')
  gtf <- rtracklayer::import.gff(args[2])
  cat('gtf file succesfully loaded\n')
  
  if (args[3] == 'main') {
    saveRDS(gtf, './annot_files/main_annotation.rds')
    cat('main_annotation.rds saved in annot_files\n')
  } else if (args[3] == 'trna') {
    saveRDS(gtf, './annot_files/trna.rds')
    cat('trna.rds saved in annot_files\n')
  } else if (args[3] == 'rep') {
    saveRDS(gtf, './annot_files/repeats.rds')
    cat('repeats.rds saved in annot_files\n')
  } else {
    cat('annot_type must be one of the main|trna|rep\n
        re-run script with the correct value\n')
  }
}
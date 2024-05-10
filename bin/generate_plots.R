suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggpubr))

args <- commandArgs(trailingOnly = TRUE)

sum_stats <- read_tsv(args[1]) %>%
  dplyr::rename('id' = '...1')

samples <- colnames(sum_stats) %>%
  .[2:(length(.)-2)] %>%
  gsub('_freq','', .) %>%
  unique()

for (i in seq_along(samples)) {
  sample_df <- sum_stats %>%
    dplyr::select(c('id',samples[i], 'analysis_step')) %>%
    dplyr::rename('counts' = samples[i]) %>%
    mutate(counts = as.numeric(counts))
  reads <- sample_df %>%
    filter(analysis_step == 'reads_processing') %>%
    mutate(frac = counts/counts[1])
  ali <- sample_df %>%
    filter(analysis_step == 'alignment') %>%
    mutate(frac = counts/counts[1])
  cl <- sample_df %>%
    filter(analysis_step == 'clusters_calling') %>%
    mutate(frac = counts/counts[1])
  
  reads_cts_plt <- ggbarplot(reads, x = 'id', y = 'counts',
                             fill = 'red', main = "reads processing") +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  reads_frac_plt <- ggbarplot(reads, x = 'id', y = 'frac',
                              fill = 'red', main = 'reads processing',
                              ylab = 'fraction') +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  ali_cts_plt <- ggbarplot(ali, x = 'id', y = 'counts',
                           fill = 'green', main = 'alignment') +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  ali_frac_plt <- ggbarplot(ali, x = 'id', y = 'frac',
                            fill = 'green', main = 'alignment',
                            ylab = 'fraction') +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  cl_cts_plt <- ggbarplot(cl, x = 'id', y = 'counts',
                          fill = 'blue', main = 'cluster calling') +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  cl_frac_plt <- ggbarplot(cl, x = 'id', y = 'frac',
                           fill = 'blue', main = 'cluster calling',
                           ylab = 'fraction') +
    rremove('xlab') +
    theme(text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  cts <- ggarrange(reads_cts_plt, ali_cts_plt, cl_cts_plt, nrow = 1, align = 'hv')
  frac <- ggarrange(reads_frac_plt, ali_frac_plt, cl_frac_plt, nrow = 1, align = 'hv')
  ggsave(cts, path = args[2], 
         filename = paste0(samples[i],'_counts_', args[3], '.pdf'),
         width = 6, height = 4, units = 'in',
         device = 'pdf')
  ggsave(frac, path = args[2], 
         filename = paste0(samples[i],'_fraction_', args[3], '.pdf'),
         width = 6, height = 4, units = 'in',
         device = 'pdf')
}


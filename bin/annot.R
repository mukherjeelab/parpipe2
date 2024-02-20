suppressMessages(library(plyranges))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
suppressMessages(library(RCAS))
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly=TRUE)
setwd(args[1])

annot_rank <- fread(args[2], header = FALSE, 
                    col.names = 'category') %>%
  mutate(rank = 1:nrow(.))

main_annot <- readRDS(args[3])
trna_annot <- readRDS(args[4])
repeats_annot <- readRDS(args[5])
merged_annot <- plyranges::bind_ranges(main_annot, trna_annot, repeats_annot)
txdbFeatures <- getTxdbFeaturesFromGRanges(main_annot) #for mRNA annotation

in_file <- fread(args[6]) %>%
  dplyr::rename(c('chromosome' = 'Chromosome', 'strand' = 'Strand', 
                  'start' = 'ClusterStart', 'end' = 'ClusterEnd')) 

in_file_gr <- makeGRangesFromDataFrame(in_file, keep.extra.columns = TRUE) #transforming tabular object into GRanges object

general_annot_tmp <-as.data.frame(queryGff(queryRegions = in_file_gr, gffData = merged_annot)) %>%
  mutate(gene_type = gsub('tRNAscan', 'tRNA', gene_type),
         annot_rank = annot_rank$rank[match(gene_type, annot_rank$category)]) %>%
  filter(is.na(annot_rank) == FALSE) %>%
  group_by(query_ClusterID) %>%
  filter(annot_rank == min(annot_rank)) %>%
  distinct(query_ClusterID, .keep_all = TRUE) %>%
  dplyr::select(-c('annot_rank')) %>%
  mutate(gene_type = factor(gene_type))

mrna_tmp <- general_annot_tmp %>%
  filter(gene_type == 'protein_coding')
mrna_gr_tmp <- in_file_gr %>%
  plyranges::filter(ClusterID %in% mrna_tmp$query_ClusterID)
mrna_ov_tmp <- lapply(txdbFeatures, function(x) GenomicRanges::findOverlaps(mrna_gr_tmp, 
                                                                            x))
cds <- as.data.frame(mrna_ov_tmp[6]) %>%
  set_colnames(c('queryHits','subjectHits')) %>%
  mutate(category = 'cds',
         rank = 1,
         transcript_id = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 9]),
         gene_name = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 10]))  
utr3 <- as.data.frame(mrna_ov_tmp[7]) %>%
  set_colnames(c('queryHits','subjectHits')) %>%
  mutate(category = 'threeUTRs',
         rank = 2,
         transcript_id = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 6]),
         gene_name = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 7]))
utr5 <- as.data.frame(mrna_ov_tmp[4]) %>%
  set_colnames(c('queryHits','subjectHits')) %>%
  mutate(category = 'fiveUTRs',
         rank = 3,
         transcript_id = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 6]),
         gene_name = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 7]))
introns <- as.data.frame(mrna_ov_tmp[5]) %>%
  set_colnames(c('queryHits','subjectHits')) %>%
  mutate(category = 'introns',
         rank = 4,
         transcript_id = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 6]),
         gene_name = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 7]))
promoters <- as.data.frame(mrna_ov_tmp[3]) %>%
  set_colnames(c('queryHits','subjectHits')) %>%
  mutate(category = 'promoters',
         rank = 5,
         transcript_id = unlist(as_tibble(txdbFeatures$promoters)[subjectHits, 7]),
         gene_name = unlist(as_tibble(txdbFeatures$promoters)[subjectHits, 8]))
df_mrna_tmp <- rbind(cds, utr3, utr5, introns, promoters) %>%
  mutate(ClusterID = mrna_gr_tmp$ClusterID[queryHits]) %>%
  group_by(ClusterID) %>%
  filter(rank == min(rank)) %>%
  distinct(ClusterID, .keep_all = TRUE) %>%
  mutate(category = factor(category, levels = c("promoters", "fiveUTRs", "cds", "introns", "threeUTRs")))

final_df <- general_annot_tmp %>%
  dplyr::select(c('queryRange', 'query_ClusterID', 'query_ClusterSequence', 'query_ReadCount',
                  'gene_id','gene_type', 'gene_name', 'query_ConversionEventCount',
                  'query_ConversionLocationCount','query_NonT2CConversionCount', 'query_logFCT2C')) %>%
  set_colnames(gsub('query|query_','', colnames(.))) %>%
  dplyr::rename(c('ClusterCoords' = 'Range')) %>%
  mutate(mRNA_region = df_mrna_tmp$category[match(ClusterID, df_mrna_tmp$ClusterID)],
         repeat_type = ifelse(gene_type == 'repeat',
                              general_annot_tmp$transcript_type[match(ClusterID, general_annot_tmp$query_ClusterID)],
                              NA)) %>%
  left_join(., in_file %>% dplyr::select(c(1,3:5,2))) %>%
  dplyr::select((ncol(.)-3):ncol(.),1:(ncol(.)-4)) %>%
  dplyr::select(c(1:11,(ncol(.)-1):ncol(.),12:(ncol(.)-2))) %>%
  arrange(desc(logFCT2C), desc(ConversionLocationCount), desc(ConversionEventCount))

fwrite(final_df, args[7], sep = '\t')

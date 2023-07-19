library(plyranges)
library(dplyr)
library(tidyverse)
library(magrittr)
library(data.table)
library(RCAS)
library(Biostrings)
library(here)

annot_rank <- fread(here('annot_files/annot_categories.txt'), header = FALSE, 
                    col.names = 'category') %>%
  mutate(rank = 1:nrow(.)) #setting annotation rank of the features
gtf <- importGtf(here('annot_files/txmax.gtf')) #main gtf
trna <- importGtf(here('annot_files/gencode.v24.tRNAs.gtf')) #tRNA gtf
repeats <- importGtf(here('annot_files/rmsk.gtf')) #repetitive elements gtf
merged_gtf <- plyranges::bind_ranges(gtf, trna, repeats) #merged gtf
txdbFeatures <- getTxdbFeaturesFromGRanges(gtf) #for mRNA annotation

in_file <- fread(here('./PARalyzer/DecapRep2_RPI2_combined.clusters')) %>%
  dplyr::rename(c('chromosome' = 'Chromosome', 'strand' = 'Strand', 
                  'start' = 'ClusterStart', 'end' = 'ClusterEnd')) %>%
  mutate(T2C_fraction = ConversionEventCount/(ConversionEventCount+NonConversionEventCount)) # reading paralyzer output, standardizing names with GRanges requirements and calculating T2C fraction

in_file_gr <- makeGRangesFromDataFrame(in_file, keep.extra.columns = TRUE) #transforming tabular object into GRanges object

general_annot_tmp <-as.data.frame(queryGff(queryRegions = in_file_gr, gffData = merged_gtf)) %>%
  mutate(gene_type = gsub('tRNAscan', 'tRNA', gene_type),
         annot_rank = annot_rank$rank[match(gene_type, annot_rank$category)]) %>%
  filter(is.na(rank) == FALSE) %>%
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

tmp <- sort(in_file_gr, ignore.strand=TRUE) %>%
  as.data.frame()
paralyzer_clusters <- rep(tmp$ClusterID, nchar(tmp$ClusterSequence))
paralyzer_strand <- rep(tmp$strand, nchar(tmp$ClusterSequence))
score_df_tmp <- fread(here("./PARalyzer/Decap_rep2.tsv"))
stopifnot(length(paralyzer_clusters) == nrow(score_df_tmp))
stopifnot(length(paralyzer_clusters) == length(paralyzer_strand))

score_df_strand_tmp <- score_df_tmp %>%
  mutate(strand = paralyzer_strand,
         ClusterID = paralyzer_clusters)

score_df_plus <- score_df_strand_tmp %>%
  filter(strand == '+') %>%
  mutate(T2C = ifelse(REF_BASE == 'T', C, 0),
         noT2C = case_when(REF_BASE == 'T' ~ A+G,
                           REF_BASE == 'A' ~ C+G+`T`,
                           REF_BASE == 'C' ~ A+G+`T`,
                           REF_BASE == 'G' ~ A+C+`T`)) %>%
  group_by(ClusterID) %>%
  mutate(total_T2C = sum(T2C, na.rm = TRUE),
         total_noT2C = sum(noT2C, na.rm = TRUE),
         total_noT2C = ifelse(total_noT2C == 0, 1, total_noT2C),
         T2C_fc = log2(total_T2C/total_noT2C),
         T2C_pos_cnt = sum(T2C>0, na.rm = TRUE))
  
score_df_minus <- score_df_strand_tmp %>%
  filter(strand == '-') %>%
  mutate(REF_BASE = case_when(REF_BASE == 'A' ~'T',
                              REF_BASE == 'C' ~'G',
                              REF_BASE == 'G' ~'C',
                              REF_BASE == 'T' ~'A')) %>%
  dplyr::rename(c('T' = 'A', 'G' = 'C', 'C' = 'G', 'A' = 'T')) %>%
  dplyr::select(1:4,8,7,6,5,9:ncol(.)) %>%
  mutate(T2C = ifelse(REF_BASE == 'T', C, 0),
         noT2C = case_when(REF_BASE == 'T' ~ A+G,
                           REF_BASE == 'A' ~ C+G+`T`,
                           REF_BASE == 'C' ~ A+G+`T`,
                           REF_BASE == 'G' ~ A+C+`T`)) %>%
  group_by(ClusterID) %>%
  mutate(total_T2C = sum(T2C, na.rm = TRUE),
         total_noT2C = sum(noT2C, na.rm = TRUE),
         total_noT2C = ifelse(total_noT2C == 0, 1, total_noT2C),
         T2C_fc = log2(total_T2C/total_noT2C),
         T2C_pos_cnt = sum(T2C>0, na.rm = TRUE))

score_df <- rbind(score_df_plus, score_df_minus) %>%
  distinct(ClusterID, .keep_all = TRUE) %>%
  arrange(desc(T2C_fc), desc(T2C_pos_cnt), desc(total_T2C)) %>%
  dplyr::select(c(16,19:ncol(.)))

final_df <- general_annot_tmp %>%
  dplyr::select(c('queryRange', 'query_ClusterID', 'query_ClusterSequence', 'query_ReadCount',
                  'gene_id','gene_type', 'gene_name')) %>%
  set_colnames(gsub('query|query_','', colnames(.))) %>%
  dplyr::rename(c('ClusterCoords' = 'Range')) %>%
  mutate(mRNA_region = df_mrna_tmp$category[match(ClusterID, df_mrna_tmp$ClusterID)],
         repeat_type = ifelse(gene_type == 'repeat',
                              general_annot_tmp$transcript_type[match(ClusterID, general_annot_tmp$query_ClusterID)],
                              NA)) %>%
  left_join(., score_df) %>%
  left_join(., in_file %>% dplyr::select(c(1,3:5,2))) %>%
  dplyr::select((ncol(.)-3):ncol(.),1:(ncol(.)-4)) %>%
  arrange(desc(T2C_fc), desc(T2C_pos_cnt), desc(total_T2C))

fwrite(final_df, here('./DecapRep2_annot_score.tsv'), sep = '\t')

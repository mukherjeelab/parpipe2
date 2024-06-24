suppressWarnings(suppressMessages(library(plyranges)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(RCAS)))
#suppressWarnings(suppressMessages(library(Biostrings)))
#suppressWarnings(suppressMessages(library(metagene2)))
suppressWarnings(suppressMessages(library(ggpubr)))

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

groups_file <- fread(args[6]) %>%
  dplyr::rename(c('chromosome' = 'Chromosome', 'strand' = 'Strand', 
                  'start' = 'GroupStart', 'end' = 'GroupEnd')) 
clust_file <- fread(args[7]) %>%
  dplyr::rename(c('chromosome' = 'Chromosome', 'strand' = 'Strand', 
                  'start' = 'ClusterStart', 'end' = 'ClusterEnd')) 

cat("Annotating groups \n")

if (nrow(groups_file) > 0) {
  groups_file_gr <- makeGRangesFromDataFrame(groups_file, keep.extra.columns = TRUE)
  groups_general_annot_tmp_gr <-queryGff(queryRegions = groups_file_gr, gffData = merged_annot)
  groups_general_annot_tmp <- as.data.frame(groups_general_annot_tmp_gr) %>%
    mutate(gene_type = gsub('tRNAscan', 'tRNA', gene_type),
           annot_rank = annot_rank$rank[match(gene_type, annot_rank$category)]) %>%
    filter(is.na(annot_rank) == FALSE) %>%
    group_by(query_GroupID) %>%
    filter(annot_rank == min(annot_rank)) %>%
    distinct(query_GroupID, .keep_all = TRUE) %>%
    dplyr::select(-c('annot_rank')) %>%
    mutate(gene_type = factor(gene_type))
  groups_unannotated_regs <- groups_file_gr %>%
    .[!.$GroupID %in% groups_general_annot_tmp_gr$query_GroupID] %>%
    as.data.frame() %>%
    mutate(GroupCoords = paste0(seqnames, ':', start, '-', end, ':', strand),
           gene_id = paste0('unannotated_', row.names(.)),
           gene_type = 'unannotated_region',
           gene_name = gene_id) %>%
    dplyr::select(c("GroupCoords","GroupID":"T_del","gene_id":"gene_name"))
  groups_mrna_tmp <- groups_general_annot_tmp %>%
    filter(gene_type == 'protein_coding')
  groups_mrna_gr_tmp <- groups_file_gr %>%
    plyranges::filter(GroupID %in% groups_mrna_tmp$query_GroupID)
  groups_mrna_ov_tmp <- lapply(txdbFeatures, function(x) GenomicRanges::findOverlaps(groups_mrna_gr_tmp, 
                                                                              x))
  groups_cds <- as.data.frame(groups_mrna_ov_tmp[6]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'cds',
           rank = 1,
           transcript_id = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 9]),
           gene_name = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 10]))  
  groups_utr3 <- as.data.frame(groups_mrna_ov_tmp[7]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'threeUTRs',
           rank = 2,
           transcript_id = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 7]))
  groups_utr5 <- as.data.frame(groups_mrna_ov_tmp[4]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'fiveUTRs',
           rank = 3,
           transcript_id = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 7]))
  groups_introns <- as.data.frame(groups_mrna_ov_tmp[5]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'introns',
           rank = 4,
           transcript_id = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 7]))
  groups_df_mrna_tmp <- rbind(groups_cds, groups_utr3, groups_utr5, groups_introns) %>%
    mutate(GroupID = groups_mrna_gr_tmp$GroupID[queryHits]) %>%
    group_by(GroupID) %>%
    filter(rank == min(rank)) %>%
    distinct(GroupID, .keep_all = TRUE) %>%
    mutate(category = factor(category, levels = c("fiveUTRs", "cds", "introns", "threeUTRs")))
  
  groups_final_df_full <- groups_general_annot_tmp %>%
    dplyr::select(c('queryRange':ncol(.),'gene_id','gene_type', 'gene_name')) %>%
    set_colnames(gsub('query|query_','', colnames(.))) %>%
    dplyr::rename('GroupCoords' = 'Range') %>%
    rbind(., groups_unannotated_regs) %>%
    mutate(mRNA_region = groups_df_mrna_tmp$category[match(GroupID, groups_df_mrna_tmp$GroupID)],
           repeat_type = ifelse(gene_type == 'repeat',
                                groups_general_annot_tmp$transcript_type[match(GroupID, groups_general_annot_tmp$query_GroupID)],
                                NA)) %>%
    left_join(., groups_file %>% dplyr::select(c(1,3:5,2))) %>%
    dplyr::select((ncol(.)-3):ncol(.),1:(ncol(.)-4)) %>%
    dplyr::select(c(1:7,(ncol(.)-4):ncol(.),8:(ncol(.)-5))) %>%
    arrange(desc(ReadCount))
  fwrite(groups_final_df_full, gsub('.tsv', '_annotated.tsv', args[6]), sep = '\t')
  groups_final_df_full %>%
    dplyr::select(-c('A_match':'T_del')) %>%
    fwrite(gsub('_conv_stats.tsv', '_annotated.tsv', gsub('stats/', 'annot/', args[6])), sep = '\t')
  
} else {
  cat("No groups to annotate \n")
  groups_final_df_full <- data.frame(matrix(nrow = 0, ncol = 33)) %>%
    set_colnames(c("chromosome","start","end","strand","GroupCoords","GroupID",
                   "GroupSequence","gene_id","gene_type","gene_name","mRNA_region","repeat_type",
                   "ReadCount","A_match","A_C","A_G","A_T","A_del",
                   "C_match","C_A","C_G","C_T","C_del","G_match",
                   "G_A","G_C","G_T","G_del","T_match","T_A",
                   "T_C","T_G","T_del"))
  fwrite(groups_final_df_full, gsub('.tsv', '_annotated.tsv', args[6]), sep = '\t')
  groups_final_df_full %>%
    dplyr::select(-c('A_match':'T_del')) %>%
    fwrite(gsub('_conv_stats.tsv', '_annotated.tsv', gsub('stats/', 'annot/', args[6])), sep = '\t')
}
    
cat("Annotating clusters \n")

if (nrow(clust_file) > 0) {
  clust_file_gr <- makeGRangesFromDataFrame(clust_file, keep.extra.columns = TRUE)
  clust_general_annot_tmp_gr <-queryGff(queryRegions = clust_file_gr, gffData = merged_annot)
  clust_general_annot_tmp <- as.data.frame(clust_general_annot_tmp_gr) %>%
    mutate(gene_type = gsub('tRNAscan', 'tRNA', gene_type),
           annot_rank = annot_rank$rank[match(gene_type, annot_rank$category)]) %>%
    filter(is.na(annot_rank) == FALSE) %>%
    group_by(query_ClusterID) %>%
    filter(annot_rank == min(annot_rank)) %>%
    distinct(query_ClusterID, .keep_all = TRUE) %>%
    dplyr::select(-c('annot_rank')) %>%
    mutate(gene_type = factor(gene_type))
  clust_unannotated_regs <- clust_file_gr %>%
    .[!.$ClusterID %in% clust_general_annot_tmp_gr$query_ClusterID] %>%
    as.data.frame() %>%
    mutate(ClusterCoords = paste0(seqnames, ':', start, '-', end, ':', strand),
           gene_id = paste0('unannotated_', row.names(.)),
           gene_type = 'unannotated_region',
           gene_name = gene_id) %>%
    dplyr::select(c(16,6:8,17:19,12,11,14,15))
  clust_mrna_tmp <- clust_general_annot_tmp %>%
    filter(gene_type == 'protein_coding')
  clust_mrna_gr_tmp <- clust_file_gr %>%
    plyranges::filter(ClusterID %in% clust_mrna_tmp$query_ClusterID)
  clust_mrna_ov_tmp <- lapply(txdbFeatures, function(x) GenomicRanges::findOverlaps(clust_mrna_gr_tmp, 
                                                                              x))
  clust_cds <- as.data.frame(clust_mrna_ov_tmp[6]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'cds',
           rank = 1,
           transcript_id = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 9]),
           gene_name = unlist(as_tibble(txdbFeatures$cds)[subjectHits, 10]))  
  clust_utr3 <- as.data.frame(clust_mrna_ov_tmp[7]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'threeUTRs',
           rank = 2,
           transcript_id = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$threeUTRs)[subjectHits, 7]))
  clust_utr5 <- as.data.frame(clust_mrna_ov_tmp[4]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'fiveUTRs',
           rank = 3,
           transcript_id = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$fiveUTRs)[subjectHits, 7]))
  clust_introns <- as.data.frame(clust_mrna_ov_tmp[5]) %>%
    set_colnames(c('queryHits','subjectHits')) %>%
    mutate(category = 'introns',
           rank = 4,
           transcript_id = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 6]),
           gene_name = unlist(as_tibble(txdbFeatures$introns)[subjectHits, 7]))
  clust_df_mrna_tmp <- rbind(clust_cds, clust_utr3, clust_utr5, clust_introns) %>%
    mutate(ClusterID = clust_mrna_gr_tmp$ClusterID[queryHits]) %>%
    group_by(ClusterID) %>%
    filter(rank == min(rank)) %>%
    distinct(ClusterID, .keep_all = TRUE) %>%
    mutate(category = factor(category, levels = c("fiveUTRs", "cds", "introns", "threeUTRs")))
  
  clust_final_df_full <- clust_general_annot_tmp %>%
    dplyr::select(c('queryRange', 'query_ClusterID', 'query_ClusterSequence', 'query_ReadCount',
                    'gene_id','gene_type', 'gene_name', 'query_ConversionEventCount',
                    'query_ConversionLocationCount','query_NonT2CConversionCount', 'query_conversion_specificity')) %>%
    set_colnames(gsub('query|query_','', colnames(.))) %>%
    dplyr::rename('ClusterCoords' = 'Range') %>%
    rbind(., clust_unannotated_regs) %>%
    mutate(mRNA_region = clust_df_mrna_tmp$category[match(ClusterID, clust_df_mrna_tmp$ClusterID)],
           repeat_type = ifelse(gene_type == 'repeat',
                                clust_general_annot_tmp$transcript_type[match(ClusterID, clust_general_annot_tmp$query_ClusterID)],
                                NA)) %>%
    left_join(., clust_file %>% dplyr::select(c(1,3:5,2))) %>%
    dplyr::select((ncol(.)-3):ncol(.),1:(ncol(.)-4)) %>%
    dplyr::select(c(1:11,(ncol(.)-1):ncol(.),12:(ncol(.)-2))) %>%
    arrange(desc(conversion_specificity), desc(ConversionLocationCount), desc(ConversionEventCount)) 
  fwrite(clust_final_df_full, gsub('.tsv', '_annotated.tsv', args[7]), sep = '\t')
  clust_final_df_full %>%
    dplyr::select(-c("ConversionEventCount":"NonT2CConversionCount")) %>%
    fwrite(gsub('_conv_stats.tsv', '_annotated.tsv', gsub('stats/', 'annot/', args[7])), sep = '\t')
  clust_final_df_full %>%
    filter(conversion_specificity > .6) %>%
    dplyr::select(c(1:3,6,17,4)) %>%
    fwrite(gsub('_conv_stats.tsv', '_filtered.bed', gsub('stats/', 'genome_viewer_files/', args[7])), 
           sep = '\t', col.names = FALSE)
  clust_final_df_full %>%
    dplyr::select(c(1:3,6,4,17)) %>%
    mutate(score = ifelse(conversion_specificity > .6, 60, 0)) %>%
    dplyr::select(c(1:4,7,5)) %>%
    fwrite(gsub('_conv_stats.tsv', '_unfiltered.bed', gsub('stats/', 'genome_viewer_files/', args[7])), 
           sep = '\t', col.names = FALSE)
} else {
  cat("No clusters to annotate \n")
  clust_final_df_full <- data.frame(matrix(nrow = 0, ncol = 17)) %>%
    set_colnames(c("chromosome","start","end","strand","ClusterCoords","ClusterID","ClusterSequence","ReadCount",
                   "gene_id","gene_type","gene_name","mRNA_region","repeat_type","ConversionEventCount",
                   "ConversionLocationCount","NonT2CConversionCount","conversion_specificity"))
  fwrite(clust_final_df_full, gsub('.tsv', '_annotated.tsv', args[7]), sep = '\t')
  clust_final_df_full %>%
    dplyr::select(-c("ConversionEventCount":"NonT2CConversionCount")) %>%
    fwrite(gsub('_conv_stats.tsv', '_annotated.tsv', gsub('stats/', 'annot/', args[7])), sep = '\t')
  clust_final_df_full %>%
    dplyr::select(c(1:3,6,17,4)) %>%
    fwrite(gsub('_conv_stats.tsv', '_filtered.bed', gsub('stats/', 'genome_viewer_files/', args[7])), 
           sep = '\t', col.names = FALSE)
  clust_final_df_full %>%
    dplyr::select(c(1:3,6,17,4)) %>%
    fwrite(gsub('_conv_stats.tsv', '_unfiltered.bed', gsub('stats/', 'genome_viewer_files/', args[7])), 
           sep = '\t', col.names = FALSE)
}
import os
import pysam
import pysamstats
import math
import argparse
import concurrent.futures

cwd = os.getcwd()

def process_group_line(line, bamfile, genome_ref):
    bam = pysam.AlignmentFile(bamfile, "rb")
    group_line = line.strip().split(',')[:7]
    if line.startswith('Chromosome'):
        group_line[6] = 'ReadCount'
        group_line.extend([
            'A_match', 'A_C', 'A_G', 'A_T', 'A_del',
            'C_match', 'C_A', 'C_G', 'C_T', 'C_del',
            'G_match', 'G_A', 'G_C', 'G_T', 'G_del',
            'T_match', 'T_A', 'T_C', 'T_G', 'T_del'
        ])
        result = '\t'.join(group_line) + '\n'
    else:
        group_line[6] = int(group_line[6])
        group_line.extend([0] * 20)
        for rec in pysamstats.stat_variation_strand(
                bam, chrom=group_line[0], start=int(group_line[2]), end=int(group_line[3]),
                pad=True, truncate=True, fafile=genome_ref, one_based=False):
            if group_line[1] == '+':
                if rec['ref'] == 'A':
                    group_line[7] += int(rec['A_fwd'])
                    group_line[8] += int(rec['C_fwd'])
                    group_line[9] += int(rec['G_fwd'])
                    group_line[10] += int(rec['T_fwd'])
                    group_line[11] += int(rec['deletions_fwd'])
                elif rec['ref'] == 'C':
                    group_line[12] += int(rec['C_fwd'])
                    group_line[13] += int(rec['A_fwd'])
                    group_line[14] += int(rec['G_fwd'])
                    group_line[15] += int(rec['T_fwd'])
                    group_line[16] += int(rec['deletions_fwd'])
                elif rec['ref'] == 'G':
                    group_line[17] += int(rec['G_fwd'])
                    group_line[18] += int(rec['A_fwd'])
                    group_line[19] += int(rec['C_fwd'])
                    group_line[20] += int(rec['T_fwd'])
                    group_line[21] += int(rec['deletions_fwd'])
                elif rec['ref'] == 'T':
                    group_line[22] += int(rec['T_fwd'])
                    group_line[23] += int(rec['A_fwd'])
                    group_line[24] += int(rec['C_fwd'])
                    group_line[25] += int(rec['G_fwd'])
                    group_line[26] += int(rec['deletions_fwd'])
            else:
                if rec['ref'] == 'T':
                    group_line[7] += int(rec['T_rev'])
                    group_line[8] += int(rec['G_rev'])
                    group_line[9] += int(rec['C_rev'])
                    group_line[10] += int(rec['A_rev'])
                    group_line[11] += int(rec['deletions_rev'])
                elif rec['ref'] == 'G':
                    group_line[12] += int(rec['G_rev'])
                    group_line[13] += int(rec['T_rev'])
                    group_line[14] += int(rec['C_rev'])
                    group_line[15] += int(rec['A_rev'])
                    group_line[16] += int(rec['deletions_rev'])
                elif rec['ref'] == 'C':
                    group_line[17] += int(rec['C_rev'])
                    group_line[18] += int(rec['T_rev'])
                    group_line[19] += int(rec['G_rev'])
                    group_line[20] += int(rec['A_rev'])
                    group_line[21] += int(rec['deletions_rev'])
                elif rec['ref'] == 'A':
                    group_line[22] += int(rec['A_rev'])
                    group_line[23] += int(rec['T_rev'])
                    group_line[24] += int(rec['G_rev'])
                    group_line[25] += int(rec['C_rev'])
                    group_line[26] += int(rec['deletions_rev'])
        result = '\t'.join(map(str, group_line)) + '\n'
    bam.close()
    return result

def process_clust_line(line, bamfile, genome_ref, conversion):
    bam = pysam.AlignmentFile(bamfile, "rb")
    clust_line = line.strip().split(',')
    if line.startswith('Chromosome'):
        clust_line[11] = 'NonConversionEventCount'
        if conversion == "T>C":
            clust_line.extend(['NonT2CConversionCount', 'conversion_specificity'])
        elif conversion == "G>A":
            clust_line.extend(['NonG2AConversionCount', 'conversion_specificity'])
        #print(f"Extended header: {clust_line}")  # Debug statement
        result = '\t'.join(clust_line) + '\n'
    else:
        clust_line[11] = int(clust_line[11])
        clust_line.append(0)
        for rec in pysamstats.stat_variation_strand(
                bam, chrom=clust_line[0], start=int(clust_line[2]), end=int(clust_line[3]),
                pad=True, truncate=True, fafile=genome_ref, one_based=False):
            if conversion == "T>C":
                if clust_line[1] == '+':
                    if rec['ref'] == 'T':
                        noTC = int(rec['mismatches_fwd']) - int(rec['C_fwd'])
                        clust_line[12] += noTC
                    else:
                        clust_line[12] += int(rec['mismatches_fwd'])
                else:
                    if rec['ref'] == 'A':
                        noTC = int(rec['mismatches_rev']) - int(rec['G_rev'])
                        clust_line[12] += noTC
                    else:
                        clust_line[12] += int(rec['mismatches_rev'])
            elif conversion == "G>A":
                if clust_line[1] == '+':
                    if rec['ref'] == 'G':
                        noGA = int(rec['mismatches_fwd']) - int(rec['A_fwd'])
                        clust_line[12] += noGA
                    else:
                        clust_line[12] += int(rec['mismatches_fwd'])
                else:
                    if rec['ref'] == 'C':
                        noGA = int(rec['mismatches_rev']) - int(rec['T_rev'])
                        clust_line[12] += noGA
                    else:
                        clust_line[12] += int(rec['mismatches_rev'])
        if clust_line[12] == 0:
            clust_line[12] += 1
        clust_line.append(math.log2(int(clust_line[10]) / clust_line[12]))
        result = '\t'.join(map(str, clust_line)) + '\n'
    bam.close()
    return result

def run_clusters_stats(bamfile, groups, clusters, genome, threads, conversion):
    max_threads = threads
    genome_ref = genome

    groupfile = groups
    output_stats_group = os.path.join(cwd, 'stats', f'{os.path.basename(groupfile).replace(".groups", "_groups_conv_stats.tsv")}')

    clustfile = clusters
    output_stats_clust = os.path.join(cwd, 'stats', f'{os.path.basename(clustfile).replace(".clusters", "_clusters_conv_stats.tsv")}')

    print(f'Adding groups conversion stats to {groupfile.replace(".groups", "")}')
    with open(groupfile, 'r') as f, open(output_stats_group, 'w') as g:
        lines = f.readlines()
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
            results = list(executor.map(lambda line: process_group_line(line, bamfile, genome_ref), lines))
        for result in results:
            g.write(result)
    print('Conversion stats added to groups')

    print(f'Adding clusters conversion stats to {clustfile.replace(".clusters", "")}')
    with open(clustfile, 'r') as f, open(output_stats_clust, 'w') as g:
        lines = f.readlines()
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
            results = list(executor.map(lambda line: process_clust_line(line, bamfile, genome_ref, conversion), lines))
        for result in results:
            g.write(result)
    print('Conversion stats added to clusters')

    return output_stats_group, output_stats_clust

def main(args):
    run_clusters_stats(args.bamfile, args.groups, args.clusters, args.genome, args.threads, args.conv)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add conversion stats to PARalyzer groups and clusters')
    parser.add_argument('--bamfile', type=str, help='PARalyzer utilized BAM', required=True)
    parser.add_argument('--groups', type=str, help='PARalyzer groups file', required=True)
    parser.add_argument('--clusters', type=str, help='PARalyzer clusters file', required=True)
    parser.add_argument('--genome', type=str, help='Genome fasta file', required=True)
    parser.add_argument('--threads', type=int, help='Number of threads', required=True)
    parser.add_argument('--conv', type=str, help='Conversion, possible values are T>C for s4U or G>A for s6G', required=True)
    args = parser.parse_args()
    main(args)


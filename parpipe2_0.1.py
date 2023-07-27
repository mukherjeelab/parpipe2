import os
import subprocess
import sys
import shutil
import argparse
import yaml
import pysam
import pysamstats
import math
from collections import OrderedDict

cwd = os.getcwd()

print("Loading config file..." )
    
config_file = sys.argv[2]

with open(config_file) as f_cfg:
    cfg = yaml.load(f_cfg, Loader=yaml.FullLoader)


print("Config file succesfully loaded")

threads = cfg['threads']
    
def trim_adapters_round1(reads, samplename):
    if not os.path.exists('reads'):
        os.mkdir('reads')

    outdir = os.path.join(cwd, 'reads')
    
    samplename = samplename.split(",") 
    reads = reads.split(",")
    
    for idx, sample in enumerate(samplename):
        
        read = reads[idx]
        output = os.path.join(cwd, 'reads', '{0}_trimmed1.fq.gz'.format(sample))
        logfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round1.log'.format(sample))
        errfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round1.err'.format(sample))
        command = ['cutadapt'] + cfg['cutadapt_params1'].split(",") +  ["-o", output, read]
        #print(command)
        
        print("Running cutadapt round 1 for {0}...".format(sample))
        #print(type(command))
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(command, stdout=f, stderr = g)
            outcode = subprocess.run(command, stdout=f, stderr=g).returncode

        if outcode == 0:     
            print("cutadapt round 1 finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during cutadapt round 1 :(\ncheck {0} for details".format(errfile)) 
            

def trim_adapters_round2(samplename):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    
    for sample in samplename:
        infile = os.path.join(cwd, 'reads', '{0}_trimmed1.fq.gz'.format(sample))
        output = os.path.join(cwd, 'reads', '{0}_trimmed_final.fq.gz'.format(sample))
        logfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round2.log'.format(sample))
        errfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round2.err'.format(sample))
        
        if cfg['cutadapt_params2'] == '' or cfg['cutadapt_params2'] == None:
            os.rename(infile, output)
            print("Skipping cutadapt round 2 for {0}...".format(sample))
            
        else:
            command = ['cutadapt'] + cfg['cutadapt_params2'].split(',') +  ['-o', output, infile]
            print("Running cutadapt round 2 for {0}...".format(sample))
            
            with open(logfile, 'w') as f, open(errfile, 'w') as g:
                subprocess.run(command, stdout=f, stderr = g)
                outcode = subprocess.run(command, stdout=f, stderr=g).returncode
             
            if outcode == 0:     
                print("cutadapt round 2 finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
            else:
                print("an ERROR occured during cutadapt round 2 :(\ncheck {0} for details".format(errfile)) 
                

def collapse_reads(samplename, umi = None):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    
    for sample in samplename:
        infile = os.path.join(cwd, 'reads', '{0}_trimmed_final.fq.gz'.format(sample))
        output = os.path.join(cwd, 'reads', '{0}_trimmed_collapsed.fa'.format(sample))
        stdoutfile =  os.path.join(cwd, 'reads', '{0}_trimmed.fa'.format(sample))
        errfile =  os.path.join(cwd, 'reads', '{0}_collapse_reads.err'.format(sample))
        
        fq_to_fa_command = ['seqtk', 'seq', '-a', infile]
        
        print("Converting FASTQ.GZ to FASTA {0}...".format(sample))
        
        with open(stdoutfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(fq_to_fa_command, stdout=f, stderr = g)
            outcode = subprocess.run(fq_to_fa_command, stdout=f, stderr=g).returncode
             
        if outcode == 0:     
            print("FASTQ.GZ to FASTA conversion finished succesfully for {0}\nFASTA file {1} created".format(sample,stdoutfile))
        else:
            print("an ERROR occured during FASTQ.GZ to FASTA conversion :(\ncheck {0} for details".format(errfile)) 
        
        print("Collapsing reads in {0}...".format(sample))
                
        if not umi:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "sed", "'s/>/>Read/g'", ">", output]
        elif umi > 0:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "fastx_trimmer", "-f", str(umi+1), "|", "sed", "'s/>/>Read/g'", ">", output]
        elif umi < 0:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "fastx_trimmer", "-t", str(-umi), "|", "sed", "'s/>/>Read/g'", ">", output]
        #print(collapse_cmd)
        os.system(" ".join(collapse_cmd))

        print("Reads collapsed succesfully for {0}\nCollapsed FASTA file {0} created".format(sample))
                 
#star_rm_reps
def star_rm_reps(samplename):
    if not os.path.exists('ali'):
        os.mkdir('ali')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    rep_idx = cfg['rep_idx']
    
    
    
    for sample in samplename:
    
        if cfg['star_rem_reps_params'] == '':
            infile = os.path.join(cwd, 'reads', '{0}_trimmed_collapsed.fa'.format(sample))
            output = os.path.join(cwd, 'reads', '{0}_rmreps.fa'.format(sample))
            shutil.copy2(infile, output)    
               
        else:        
            infile = os.path.join(cwd, 'reads', '{0}_trimmed_collapsed.fa'.format(sample))
            output = os.path.join(cwd, 'reads', '{0}_rmreps.fa'.format(sample))
            star_prefix = os.path.join(cwd, 'ali', '{0}_rmreps'.format(sample))
            logfile =  os.path.join(cwd, 'ali', '{0}_rmreps.log'.format(sample))
            errfile =  os.path.join(cwd, 'ali', '{0}_rmreps.err'.format(sample))
        
            command = ['STAR'] + cfg['star_rem_reps_params'].split(',') + ['--runThreadN', str(threads), '--genomeDir', rep_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile]
        
            print("Running STAR to remove repetitive sequences for {0}...".format(sample))
        
            with open(logfile, 'w') as f, open(errfile, 'w') as g:
                subprocess.run(command, stdout=f, stderr = g)
                outcode = subprocess.run(command, stdout=f, stderr=g).returncode
             
            if outcode == 0:     
                print("STAR finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
            else:
                print("an ERROR occured during STAR run :(\ncheck {0} for details".format(errfile)) 
            
            os.rename('{0}Unmapped.out.mate1'.format(star_prefix), output)
            os.rename('{0}Aligned.out.bam'.format(star_prefix), os.path.join(cwd, 'ali', '{0}.reps.bam'.format(sample)))
            
#star_genome_map 
def star_genome_map(samplename):

    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    genome_idx = cfg['genome_idx']
    
    for sample in samplename:
        infile = os.path.join(cwd, 'reads', '{0}_rmreps.fa'.format(sample))
        output = os.path.join(cwd, 'ali', '{0}.aligned.mapped.sorted.sam'.format(sample))
        star_prefix = os.path.join(cwd, 'ali', '{0}'.format(sample))
        tmp_bam = os.path.join(cwd, 'ali', '{0}.aligned.bam'.format(sample))
        tmp_sam = os.path.join(cwd, 'ali', '{0}.aligned.mapped.sam'.format(sample))
        logfile =  os.path.join(cwd, 'ali', '{0}_genome_map.log'.format(sample))
        errfile =  os.path.join(cwd, 'ali', '{0}_genome_map.err'.format(sample))
        
        star_command = ['STAR'] + cfg['star_map_params'].split(',') + ['--runThreadN', str(threads), '--genomeDir', genome_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile]
        
        print("Running STAR to align to genome for {0}...".format(sample))
        
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(star_command, stdout=f, stderr = g)
            outcode = subprocess.run(star_command, stdout=f, stderr=g).returncode
             
        if outcode == 0:     
            print("STAR finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during STAR run :(\ncheck {0} for details".format(errfile)) 
            
        os.rename('{0}Aligned.out.bam'.format(star_prefix), tmp_bam)
        bam_to_mapped_sam_command = ['samtools', 'view', '-F', '4', '-h', '-o', tmp_sam, tmp_bam]
        sorting_sam_command = ['samtools', 'sort', '-o', output, tmp_sam]
        
        print("Filtering mapped reads and converting BAM to SAM for {0}...".format(sample))
        subprocess.run(bam_to_mapped_sam_command)
        print("Filtering and conversion finished succesfully for {0}...\nSorting SAM file...".format(sample))
        subprocess.run(sorting_sam_command)
        print('Sorted SAM file created for {0}'.format(sample))
        
        if cfg['star_rem_reps_params'] == '':
            os.remove(infile) #rm_rep step was skipped and this is a copy of collapsed fasta file
        
#config_and_run_PARalyzer  
def config_and_run_PARalyzer(samplename):

    if not os.path.exists('PARalyzer'):
        os.mkdir('PARalyzer')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'PARalyzer')
    paralyzer_mem = cfg['paralyzer_memory']
    
    for sample in samplename:
        sam_file = os.path.join(cwd,'ali', '{0}.aligned.mapped.sorted.sam'.format(sample))
        logfile =  os.path.join(cwd,outdir, '{0}_PARalyzer.log'.format(sample))
        errfile =  os.path.join(cwd,outdir, '{0}__PARalyzer.err'.format(sample))
        config_out = os.path.join(cwd,outdir, '{0}.ini'.format(sample))
        
        with open(config_out, 'w') as ini_file:
            ini_file.write('#\nBANDWIDTH={0}\nCONVERSION={1}\nMINIMUM_READ_COUNT_PER_GROUP={2}\nMINIMUM_READ_COUNT_PER_CLUSTER={3}\nMINIMUM_READ_COUNT_FOR_KDE={4}\nMINIMUM_CLUSTER_SIZE={5}\nMINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER={6}\nMINIMUM_CONVERSION_COUNT_FOR_CLUSTER={7}\nMINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION={8}\nMINIMUM_READ_LENGTH={9}\nMAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES={10}\n\n'.format(cfg['BANDWIDTH'],cfg['CONVERSION'],cfg['MINIMUM_READ_COUNT_PER_GROUP'],            cfg['MINIMUM_READ_COUNT_PER_CLUSTER'],cfg['MINIMUM_READ_COUNT_FOR_KDE'],cfg['MINIMUM_CLUSTER_SIZE'],cfg['MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER'],            cfg['MINIMUM_CONVERSION_COUNT_FOR_CLUSTER'],cfg['MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION'],cfg['MINIMUM_READ_LENGTH'],cfg['MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES']))
            if cfg['signal_extension'] == 'EXTEND_BY_READ':
                ini_file.write('EXTEND_BY_READ\n\n')
            elif cfg['signal_extension'] == 'HAFNER_APPROACH':
                ini_file.write('HAFNER_APPROACH\n\n')
            elif cfg['signal_extension'] == 'ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL':
                ini_file.write('ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL={0}\n\n'.format(cfg['ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL']))
            ini_file.write('GENOME_2BIT_FILE={0}\nSAM_FILE={1}\n\nOUTPUT_DISTRIBUTIONS_FILE={2}.distribution\nOUTPUT_GROUPS_FILE={2}.groups\nOUTPUT_CLUSTERS_FILE={2}.clusters\nOUTPUT_READS_FILE={2}_PARalyzer_Utilized.sam\n\n'.format(cfg['GENOME_2BIT_FILE'],sam_file,os.path.join(cwd,outdir,sample)))
            
        paralyzer_cmd = ['PARalyzer', paralyzer_mem, config_out]   
        
        print("Running PARalyzer for {0}...".format(sample)) 
        
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(paralyzer_cmd, stdout=f, stderr = g)
            outcode = subprocess.run(paralyzer_cmd, stdout=f, stderr=g).returncode
            
        if outcode == 0:     
            print("PARalyzer finished succesfully for {0}".format(sample))
        else:
            print("an ERROR occured during PARalyzer run :(\ncheck {0} for details".format(errfile)) 
        

#run_clusters_stats
def run_clusters_stats(samplename):

    samplename = samplename.split(",") 
    
    for sample in samplename:
        clustfile = os.path.join(cwd,'PARalyzer', '{0}.clusters'.format(sample))
        samfile = os.path.join(cwd,'PARalyzer', '{0}_PARalyzer_Utilized.sam'.format(sample))
        bamfile = os.path.join(cwd,'PARalyzer', '{0}_PARalyzer_Utilized.sorted.bam'.format(sample))
        output = os.path.join(cwd,'PARalyzer', '{0}_clusters_conv_stats.tsv'.format(sample))
        genome_ref = cfg['genome_fasta']
 
        samtobam = ['samtools', 'sort', '-o', bamfile, samfile]
        bam_idx = ['samtools', 'index', bamfile]
        
        print("Converting PARalyzer output *.sam to sorted *.bam for {0}...".format(sample))
        subprocess.run(samtobam)
        print("PARalyzer output converted to *.bam for {0}".format(sample))
        print("Indexing {0} bam file...".format(sample))
        subprocess.run(bam_idx)
        print("{0} bam file indexed".format(sample))
        print("Adding stats to PARalyzer output for {0}...".format(sample))
        
        bam = pysam.AlignmentFile(bamfile, "rb")
        
        with open(clustfile, 'r') as f, open(output, 'w') as g:
           
           for line in f:
                clust_line = line.split(',')
                if line.startswith('Chromosome'):
                    clust_line[11] = 'NonConversionEventCount' #remove new line character
                    clust_line.extend(('NonT2CConversionCount', 'logFCT2C'))
                    g.write('\t'.join(clust_line) +'\n')
                else:
                    clust_line[11] = int(clust_line[11]) #remove new line character
                    clust_line.append(0)
                    for rec in pysamstats.stat_variation_strand(bam, chrom = clust_line[0], start=int(clust_line[2]), end = int(clust_line[3]), pad=True, truncate=True, fafile = genome_ref, one_based=False):
       
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
                    if clust_line[12] == 0:
                        clust_line[12] += 1
                    clust_line.append(math.log2(int(clust_line[10])/clust_line[12]))
                    g.write('\t'.join(str(x) for x in clust_line) +'\n')

        bam.close()
                         
        print("Stats added to PARalyzer output for {0}...".format(sample))
        
          
#run_annotR     
def run_annot(samplename):

    if not os.path.exists('output'):
        os.mkdir('output')
        
    samplename = samplename.split(",") 
    parpipe2_dir = cfg['parpipe2_dir']
    annot_rank_file = cfg['annot_rank']
    annotR = os.path.join(parpipe2_dir, 'annot.R')
    
    for sample in samplename:
        infile_clust = os.path.join(cwd,'PARalyzer', '{0}_clusters_conv_stats.tsv'.format(sample))
        logfile =  os.path.join(cwd,'output', '{0}_annotation.log'.format(sample))
        errfile =  os.path.join(cwd,'output', '{0}_annotation.err'.format(sample))
        
        command = ['Rscript', annotR, cwd, annot_rank_file, infile_clust, sample]
        
        print("Running annotation for {0}...".format(sample)) 
        
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(command, stdout=f, stderr = g)
            outcode = subprocess.run(command, stdout=f, stderr=g).returncode
            
        if outcode == 0:     
            print("Annotation finished succesfully for {0}".format(sample))
        else:
            print("an ERROR occured during annotation :(\ncheck {0} for details".format(errfile))
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Analyze PAR-CLIP')
    parser.add_argument('--config_file', type = str, help = 'Yaml file with cutadapt and STAR parameters', required = True)
    parser.add_argument('--reads', type = str, help = 'Reads containing RNA sequence (dependnent on library protocol may be either read1 or read2, if > 1, comma separated)). Gzipped fastq.', required = True)
    parser.add_argument('--samplename', type = str, help = 'Sample name (if > 1, comma separated). Will be appended to output files.', required = True)
    parser.add_argument('--umi', type = int, help = "UMI location, if 5' UMI - positive number, 3' UMI negative number,  default None", required = False)
    args = parser.parse_args()   

    config_file = os.path.abspath(args.config_file)
    reads = os.path.abspath(args.reads)    
    samplename = args.samplename
    if args.umi:
        umi = args.umi
    elif not args.umi:
        umi = None

#read_configfile(config_file)
trim_adapters_round1(reads, samplename)
trim_adapters_round2(samplename)
collapse_reads(samplename, umi)
star_rm_reps(samplename)
star_genome_map(samplename)
config_and_run_PARalyzer(samplename)
run_clusters_stats(samplename)
run_annot(samplename)

#TODO   
#consider using pysam instead samtools
#replace sys.argv with argparse for config file
#deleting temp or useless files after each step
#map to transcriptome (manually curated and annotated) instead of genome (salmon + postmaster + curated sequences from Tom)
#write manual

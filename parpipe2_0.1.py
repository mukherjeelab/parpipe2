import os
import subprocess
import sys
import shutil
import argparse
import yaml
from collections import OrderedDict

cwd = os.getcwd()

print("Loading config file..." )
    
config_file = os.path.join(cwd, sys.argv[2])

with open(config_file) as f_cfg:
    cfg = yaml.load(f_cfg, Loader=yaml.FullLoader)
f_cfg.close()

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
        output = outdir + '/' + sample + '_trimmed1.fq.gz'
        logfile =  str(outdir) + '/' + str(sample) + '_cutadapt_round1.log'
        errfile =  str(outdir) + '/' + str(sample) + '_cutadapt_round1.err'
        command = ['cutadapt'] + cfg['cutadapt_params1'].split(",") +  ["-o", output, read]
        print(command)
        
        print("Running cutadapt round 1 for {0}...".format(sample))
        print(type(command))
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(command, stdout=f, stderr = g)
            outcode = subprocess.run(command, stdout=f, stderr=g).returncode
        f.close()
        g.close()
             
        if outcode == 0:     
            print("cutadapt round 1 finished succesfully for {0}".format(sample) +"\n"+ "check {0} for details".format(logfile))
        else:
            print("an ERROR occured during cutadapt round 1 :("+"\n"+ "check {0} for details".format(errfile)) 
            

def trim_adapters_round2(samplename):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    
    for sample in samplename:
        infile = str(outdir) + '/' + str(sample) + '_trimmed1.fq.gz'
        output = str(outdir) + '/' + str(sample) + '_trimmed_final.fq.gz'
        logfile =  str(outdir) + '/' + str(sample) + '_cutadapt_round2.log'
        errfile =  str(outdir) + '/' + str(sample) + '_cutadapt_round2.err'
        
        if cfg['cutadapt_params2'].split(',') == []:
            command = ['mv', infile, output]
            subprocess.run(command)
            
        else:
            command = ['cutadapt'] + cfg['cutadapt_params2'].split(',') +  ['-o', output, infile]
            print("Running cutadapt round 2 for {0}...".format(sample))
            #print(command)
            
            with open(logfile, 'w') as f, open(errfile, 'w') as g:
                subprocess.run(command, stdout=f, stderr = g)
                outcode = subprocess.run(command, stdout=f, stderr=g).returncode
            f.close()
            g.close()
             
            if outcode == 0:     
                print("cutadapt round 2 finished succesfully for {0}".format(sample) +"\n"+ "check {0} for details".format(logfile))
            else:
                print("an ERROR occured during cutadapt round 2 :("+"\n"+ "check {0} for details".format(errfile)) 
                

def collapse_reads(samplename, umi = None):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    
    for sample in samplename:
        infile = str(outdir) + '/' + str(sample) + '_trimmed_final.fq.gz'
        output = str(outdir) + '/' + str(sample) + '_trimmed_collapsed.fa'
        stdoutfile =  str(outdir) + '/' + str(sample) + '_trimmed.fa'
        errfile =  str(outdir) + '/' + str(sample) + '_collapse_reads.err'
        
        fq_to_fa_command = ['seqtk', 'seq', '-a', infile]
        
        print("Converting FASTQ.GZ to FASTA {0}...".format(sample))
        print(fq_to_fa_command)
        
        with open(stdoutfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(fq_to_fa_command, stdout=f, stderr = g)
            outcode = subprocess.run(fq_to_fa_command, stdout=f, stderr=g).returncode
        f.close()
        g.close()
             
        if outcode == 0:     
            print("FASTQ.GZ to FASTA conversion finished succesfully for {0}".format(sample) +"\n"+ "FASTA file {0} created".format(stdoutfile))
        else:
            print("an ERROR occured during FASTQ.GZ to FASTA conversion :("+"\n"+ "check {0} for details".format(errfile)) 
        
        print("Collapsing reads in {0}...".format(sample))
                
        counts = []
        reads = []
        
        with open(stdoutfile, 'r') as f, open(output, 'w') as g:
            for line in f:
               if not line.startswith('>') and not line in reads:
                   reads.append(line)
                   counts.append(1)
               elif not line.startswith('>') and line in reads:
                   counts[reads.index(line)] += 1

            d = OrderedDict(zip(reads, counts))
            reads_tmp = list(OrderedDict(sorted(d.items(), key=lambda x: x[1])))
            reads_final = reads_tmp[::-1]
            counts.sort(reverse = True)
            
            for idx, line in enumerate(reads_final):
                g.write('>Read' + str(idx+1) + '-' + str(counts[idx]) + '\n')
                
                if not umi:
                    g.write(line)
                elif umi > 0:
                    g.write(line[umi:])
                elif umi < 0:
                    g.write(line[0:umi-1] +'\n')
        f.close()
        g.close()
        print("Reads collapsed succesfully for {0}".format(sample)  +"\n"+ "Collapsed FASTA file {0} created".format(output))
                 
#star_rm_reps
def star_rm_reps(samplename):
    if not os.path.exists('ali'):
        os.mkdir('ali')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    rep_idx = cfg['rep_idx']
    
    for sample in samplename:
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
            print("STAR finished succesfully for {0}".format(sample) +"\n"+ "check {0} for details".format(logfile))
        else:
            print("an ERROR occured during STAR run :("+"\n"+ "check {0} for details".format(errfile)) 
            
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
            print("STAR finished succesfully for {0}".format(sample) +"\n"+ "check {0} for details".format(logfile))
        else:
            print("an ERROR occured during STAR run :("+"\n"+ "check {0} for details".format(errfile)) 
            
        os.rename('{0}Aligned.out.bam'.format(star_prefix), tmp_bam)
        bam_to_mapped_sam_command = ['samtools', 'view', '-F', '4', '-h', '-o', tmp_sam, tmp_bam]
        sorting_sam_command = ['samtools', 'sort', '-o', output, tmp_sam]
        
        print("Filtering mapped reads and converting BAM to SAM for {0}...".format(sample))
        subprocess.run(bam_to_mapped_sam_command)
        print("Filtering and conversion finished succesfully for {0}...\nSorting SAM file...".format(sample))
        subprocess.run(sorting_sam_command)
        print('Sorted SAM file created for {0}'.format(sample))
        
#config_and_run_PARalyzer  
def config_and_run_PARalyzer(samplename):

    if not os.path.exists('PARalyzer'):
        os.mkdir('PARalyzer')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'PARalyzer')
    paralyzer_mem = cfg['paralyzer_memory']
    
    for sample in samplename:
        sam_file = os.path.join(cwd,'ali', '{0}.aligned.mapped.sorted.sam'.format(sample))
        logfile =  str(outdir) + '/' + str(sample) + '_PARalyzer.log'
        errfile =  str(outdir) + '/' + str(sample) + '_PARalyzer.err'
        config_out = str(outdir) + '/' + str(sample) + '.ini'
        
        with open(config_out, 'w') as ini_file:
            ini_file.write('#\nBANDWIDTH={0}\nCONVERSION={1}\nMINIMUM_READ_COUNT_PER_GROUP={2}\nMINIMUM_READ_COUNT_PER_CLUSTER={3}\nMINIMUM_READ_COUNT_FOR_KDE={4}\nMINIMUM_CLUSTER_SIZE={5}\nMINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER={6}\nMINIMUM_CONVERSION_COUNT_FOR_CLUSTER={7}\nMINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION={8}\nMINIMUM_READ_LENGTH={9}\nMAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES={10}\n\n'.format(cfg['BANDWIDTH'],cfg['CONVERSION'],cfg['MINIMUM_READ_COUNT_PER_GROUP'],            cfg['MINIMUM_READ_COUNT_PER_CLUSTER'],cfg['MINIMUM_READ_COUNT_FOR_KDE'],cfg['MINIMUM_CLUSTER_SIZE'],cfg['MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER'],            cfg['MINIMUM_CONVERSION_COUNT_FOR_CLUSTER'],cfg['MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION'],cfg['MINIMUM_READ_LENGTH'],cfg['MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES']))
            if cfg['signal_extension'] == 'EXTEND_BY_READ':
                ini_file.write('EXTEND_BY_READ\n\n')
            elif cfg['signal_extension'] == 'HAFNER_APPROACH':
                ini_file.write('HAFNER_APPROACH\n\n')
            elif cfg['signal_extension'] == 'ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL':
                ini_file.write('ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL={0}\n\n'.format(cfg['ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL']))
            ini_file.write('GENOME_2BIT_FILE={0}\nSAM_FILE={1}\n\nOUTPUT_DISTRIBUTIONS_FILE={2}.distribution\nOUTPUT_GROUPS_FILE={2}.groups\nOUTPUT_CLUSTERS_FILE={2}.clusters\nOUTPUT_READS_FILE={2}_PARalyzer_Utilized.sam\n\n'.format(cfg['GENOME_2BIT_FILE'],sam_file,sample))
            
        paralyzer_cmd = ['PARalyzer', paralyzer_mem, config_out]   
        
        print("Running PARalyzer for {0}...".format(sample)) 
        
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(paralyzer_cmd, stdout=f, stderr = g)
            outcode = subprocess.run(paralyzer_cmd, stdout=f, stderr=g).returncode
            
        if outcode == 0:     
            print("PARalyzer finished succesfully for {0}".format(sample))
        else:
            print("an ERROR occured during PARalyzer run :("+"\n"+ "check {0} for details".format(errfile)) 
        

#run perbase
def run_perbase(samplename):

    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'PARalyzer')
    
    for sample in samplename:
        clustfile = str(outdir) + '/' + str(sample) + '.clusters'
        samfile = str(outdir) + '/' + str(sample) + '_PARalyzer_Utilized.sam'
        bamfile = str(outdir) + '/' + str(sample) + '_PARalyzer_Utilized.sorted.bam'
        bedfile = str(outdir) + '/' + str(sample) + '_perbase.bed'
        output = str(outdir) + '/' + str(sample) + '_perbase.tsv'
        logfile =  str(outdir) + '/' + str(sample) + '_perbase.log'
        errfile =  str(outdir) + '/' + str(sample) + '_perbase.err'
        genome_file = cfg['genome_dir'].split(',')
 
        samtobam = ['samtools', 'sort', '-o', bamfile, samfile]
        bam_idx = ['samtools', 'index', bamfile]
        
        print("Converting PARalyzer output *.sam to sorted *.bam for {0}...".format(sample))
        subprocess.run(samtobam)
        print("PARalyzer output converted to *.bam for {0}".format(sample))
        print("Indexing {0} bam file...".format(sample))
        subprocess.run(bam_idx)
        print("{0} bam file indexed".format(sample))
        print("Creating perbase bed from PARalyzer clusters for {0}...".format(sample))
        
        with open(clustfile, 'r') as f, open(bedfile, 'w') as g:
            
            for line in f:
                if not line.startswith('Chromosome'):
                    bed_line = []
                    clust_line = line.split(',')
                    bed_line.extend((clust_line[0], str(int(clust_line[2])-1),clust_line[3],clust_line[4],clust_line[6],clust_line[1]))
                    g.write('\t'.join(bed_line) +'\n')
                    
        print("Perbase bed created for {0}...".format(sample))
        
        perbase_cmd = ['perbase', 'base-depth', '-b', bedfile, '-r', genome_file, bamfile, '-o', output]
        
        print("Running perbase for {0}...".format(sample)) 
        
        with open(logfile, 'w') as f, open(errfile, 'w') as g:
            subprocess.run(perbase_cmd, stdout=f, stderr = g)
            outcode = subprocess.run(perbase_cmd, stdout=f, stderr=g).returncode
            
        if outcode == 0:     
            print("Perbase finished succesfully for {0}".format(sample))
        else:
            print("an ERROR occured during perbase run :("+"\n"+ "check {0} for details".format(errfile)) 

#TODO             
#run_annotR     

                
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

trim_adapters_round1(reads, samplename)
trim_adapters_round2(samplename)
collapse_reads(samplename, umi)
star_rm_reps(samplename)
star_genome_map(samplename)
config_and_run_PARalyzer(samplename)
run_perbase(samplename)

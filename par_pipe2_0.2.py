import os
import sys
import shutil
import argparse
import yaml
import concurrent.futures
import pysam
import pysamstats
import math
import wget
import time
import gzip
import pandas as pd

cwd = os.getcwd()
ts = str(time.time())

print('YOUR JOB ID IS {0}'.format(ts))
    
def read_configfile(config_file):

    print("Loading config file..." )

    with open(config_file) as f_cfg:
        cfg = yaml.load(f_cfg, Loader=yaml.FullLoader)
    print("Config file succesfully loaded")
   
    return cfg


def prepare_genome_and_annotation(configuration):
  
    if not configuration["genome_fasta"] == "" or configuration["genome_fasta"] == None: 
        genome = configuration["genome_fasta"]
        annotation = configuration["gtf"]
        outdir = os.path.dirname(os.path.abspath(genome))
        print("Genome and annotation files are provided by user")
    
    elif configuration["genome_fasta"] == "" or configuration["genome_fasta"] == None: 
        if not os.path.exists('genome'):
            os.mkdir('genome')
        outdir = os.path.join(cwd, 'genome')
        genome_file = "GRCh38.primary_assembly.genome.fa"   
        if configuration['genome_version'] == '' or configuration['genome_version'] == None:
            genome_ver = '45'    
        else:
             genome_ver = configuration['genome_version']
        if int(genome_ver) > 41:
            gtf_file = "gencode.v{0}.primary_assembly.basic.annotation.gtf".format(genome_ver)
        else:
            gtf_file = "gencode.v{0}.primary_assembly.annotation.gtf".format(genome_ver)
                 
        print("Selected genome version is {0}".format(genome_ver))
        genome = os.path.join(cwd, 'genome', genome_file)
        annotation = os.path.join(cwd, 'genome', gtf_file)
        
        if not os.path.isfile(annotation) or not os.path.isfile(genome):
            genome_gz = os.path.join(cwd, 'genome', '{0}.gz'.format(genome_file))
            gtf_gz = os.path.join(cwd, 'genome', '{0}.gz'.format(gtf_file))
            logfile =  os.path.join(cwd, 'genome', 'genome_and_gtf_unpacking.{0}.log'.format(ts))
            errfile =  os.path.join(cwd, 'genome', 'genome_and_gtf_unpacking.{0}.err'.format(ts))  
            url_genome = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{0}/{1}.gz".format(genome_ver, genome_file)          
            url_gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{0}/{1}.gz".format(genome_ver, gtf_file)
            print("Downloading GRCh38 v.{0} genome from Gencode...\n".format(genome_ver))
            wget.download(url_genome, out = outdir)
            print("\nDownloading genome annotation file...\n")
            wget.download(url_gtf, out = outdir)
            print("\nUnpacking genome and gtf files...")
            gunzip_genome_cmd = ['gunzip', genome_gz, '>', logfile, '2>', errfile] 
            gunzip_gtf_cmd = ['gunzip', gtf_gz, '>>', logfile, '2>>', errfile]
            outcode_genome = os.WEXITSTATUS(os.system(" ".join(gunzip_genome_cmd)))
            outcode_gtf = os.WEXITSTATUS(os.system(" ".join(gunzip_gtf_cmd)))
            if outcode_genome == 0:     
                print("Genome unpacked succesfully")
            else:
                print("an ERROR occured during genome unpacking :(\ncheck {0} for details".format(errfile))  
            if outcode_gtf == 0:     
                print("Gtf unpacked succesfully")
            else:
                print("An ERROR occured during gtf file unpacking :(\ncheck {0} for details".format(errfile)) 
        else:
            print('Genome and annotation files already exist')
                      
    pp2_path = configuration['parpipe2_dir'] 
    main_annot_rds = os.path.join(outdir,'main_annotation.rds')
    
    if os.path.isfile(main_annot_rds) and os.path.getmtime(main_annot_rds) > os.path.getmtime(annotation):
        print('main_annotation.rds file from current gtf already exist')
    else:
        print("Preparing main_annotation.rds file")
        logfile =  os.path.join(outdir, 'prep_main_rds_from_gtf.{0}.log'.format(ts))
        errfile =  os.path.join(outdir, 'prep_main_rds_from_gtf.{0}.err'.format(ts))  
        ann_script = os.path.join(pp2_path,"bin","prepareAnnotFilesFromGtf.R")
        ann_cmd = ['Rscript', ann_script, outdir, annotation, 'main', '>', logfile, '2>', errfile]
        outcode = os.WEXITSTATUS(os.system(" ".join(ann_cmd)))
    
        if outcode == 0:     
            print("main_annotation.rds file created succesfully\ncheck {0} for details".format(logfile))
        else:
            print("an ERROR occured during creating main_annotation.rds file :(\ncheck {0} for details".format(errfile)) 
    
    genome_2bit = os.path.normpath(str(genome).replace('fasta','2bit').replace('fa','2bit'))        
    
    if os.path.isfile(genome_2bit) and os.path.isfile('{0}.fai'.format(genome)):
        print('Index and 2bit files exist. Moving to the next step')
    else:
        logfile =  os.path.join(cwd, 'genome', 'genome_indexing_and_2bit_preparation.{0}.log'.format(ts))
        errfile =  os.path.join(cwd, 'genome', 'genome_indexing_and_2bit_preparation.{0}.err'.format(ts))  
        index_genome_cmd = ['samtools', 'faidx', genome, '>', logfile, '2>', errfile]
        command_2bit = ['faToTwoBit', genome, genome_2bit, '>>', logfile, '2>>', errfile]
        outcode = os.WEXITSTATUS(os.system(" ".join(index_genome_cmd)))
        if outcode == 0:     
            print("Genome file indexed succesfully\ncheck {0} for details".format(logfile))
        else:
            print("an ERROR occured during genome indexing :(\ncheck {0} for details".format(errfile))
        print("Preparing genome 2bit file...")
        outcode = os.WEXITSTATUS(os.system(" ".join(command_2bit)))
        if outcode == 0:     
            print("Genome 2bit file prepared succesfully\ncheck {0} for details".format(logfile))
        else:
            print("an ERROR occured during genome 2bit preparation :(\ncheck {0} for details".format(errfile))
            
    return genome, annotation, genome_2bit, main_annot_rds
    
    
def prepare_star_index(configuration, genome, annotation, threads):
    
    genome = genome
    annotation = annotation
    if configuration['star_index_dir'] =="" or configuration['star_index_dir'] == None:
        if not os.path.exists('star_index'):
            os.mkdir('star_index')
        star_idx_outdir = os.path.join(cwd, 'star_index')     
        logfile =  os.path.join(cwd, 'star_index', 'star_index.{0}.log'.format(ts))
        errfile =  os.path.join(cwd, 'star_index', 'star_index{0}.err'.format(ts)) 
        if configuration['star_index_params'] == '' or configuration['star_index_params'] == None:
            command = ['STAR', '--runMode genomeGenerate --runThreadN', str(threads), '--genomeDir', star_idx_outdir, '--genomeFastaFiles', genome, '--sjdbGTFfile', annotation, '--sjdbOverhang 100', '>', logfile, '2>', errfile]
        else:
            command = ['STAR', '--runMode genomeGenerate --runThreadN',str(threads),'--genomeDir', star_idx_outdir, '--genomeFastaFiles', genome, '--sjdbGTFfile', annotation] + cfg['star_index_params'].split(",") + ['>', logfile, '2>', errfile]
            
        #print(command)
    
        print("Generating STAR index...")    
        outcode = os.WEXITSTATUS(os.system(" ".join(command)))
        if outcode == 0:     
            print("STAR index generated succesfully\ncheck {0} for details".format(logfile))
        else:
            print("an ERROR occured during STAR index generation :(\ncheck {0} for details".format(errfile))
    else:
        print("STAR index provided by user")
        star_idx_outdir = configuration['star_index_dir']
            
    return star_idx_outdir
    
    
def prepare_bowtie_index(configuration, genome, threads):

    genome = genome
    if configuration['bowtie_index_dir'] =="" or configuration['bowtie_index_dir'] == None:
        if not os.path.exists('bowtie_index'):
            os.mkdir('bowtie_index')
        bowtie_idx_outdir = os.path.join(cwd, 'bowtie_index')     
        logfile =  os.path.join(cwd, 'bowtie_index', 'bowtie_index.{0}.log'.format(ts))
        errfile =  os.path.join(cwd, 'bowtie_index', 'bowtie_index{0}.err'.format(ts)) 
        if configuration['bowtie_index_params'] == '' or configuration['bowtie_index_params'] == None:
            command = ['bowtie-build', '--threads', str(threads), '-f', genome, '{0}/bwt'.format(str(bowtie_idx_outdir)), '>', logfile, '2>', errfile]
        else:
            command = ['bowtie-build'] + cfg['bowtie_index_params'].split(",") + ['-f', genome, '{0}/bwt'.format(str(bowtie_idx_outdir)), '>', logfile, '2>', errfile]
        
        print("Generating bowtie index...")    
        outcode = os.WEXITSTATUS(os.system(" ".join(command)))
        if outcode == 0:     
            print("bowtie index generated succesfully\ncheck {0} for details".format(logfile))
        else:
            print("an ERROR occured during bowtie index generation :(\ncheck {0} for details".format(errfile))
    else:
        print("bowtie index provided by user")
        bowtie_idx_outdir = configuration['bowtie_index_dir']
            
    return bowtie_idx_outdir
    
    
def trim_adapters_round1(configuration, reads, samplename):
    if not os.path.exists('reads'):
        os.mkdir('reads')

    outdir = os.path.join(cwd, 'reads')
    
    samplename = samplename.split(",") 
    reads = reads.split(",")
    out_trim_reads1_list = []
    cutadapt_log_list = []
    
    for idx, sample in enumerate(samplename):
        
        read = reads[idx]
        out_trim_reads1 = os.path.join(cwd, 'reads', '{0}_{1}.trimmed1.fq.gz'.format(sample,ts))
        out_trim_reads1_list.append(out_trim_reads1)
        logfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round1.{1}.log'.format(sample,ts))
        cutadapt_log_list.append(logfile)
        errfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round1.{1}.err'.format(sample,ts))
        if configuration['nextseq'] in ['', None, False, 'False', 0, '0']: 
            if not configuration['three_prime_adapter'] in ['', None, False, 'False', 0, '0'] and not configuration['five_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-a', configuration['three_prime_adapter'], '-g', configuration['five_prime_adapter']] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            elif not configuration['three_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-a', configuration['three_prime_adapter']] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            elif not configuration['five_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-g', configuration['five_prime_adapter']] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            else:
                command = ['cutadapt'] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
        else:
            if not configuration['three_prime_adapter'] in ['', None, False, 'False', 0, '0'] and not configuration['five_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-a', configuration['three_prime_adapter'], '-g', configuration['five_prime_adapter'], '--nextseq-trim=15'] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            elif not configuration['three_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-a', configuration['three_prime_adapter'], '--nextseq-trim=15'] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            elif not configuration['five_prime_adapter'] in ['', None, False, 'False', 0, '0']:
                command = ['cutadapt', '-g', configuration['five_prime_adapter'], '--nextseq-trim=15'] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
            else:
                command = ['cutadapt', '--nextseq-trim=15'] + configuration['cutadapt_params1'].split(",") +  ["-o", out_trim_reads1, read, '>', logfile, '2>', errfile]
        
        print("Running cutadapt round 1 for {0}...".format(sample))
        outcode = os.WEXITSTATUS(os.system(" ".join(command)))
        if outcode == 0:     
            print("cutadapt round 1 finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during cutadapt round 1 :(\ncheck {0} for details".format(errfile)) 
      
    return out_trim_reads1_list, cutadapt_log_list
    
    
def trim_adapters_round2(configuration, samplename, out_trim_reads1_list):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    out_trim_reads2_list = []
    
    for idx, sample in enumerate(samplename):
        infile = out_trim_reads1_list[idx]
        out_trim_reads2 = os.path.join(cwd, 'reads', '{0}_{1}.trimmed_final.fq.gz'.format(sample,ts))
        out_trim_reads2_list.append(out_trim_reads2)
        logfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round2.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd, 'reads', '{0}_cutadapt_round2.{1}.err'.format(sample,ts))
        
        if configuration['cutadapt_params2'] == '' or configuration['cutadapt_params2'] == None:
            os.rename(infile, out_trim_reads2)
            print("Skipping cutadapt round 2 for {0}...".format(sample))
            
        else:
            command = ['cutadapt'] + configuration['cutadapt_params2'].split(',') +  ['-o', out_trim_reads2, infile, '>', logfile, '2>', errfile]
            print("Running cutadapt round 2 for {0}...".format(sample))
            
            outcode = os.WEXITSTATUS(os.system(" ".join(command))) 
            if outcode == 0:     
                print("cutadapt round 2 finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
            else:
                print("an ERROR occured during cutadapt round 2 :(\ncheck {0} for details".format(errfile)) 
            
    return out_trim_reads2_list
    
def collapse_reads(samplename, out_trim_reads2_list, umi = None):
    
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'reads')
    collapsed_reads_list = []
    
    for idx, sample in enumerate(samplename):
        infile = out_trim_reads2_list[idx]
        collapsed_reads = os.path.join(cwd, 'reads', '{0}_{1}.trimmed_collapsed.fa'.format(sample,ts))
        collapsed_reads_list.append(collapsed_reads)
        stdoutfile =  os.path.join(cwd, 'reads', '{0}_{1}.trimmed.fa'.format(sample,ts))
        errfile =  os.path.join(cwd, 'reads', '{0}_collapsed_reads.{1}.err'.format(sample, ts))
        
        fq_to_fa_command = ['seqtk', 'seq', '-a', infile, '>', stdoutfile, '2>', errfile]
        
        print("Converting FASTQ.GZ to FASTA for {0}...".format(sample))
        
        outcode = os.WEXITSTATUS(os.system(" ".join(fq_to_fa_command)))    
        if outcode == 0:     
            print("FASTQ.GZ to FASTA conversion finished succesfully for {0}\nFASTA file {1} created".format(sample,stdoutfile))
        else:
            print("an ERROR occured during FASTQ.GZ to FASTA conversion :(\ncheck {0} for details".format(errfile)) 
        
        print("Collapsing reads in {0}...".format(sample))
                
        if not umi:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "sed", "'s/>/>Read/g'", ">", collapsed_reads]
        elif umi > 0:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "fastx_trimmer", "-f", str(umi+1), "|", "sed", "'s/>/>Read/g'", ">", collapsed_reads, '2>', errfile]
        elif umi < 0:
            collapse_cmd = ["fastx_collapser", "-i", stdoutfile, "|", "fastx_trimmer", "-t", str(-umi), "|", "sed", "'s/>/>Read/g'", ">", collapsed_reads, '2>', errfile]
        
        outcode = os.WEXITSTATUS(os.system(" ".join(collapse_cmd)))  
        if outcode == 0:
            print("Reads collapsed succesfully for {0}\nCollapsed FASTA file for {0} created".format(sample))
        else:
            print("an ERROR occured during reads collapsing :(\ncheck {0} for details".format(errfile)) 
        
    return collapsed_reads_list
        
        
def star_rm_reps(configuration, samplename, collapsed_reads_list, threads):
    if not os.path.exists('ali'):
        os.mkdir('ali')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    rm_rep_reads_list = []
        
    for idx, sample in enumerate(samplename):
    
        infile = collapsed_reads_list[idx]
        rm_rep_reads = os.path.join(cwd, 'reads', '{0}_{1}.rmreps.fa'.format(sample,ts))
        rm_rep_reads_list.append(rm_rep_reads)
        star_prefix = os.path.join(cwd, 'ali', '{0}_{1}.rmreps'.format(sample,ts))
    
        if configuration['rem_rep'] in [False, 'False', '0', 0, 'no']:
            print('Skipping removing repetitive elements')
            shutil.copy2(infile, rm_rep_reads)    
               
        else:        
            logfile =  os.path.join(cwd, 'reads', '{0}_rmreps.{1}.log'.format(sample,ts))
            errfile =  os.path.join(cwd, 'reads', '{0}_rmreps.{1}.err'.format(sample,ts))
            
            if configuration['rep_idx'] == '' or configuration['rep_idx'] == None:
                pp2_path = configuration['parpipe2_dir'] 
                rep_idx = os.path.join(pp2_path, 'annotation', 'rep_idx')
                if os.path.isfile('{0}.tar.gz'.format(rep_idx)):
                    gz_cmd = ['tar', '-xvzf', '{0}.tar.gz'.format(rep_idx), '-C', os.path.join(pp2_path, 'annotation'), '&&', 'rm', '{0}.tar.gz'.format(rep_idx),  '>>', logfile, '2>>', errfile]
                    print('Unzipping rep index...')
                    outcode = os.WEXITSTATUS(os.system(" ".join(gz_cmd)))
                    if outcode == 0:     
                        print("Rep index unzipped succesfully")
                    else:
                        print("an ERROR occured during rep index unzipping :(\ncheck {0} for details".format(errfile)) 
            else:
                rep_idx = configuration['rep_idx'] 
                
            if configuration['star_rem_reps_params'] == '' or configuration['star_rem_reps_params'] == None:
                command = ['STAR'] + "--runMode,alignReads,--genomeLoad,NoSharedMemory,--alignEndsType,EndToEnd,--outSAMunmapped,Within,--outFilterMultimapNmax,30,--outFilterMultimapScoreRange,1,--outSAMtype,BAM,Unsorted,--outFilterType,BySJout,--outBAMcompression,10,--outReadsUnmapped,Fastx,--outFilterScoreMin,10,--outSAMattrRGline,ID:foo,--outSAMattributes,All,--outSAMmode,Full".split(",") + ['--runThreadN', str(threads), '--genomeDir', rep_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile, '>>', logfile, '2>>', errfile]
                
            else:
                command = ['STAR'] + configuration['star_rem_reps_params'].split(',') + ['--runThreadN', str(threads), '--genomeDir', rep_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile, '>>', logfile, '2>>', errfile]
        
            print("Running STAR to remove repetitive sequences for {0}...".format(sample))
        
            outcode = os.WEXITSTATUS(os.system(" ".join(command)))
             
            if outcode == 0:     
                print("STAR finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
            else:
                print("an ERROR occured during STAR run :(\ncheck {0} for details".format(errfile)) 
            
            os.rename('{0}Unmapped.out.mate1'.format(star_prefix), rm_rep_reads)
            os.rename('{0}Aligned.out.bam'.format(star_prefix), os.path.join(cwd, 'ali', '{0}.{1}.reps.bam'.format(sample,ts)))
            
    return rm_rep_reads_list
            
    
def star_genome_map(configuration, samplename, star_index_outdir, rm_rep_reads_list, threads):

    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    genome_idx = star_index_outdir
    output_sam_list = []
    ali_log_list = []
    
    for idx, sample in enumerate(samplename):
        infile = rm_rep_reads_list[idx]
        output_sam = os.path.join(cwd, 'ali', '{0}.{1}.aligned.mapped.sorted.sam'.format(sample,ts))
        output_sam_list.append(output_sam)
        star_prefix = os.path.join(cwd, 'ali', '{0}.{1}.'.format(sample,ts))
        tmp_bam = os.path.join(cwd, 'ali', '{0}.{1}.aligned.bam'.format(sample,ts))
        tmp_sam = os.path.join(cwd, 'ali', '{0}.{1}.aligned.mapped.sam'.format(sample,ts))
        logfile =  os.path.join(cwd, 'ali', '{0}_genome_map.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd, 'ali', '{0}_genome_map.{1}.err'.format(sample,ts))
        log_star = os.path.join(cwd, 'ali', '{0}.{1}.Log.final.out'.format(sample,ts))
        ali_log_list.append(log_star)
        
        if configuration['star_map_params'] == '' or configuration['star_map_params'] == None:
            star_command = ['STAR'] + "--runMode,alignReads,--genomeLoad,NoSharedMemory,--alignEndsType,EndToEnd,--outSAMunmapped,Within,--outFilterMultimapNmax,10,--outSAMmultNmax,1,--outSAMattributes,All,--outSAMtype,BAM,Unsorted,--outFilterType,BySJout,--outFilterMismatchNoverReadLmax,0.1,--outSAMattrRGline,ID:foo,--outStd,Log,--outBAMcompression,10,--outSAMmode,Full".split(",") + ['--runThreadN', str(threads), '--genomeDir', genome_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile, '>', logfile, '2>', errfile]
        else:
            star_command = ['STAR'] + configuration['star_map_params'].split(',') + ['--runThreadN', str(threads), '--genomeDir', genome_idx, '--outFileNamePrefix', star_prefix, '--readFilesIn', infile, '>', logfile, '2>', errfile]
        
        print("Running STAR to align to genome for {0}...".format(sample))
        
        outcode = os.WEXITSTATUS(os.system(" ".join(star_command)))    
        if outcode == 0:     
            print("STAR finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during STAR run :(\ncheck {0} for details".format(errfile)) 
            
        os.rename('{0}Aligned.out.bam'.format(star_prefix), tmp_bam)
        bam_to_mapped_sam_command = ['samtools', 'view', '-F', '4', '-h', '-o', tmp_sam, tmp_bam, '>>', logfile, '2>>', errfile]
        sorting_sam_command = ['samtools', 'sort', '-o', output_sam, tmp_sam, '>>', logfile, '2>>', errfile]
        
        print("Filtering mapped reads and converting BAM to SAM for {0}...".format(sample))
        outcode = os.WEXITSTATUS(os.system(" ".join(bam_to_mapped_sam_command)))    
        if outcode == 0:     
            print("Mapped SAM file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during mapped SAM file creation :(\ncheck {0} for details".format(errfile)) 
        outcode = os.WEXITSTATUS(os.system(" ".join(sorting_sam_command))) 
        if outcode == 0:     
            print("Mapped SAM file sorted succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during mapped SAM file sorting :(\ncheck {0} for details".format(errfile)) 
        
    return output_sam_list, ali_log_list 
    
    
def bowtie_genome_map(configuration, samplename, bowtie_index_outdir, rm_rep_reads_list, threads):

    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'ali')
    genome_idx = '{0}/bwt'.format(str(bowtie_index_outdir))
    output_sam_list = []
    ali_log_list = []
    
    for idx, sample in enumerate(samplename):
        infile = rm_rep_reads_list[idx]
        output_sam = os.path.join(cwd, 'ali', '{0}.{1}.aligned.mapped.sorted.sam'.format(sample,ts))
        output_sam_list.append(output_sam)
        logfile =  os.path.join(cwd, 'ali', '{0}_genome_map.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd, 'ali', '{0}_genome_map.{1}.err'.format(sample,ts))
        ali_log_list.append(errfile)
        
        if configuration['bowtie_map_params'] == '' or configuration['bowtie_map_params'] == None:
            bowtie_command = ['bowtie', genome_idx, '-v 1 -m 10 --best --strata -p', str(threads), '-S -f', infile, '--quiet', '2>>', errfile, '| samtools view -hS -F 4 - | samtools sort -O SAM - -o', output_sam, '>>', logfile, '2>>', errfile]
        else:
            bowtie_command = ['bowtie', genome_idx] + configuration['bowtie_map_params'].split(',') + ['-f', infile,  '2>>', errfile, '| samtools view -hS -F 4 - | samtools sort -O SAM - -o', output_sam, '>>', logfile, '2>>', errfile]
            
        print("Running bowtie to align to genome for {0}...".format(sample))
        
        outcode = os.WEXITSTATUS(os.system(" ".join(bowtie_command)))    
        if outcode == 0:     
            print("Bowtie finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during bowtie run :(\ncheck {0} for details".format(errfile)) 
            
    return output_sam_list, ali_log_list
    
def annotate_mapped_reads(configuration, samplename, output_sam_list, main_annot_rds):
    
    samplename = samplename.split(",") 
    if not os.path.exists('total_conv_stats'):
        os.mkdir('total_conv_stats')
    outdir = os.path.join(cwd, 'total_conv_stats')
    #d0_list = []
    #d1_list = []
    #d2_list = []
    full_bed_mm_list = []
    merged_sam_mm_list = []
    merged_sam_list = []
    out_genome_annot_list = []
    #d0_list_bed = []
    #d1_list_bed = []
    #d2_list_bed = []
    genome_annotR = os.path.join(pp2_dir, 'bin', 'genome_annot.R')
    pp2_dir = configuration['parpipe2_dir']
    if not configuration['annot_rank'] in ['', None, False, 'False', 0, '0']: 
        annot_rank_file = configuration['annot_rank']
    else:
        annot_rank_file = os.path.join(pp2_dir, 'annotation', 'annot_categories.txt')
    annotR = os.path.join(pp2_dir, 'bin', 'annot.R')
    main_annot = main_annot_rds
    trna_annot = os.path.join(pp2_dir, 'annotation', 'rds_files','trna.rds')
    repeats_annot = os.path.join(pp2_dir, 'annotation', 'rds_files','repeats.rds')
    
    for idx, sample in enumerate(samplename):
        infile = output_sam_list[idx]
        logfile =  os.path.join(cwd, 'total_conv_stats', '{0}_annot_mapped_reads.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd, 'total_conv_stats', '{0}_annot_mapped_reads.{1}.err'.format(sample,ts))
        out_d0 = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d0.sam'.format(sample,ts))
        out_d1 = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d1.sam'.format(sample,ts))
        out_d2 = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d2.sam'.format(sample,ts))
        out_d0_bed = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d0.bed'.format(sample,ts))
        out_d1_bed = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d1.bed'.format(sample,ts))
        out_d2_bed = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.d2.bed'.format(sample,ts))
        full_bed_mm = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.mm.bed'.format(sample,ts))
        merged_sam_file = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.merged_sorted.sam'.format(sample,ts))
        merged_sam_mm_file = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.mm_merged_sorted.sam'.format(sample,ts))
        out_genome_annot = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.genome_annot.tsv'.format(sample,ts))
        #d0_list.append(out_d0)
        #d1_list.append(out_d1)
        #d2_list.append(out_d2)
        full_bed_list.append(full_bed)
        merged_sam_list.append(merged_sam_file)
        merged_sam_mm_list.append(merged_sam_mm_file)
        out_genome_annot_list.append(out_genome_annot)
        #d0_list_bed.append(out_d0_bed)
        #d1_list_bed.append(out_d1_bed)
        #d2_list_bed.append(out_d2_bed)
        command_d0 = ["samtools view -h -e '[NM]==0' -O SAM -o", out_d0, infile, '>', logfile, '2>', errfile]
        command_d1 = ["samtools view -h -e '[NM]==1' -O SAM -o", out_d1, infile, '>>', logfile, '2>>', errfile]
        command_d2 = ["samtools view -h -e '[NM]==2' -O SAM -o", out_d2, infile, '>>', logfile, '2>>', errfile]
        command_d0_bed = ["bedtools bamtobed -i", out_d0, '>', out_d0_bed , '2>>', errfile]
        command_d1_bed = ["bedtools bamtobed -i", out_d1, '>', out_d1_bed , '2>>', errfile]
        command_d2_bed = ["bedtools bamtobed -i", out_d2, '>', out_d2_bed , '2>>', errfile]
        command_merge_sam = ['samtools merge -', out_d0, out_d1, out_d2, '| samtools sort -O sam -o', merged_sam_file]
        command_merge_sam_mm = ['samtools merge -', out_d1, out_d2, '| samtools sort -O sam -o', merged_sam_mm_file]
        command_full_bed = ["bedtools bamtobed -i", merged_sam_mm_file, '| bedtools merge -s -c 6 -o distinct >', full_bed_mm , '2>>', errfile]
        command_annot = ['Rscript', genome_annotR, outdir, annot_rank_file, main_annot, trna_annot, repeats_annot, out_d0_bed, out_d1_bed, out_d2_bed, out_genome_annot,'>', logfile, '2>', errfile]
        print("Spliting reads according to edit distance for {0}...".format(sample))
        outcode_d0 = os.WEXITSTATUS(os.system(" ".join(command_d0)))    
        if outcode_d0 == 0:     
            print("Distance 0 SAM file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 0 SAM file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_d0_bed = os.WEXITSTATUS(os.system(" ".join(command_d0_bed)))    
        if outcode_d0_bed == 0:     
            print("Distance 0 BED file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 0 BED file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_d1 = os.WEXITSTATUS(os.system(" ".join(command_d1)))    
        if outcode_d1 == 0:     
            print("Distance 1 SAM file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 1 SAM file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_d1_bed = os.WEXITSTATUS(os.system(" ".join(command_d1_bed)))    
        if outcode_d1_bed == 0:     
            print("Distance 1 BED file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 1 BED file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_d2 = os.WEXITSTATUS(os.system(" ".join(command_d2)))    
        if outcode_d2 == 0:     
            print("Distance 2 SAM file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 2 SAM file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_d2_bed = os.WEXITSTATUS(os.system(" ".join(command_d2_bed)))    
        if outcode_d2_bed == 0:     
            print("Distance 2 BED file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during distance 2 BED file creation :(\ncheck {0} for details".format(errfile)) 
        outcode_merge_sam = os.WEXITSTATUS(os.system(" ".join(command_merge_sam)))
        if outcode_merge_sam == 0:     
            print("Merged SAM file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during SAM file merging :(\ncheck {0} for details".format(errfile))  
        outcode_full_bed = os.WEXITSTATUS(os.system(" ".join(command_full_bed)))   
        if outcode_full_bed == 0:     
            print("Merged BED file created succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during BED file merging :(\ncheck {0} for details".format(errfile)) 
        print("Running genome annotation for {0}...".format(sample))
        outcode_ga = os.WEXITSTATUS(os.system(" ".join(command_annot))) 
            
    return merged_sam_list, full_bed_list, out_genome_annot_list #,annot_genome (from r script)
    #merge samfiles and use it as input
        

        
def genome_conv_stats(samplename, merged_sam_list, full_bed_list, out_genome_annot_list, genome):#, annot_genome):
    samplename = samplename.split(",") 
    total_conv_list = []
    
    for idx, sample in enumerate(samplename):
        samfile = merged_sam_list[idx]
        bamfile = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.merged_sorted.bam'.format(sample,ts))
        bedfile = full_bed_list[idx]
        output_total_conv = os.path.join(cwd, 'total_conv_stats', '{0}.{1}.total_conv_stats.tsv'.format(sample,ts))
        total_conv_list.append(output_total_conv)
        logfile =  os.path.join(cwd, 'total_conv_stats', '{0}_total_conv_stats.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd, 'total_conv_stats', '{0}_total_conv_stats.{1}.err'.format(sample,ts))
        genome_ref = genome
        
        command_sam_bam = ["samtools view -bS  -o", bamfile, samfile, '>', logfile, '2>', errfile]
        command_bam_idx = ["samtools index", bamfile, '>>', logfile, '2>>', errfile]
        print("Converting meged SAM to merged BAM for {0}".format(sample))
        outcode_sam_bam = os.WEXITSTATUS(os.system(" ".join(command_sam_bam)))   
        if outcode_sam_bam == 0:     
            print("Coversion finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during SAM to BAM conversion :(\ncheck {0} for details".format(errfile)) 
        print("Indexing merged BAM for {0}".format(sample))
        outcode_bam_idx = os.WEXITSTATUS(os.system(" ".join(command_bam_idx)))   
        if outcode_bam_idx == 0:     
            print("Indexing finished succesfully for {0}\ncheck {1} for details".format(sample,logfile))
        else:
            print("an ERROR occured during BAM indexing :(\ncheck {0} for details".format(errfile)) 
        
        
        conv_stats = {"AC": 0,
        "AG" : 0,
        "AT" : 0,
        "Adel": 0,
        "CA" : 0,
        "CG" : 0,
        "CT": 0,
        "Cdel" : 0,
        "GA" : 0,
        "GC" : 0,
        "GT" : 0,
        "Gdel" : 0,
        "TA" : 0,
        "TC" : 0,
        "TG" : 0,
        "Tdel": 0,
        "other" : 0
        }
        
        bam = pysam.AlignmentFile(bamfile, "rb")
        print("Generating genome-wide conversion statistics... (this step may be quite long)...")
        with open (bedfile, 'r') as f:
            for line in f:
                bed_line = line.split('\t')
                
                for rec in pysamstats.stat_variation_strand(bam, chrom = bed_line[0], start=int(bed_line[1]), end = int(bed_line[2]), pad=True, truncate=True, fafile = genome_ref, one_based=False):
                    if bed_line[3] == "+":
                        if rec['ref'] == 'A':
                            conv_stats["AC"] += int(rec['C_fwd'])
                            conv_stats["AG"] += int(rec['G_fwd']) 
                            conv_stats["AT"] += int(rec['T_fwd'])
                            conv_stats["Adel"] += int(rec['deletions_fwd'])
                        elif rec['ref'] == 'C':
                            conv_stats["CA"] += int(rec['A_fwd'])
                            conv_stats["CG"] += int(rec['G_fwd']) 
                            conv_stats["CT"] += int(rec['T_fwd'])
                            conv_stats["Cdel"] += int(rec['deletions_fwd'])
                        elif rec['ref'] == 'G':
                            conv_stats["GA"] += int(rec['A_fwd'])
                            conv_stats["GC"] += int(rec['C_fwd']) 
                            conv_stats["GT"] += int(rec['T_fwd'])
                            conv_stats["Gdel"] += int(rec['deletions_fwd'])
                        elif rec['ref'] == 'T':
                            conv_stats["TA"] += int(rec['A_fwd'])
                            conv_stats["TC"] += int(rec['C_fwd']) 
                            conv_stats["TG"] += int(rec['G_fwd'])
                            conv_stats["Tdel"] += int(rec['deletions_fwd'])
                    else:
                        if rec['ref'] == 'T':
                            conv_stats["AC"] += int(rec['G_rev'])
                            conv_stats["AG"] += int(rec['C_rev']) 
                            conv_stats["AT"] += int(rec['A_rev'])
                            conv_stats["Adel"] += int(rec['deletions_rev'])
                        elif rec['ref'] == 'G':
                            conv_stats["CA"] += int(rec['T_rev'])
                            conv_stats["CG"] += int(rec['C_rev']) 
                            conv_stats["CT"] += int(rec['A_rev'])
                            conv_stats["Cdel"] += int(rec['deletions_rev'])
                        elif rec['ref'] == 'C':
                            conv_stats["GA"] += int(rec['T_rev'])
                            conv_stats["GC"] += int(rec['G_rev']) 
                            conv_stats["GT"] += int(rec['A_rev'])
                            conv_stats["Gdel"] += int(rec['deletions_rev'])
                        elif rec['ref'] == 'A':
                            conv_stats["TA"] += int(rec['T_rev'])
                            conv_stats["TC"] += int(rec['G_rev']) 
                            conv_stats["TG"] += int(rec['C_rev'])
                            conv_stats["Tdel"] += int(rec['deletions_rev'])



        bam.close()

         
        conv_tab = pd.DataFrame.from_dict(conv_stats, orient = 'index', columns = ['counts'])
        conv_tab.index.name = 'conversion'
        conv_tab.reset_index(inplace=True)
        conv_tab.to_csv(output_total_conv, sep = "\t", index=False)

        print("Genome-wide conversion statistics generated succesfully. Check {0} for details".format(output_total_conv))
def config_and_run_PARalyzer(configuration, samplename, output_sam_list, genome_2bit):

    if not os.path.exists('PARalyzer'):
        os.mkdir('PARalyzer')
    samplename = samplename.split(",") 
    outdir = os.path.join(cwd, 'PARalyzer')
    paralyzer_mem = configuration['paralyzer_memory']
    output_clust_list = []
    output_groups_list = []
    
    for idx, sample in enumerate(samplename):
        sam_file = output_sam_list[idx]
        output_clust = os.path.join(cwd,outdir, '{0}.{1}.clusters'.format(sample,ts))
        output_clust_list.append(output_clust)
        output_group = os.path.join(cwd,outdir, '{0}.{1}.groups'.format(sample,ts))
        output_groups_list.append(output_group)
        genome_2bit = genome_2bit
        logfile =  os.path.join(cwd,outdir, '{0}_PARalyzer.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd,outdir, '{0}_PARalyzer.{1}.err'.format(sample,ts))
        config_out = os.path.join(cwd,outdir, '{0}.{1}.ini'.format(sample,ts))
        
        with open(config_out, 'w') as ini_file:
            ini_file.write('#\nBANDWIDTH={0}\nCONVERSION={1}\nMINIMUM_READ_COUNT_PER_GROUP={2}\nMINIMUM_READ_COUNT_PER_CLUSTER={3}\nMINIMUM_READ_COUNT_FOR_KDE={4}\nMINIMUM_CLUSTER_SIZE={5}\nMINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER={6}\nMINIMUM_CONVERSION_COUNT_FOR_CLUSTER={7}\nMINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION={8}\nMINIMUM_READ_LENGTH={9}\nMAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES={10}\n\n'.format(configuration['BANDWIDTH'],configuration['CONVERSION'],configuration['MINIMUM_READ_COUNT_PER_GROUP'],            configuration['MINIMUM_READ_COUNT_PER_CLUSTER'],configuration['MINIMUM_READ_COUNT_FOR_KDE'],configuration['MINIMUM_CLUSTER_SIZE'],configuration['MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER'],            configuration['MINIMUM_CONVERSION_COUNT_FOR_CLUSTER'],configuration['MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION'],configuration['MINIMUM_READ_LENGTH'],configuration['MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES']))
            if configuration['signal_extension'] == 'EXTEND_BY_READ':
                ini_file.write('EXTEND_BY_READ\n\n')
            elif configuration['signal_extension'] == 'HAFNER_APPROACH':
                ini_file.write('HAFNER_APPROACH\n\n')
            elif configuration['signal_extension'] == 'ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL':
                ini_file.write('ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL={0}\n\n'.format(configuration['ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL']))
            if not configuration['use_pcr_duplicates'] in ["TRUE", "True", True, "true", "T"]:
                ini_file.write('GENOME_2BIT_FILE={0}\nSAM_FILE={1}\n\nOUTPUT_DISTRIBUTIONS_FILE={2}.distribution\nOUTPUT_GROUPS_FILE={2}.groups\nOUTPUT_CLUSTERS_FILE={2}.clusters\nOUTPUT_READS_FILE={2}_PARalyzer_Utilized.sam\n\n'.format(genome_2bit,sam_file,str(output_clust).replace('.clusters','')))
            else: 
                ini_file.write('GENOME_2BIT_FILE={0}\nSAM_FILE={1}=COLLAPSED\n\nOUTPUT_DISTRIBUTIONS_FILE={2}.distribution\nOUTPUT_GROUPS_FILE={2}.groups\nOUTPUT_CLUSTERS_FILE={2}.clusters\nOUTPUT_READS_FILE={2}_PARalyzer_Utilized.sam\n\n'.format(genome_2bit,sam_file,str(output_clust).replace('.clusters','')))
            
        paralyzer_cmd = ['PARalyzer', paralyzer_mem, config_out, '>', logfile, '2>', errfile]   
        
        print("Running PARalyzer for {0}...".format(sample)) 
        
        outcode = os.WEXITSTATUS(os.system(" ".join(paralyzer_cmd)))
            
        if outcode == 0:     
            print("PARalyzer finished succesfully for {0}".format(sample))
        else:
            print("an ERROR occured during PARalyzer run :(\ncheck {0} for details".format(errfile))   
            
    return output_groups_list, output_clust_list
    
    
def run_clusters_stats(samplename, output_sam_list, output_clust_list, genome):

    samplename = samplename.split(",") 
    output_stats_list = []
    
    for idx, sample in enumerate(samplename):
        clustfile = output_clust_list[idx]
        samfile = output_sam_list[idx]
        bamfile = os.path.join(cwd,'PARalyzer', '{0}_{1}.PARalyzer_Utilized.sorted.bam'.format(sample,ts))
        output_stats = os.path.join(cwd,'PARalyzer', '{0}_{1}.clusters_conv_stats.tsv'.format(sample,ts))
        output_stats_list.append(output_stats)
        logfile =  os.path.join(cwd,'PARalyzer', '{0}_bam_prep.{1}.log'.format(sample,ts))
        errfile =  os.path.join(cwd,'PARalyzer', '{0}_bam_prep.{1}.err'.format(sample,ts))
        genome_ref = genome
 
        samtobam = ['samtools', 'sort', '-o', bamfile, samfile, '>', logfile, '2>', errfile]
        bam_idx = ['samtools', 'index', bamfile, '>>', logfile, '2>>', errfile]
        
        print("Converting PARalyzer output *.sam to sorted *.bam for {0}...".format(sample))
        outcode = os.WEXITSTATUS(os.system(" ".join(samtobam)))
        if outcode == 0:     
            print("PARalyzer outbut SAM succesfully converted to BAM for {0}".format(sample))
        else:
            print("an ERROR occured during PARalyzer outbut SAM conversion to BAM :(\ncheck {0} for details".format(errfile)) 
        print("Indexing {0} bam file...".format(sample))
        outcode = os.WEXITSTATUS(os.system(" ".join(bam_idx)))
        if outcode == 0:     
            print("BAM file succesfully indexed for {0}".format(sample))
        else:
            print("an ERROR occured during BAM file indexing :(\ncheck {0} for details".format(errfile)) 
        print("Adding stats to PARalyzer output for {0}...".format(sample))
        
        bam = pysam.AlignmentFile(bamfile, "rb")
        
        with open(clustfile, 'r') as f, open(output_stats, 'w') as g:
           
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
        
    return output_stats_list
    
    
def run_annot(configuration, samplename, output_stats_list, main_annot_rds):

    if not os.path.exists('output'):
        os.mkdir('output')
    outdir = os.path.join(cwd, 'output')
    samplename = samplename.split(",") 
    pp2_dir = configuration['parpipe2_dir']
    if not configuration['annot_rank'] in ['', None, False, 'False', 0, '0']: 
        annot_rank_file = configuration['annot_rank']
    else:
        annot_rank_file = os.path.join(pp2_dir, 'annotation', 'annot_categories.txt')
    annotR = os.path.join(pp2_dir, 'bin', 'annot.R')
    main_annot = main_annot_rds
    trna_annot = os.path.join(pp2_dir, 'annotation', 'rds_files','trna.rds')
    repeats_annot = os.path.join(pp2_dir, 'annotation', 'rds_files','repeats.rds')
    annot_list = []
    
    
    for idx, sample in enumerate(samplename):
        infile_clust = output_stats_list[idx]
        annot = os.path.join(outdir, '{0}_{1}.annotation.tsv'.format(sample,ts))
        annot_list.append(annot)
        logfile =  os.path.join(outdir, '{0}_annotation.{1}.log'.format(sample,ts))
        errfile =  os.path.join(outdir, '{0}_annotation.{1}.err'.format(sample,ts))
        
        with open(infile_clust, 'r') as f:
            if sum(1 for line in f) > 1:
                command = ['Rscript', annotR, outdir, annot_rank_file, main_annot, trna_annot, repeats_annot, infile_clust, annot, '>', logfile, '2>', errfile]
        
                print("Running annotation for {0}...".format(sample)) 
                outcode = os.WEXITSTATUS(os.system(" ".join(command)))    
                if outcode == 0:     
                    print("Annotation finished succesfully for {0}".format(sample))
                else:
                    print("an ERROR occured during annotation :(\ncheck {0} for details".format(errfile))
            else:
                command = ['echo No clusters for annotation > ', annot] 
                print("No clusters to annotate for {0}...".format(sample)) 
                os.system(" ".join(command))
            
    return annot_list
    
    
def get_summary_stats(configuration, samplename, reads, cutadapt_log_list, collapsed_reads_list, rm_rep_reads_list, ali_log_list, output_groups_list, output_clust_list, annot_list):

    sum_stats = {"raw_reads": [],
    "reads_with_adapters" : [],
    "too_short_reads" : [],
    "passed_reads":[],
    "collapsed_reads" : [],
    "rm_rep_reads" : [],
    "mapped_reads_total": [],
    'uniquely_mapped_reads' : [],
    'multimapped_reads' : [],
    "unmapped_reads_total" : [],
    "unmapped_reads_too_many_multimaps" : [],
    "unmapped_reads_too_many_mismatches" : [],
    "unmapped_reads_too_short" : [],
    "unmapped_reads_other" : [],
    "groups" : [],
    "clusters": [],
    "annotated_clusters": []
    }
    samplename = samplename.split(",")
    reads = reads.split(",")
    if not configuration["sum_stats_prefix"] in ['', None, False, 'False', '0', 0]:
        out_file = os.path.join(cwd, 'output', '{0}_{1}.summary_stats.tsv'.format(configuration["sum_stats_prefix"], ts))
    else:
        out_file = os.path.join(cwd, 'output', '{0}_{1}.summary_stats.tsv'.format('_'.join(samplename),ts))

    for idx, sample in enumerate(samplename):
        raw_read = reads[idx]
        with gzip.open(raw_read, 'rt') as f:
            raw_read_cnt = sum(1 for line in f)/4
        sum_stats['raw_reads'].append(raw_read_cnt)
        cutadapt_log = cutadapt_log_list[idx]
        with open(cutadapt_log, 'r') as f:
            for line in f:
                spl_line = line.split()
                if line.startswith("Reads with adapters"):
                    sum_stats["reads_with_adapters"].append(int(spl_line[3].replace(",","")))
                elif line.startswith("Reads that were too short"):
                    sum_stats["too_short_reads"].append(int(spl_line[5].replace(",","")))
                elif line.startswith("Reads written (passing filters)"):
                    sum_stats["passed_reads"].append(int(spl_line[4].replace(",","")))
                    
        collapsed_read = collapsed_reads_list[idx]
        with open(collapsed_read, 'r') as f:
            collapsed_read_cnt = sum(1 for line in f)/2
        sum_stats['collapsed_reads'].append(collapsed_read_cnt)
        if not configuration['rem_rep'] in [False, 'False', '0', 0, 'no']:
            rm_rep_read = rm_rep_reads_list[idx]    
            with open(rm_rep_read, 'r') as f:
                rm_rep_read_cnt = sum(1 for line in f)/2
            sum_stats['rm_rep_reads'].append(rm_rep_read_cnt)
        ali_stats = ali_log_list[idx]
        mapped_reads_total = 0
        unmapped_reads_total = 0
        if configuration['aligner'] in ['S', 's', 'STAR', 'star', 'Star']:
            with open(ali_stats, 'r') as f:
                for line in f:
                    spl_line = line.split()
                    if line.lstrip().startswith('Uniquely mapped reads number'):
                        sum_stats['uniquely_mapped_reads'].append(int(spl_line[5]))
                        mapped_reads_total += int(spl_line[5])
                    elif line.lstrip().startswith('Number of reads mapped to multiple loci'):
                        sum_stats['multimapped_reads'].append(int(spl_line[8]))
                        mapped_reads_total += int(spl_line[8])
                    elif line.lstrip().startswith("Number of reads mapped to too many loci"):
                        sum_stats["unmapped_reads_too_many_multimaps"].append(int(spl_line[9]))
                        unmapped_reads_total += int(spl_line[9])
                    elif line.lstrip().startswith("Number of reads unmapped: too many mismatches"):
                        sum_stats["unmapped_reads_too_many_mismatches"].append(int(spl_line[8]))
                        unmapped_reads_total += int(spl_line[8])
                    elif line.lstrip().startswith("Number of reads unmapped: too short"):
                        sum_stats["unmapped_reads_too_short"].append(int(spl_line[7]))
                        unmapped_reads_total += int(spl_line[7])
                    elif line.lstrip().startswith("Number of reads unmapped: other"):
                        sum_stats["unmapped_reads_other"].append(int(spl_line[6]))    
                        unmapped_reads_total += int(spl_line[6])
            sum_stats["mapped_reads_total"].append(mapped_reads_total)
            sum_stats["unmapped_reads_total"].append(unmapped_reads_total)
        else:
            with open(ali_stats, 'r') as f:
                for line in f:
                    spl_line = line.split()
                    if line.startswith('# reads with at least one alignment'):
                        sum_stats['mapped_reads_total'].append(int(spl_line[7]))
                    elif line.startswith("# reads that failed to align"):
                        sum_stats["unmapped_reads_total"].append(int(spl_line[6]))
                    elif line.startswith("# reads with alignments suppressed due to -m"):
                        sum_stats["unmapped_reads_too_many_multimaps"].append(int(spl_line[8]))
        output_group = output_groups_list[idx]
        with open(output_group, 'r') as f:
            output_group_cnt = sum(1 for line in f)-1
        sum_stats['groups'].append(output_group_cnt)
        output_clust = output_clust_list[idx]    
        with open(output_clust, 'r') as f:
            output_clust_cnt = sum(1 for line in f)-1
        sum_stats['clusters'].append(output_clust_cnt)
        annot = annot_list[idx]
        with open(annot, 'r') as f:
            annot_clust_cnt = sum(1 for line in f)-1
        sum_stats['annotated_clusters'].append(annot_clust_cnt)
    
    freq_dict = {
    "raw_reads": [100] * len(samplename),
    "reads_with_adapters": [x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["reads_with_adapters"], sum_stats["raw_reads"])],
    "too_short_reads": [x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["too_short_reads"], sum_stats["raw_reads"])],
    "passed_reads": [x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["passed_reads"], sum_stats["raw_reads"])],
    "collapsed_reads": [x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["collapsed_reads"], sum_stats["raw_reads"])],
    "rm_rep_reads": [x / y * 100 if y != 0 else [] for x, y in zip(sum_stats["rm_rep_reads"], sum_stats["raw_reads"])],
    "mapped_reads_total": [x / y * 100 if all(y != 0 for y in sum_stats["collapsed_reads"]) and not sum_stats["rm_rep_reads"] else x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["mapped_reads_total"], sum_stats["rm_rep_reads"])],
    "uniquely_mapped_reads": [x / y * 100 if all(y != 0 for y in sum_stats["mapped_reads_total"]) and sum_stats["uniquely_mapped_reads"] and sum_stats["mapped_reads_total"] else [] for x, y in zip(sum_stats["uniquely_mapped_reads"], sum_stats["mapped_reads_total"])],
    "multimapped_reads": [x / y * 100 if all(y != 0 for y in sum_stats["mapped_reads_total"]) else [] for x, y in zip(sum_stats["multimapped_reads"], sum_stats["mapped_reads_total"])],
    "unmapped_reads_total": [x / y * 100 if all(y != 0 for y in sum_stats["collapsed_reads"]) and not sum_stats["rm_rep_reads"] else x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["unmapped_reads_total"], sum_stats["rm_rep_reads"])],
    "unmapped_reads_too_many_multimaps": [x / y * 100 if y != 0 else 0 for x, y in zip(sum_stats["unmapped_reads_too_many_multimaps"], sum_stats["unmapped_reads_total"])],
    "unmapped_reads_too_many_mismatches": [x / y * 100 if all(y != 0 for y in sum_stats["unmapped_reads_total"]) else [] for x, y in zip(sum_stats["unmapped_reads_too_many_mismatches"], sum_stats["unmapped_reads_total"])],
    "unmapped_reads_too_short": [x / y * 100 if all(y != 0 for y in sum_stats["unmapped_reads_total"]) else [] for x, y in zip(sum_stats["unmapped_reads_too_short"], sum_stats["unmapped_reads_total"])],
    "unmapped_reads_other": [x / y * 100 if all(y != 0 for y in sum_stats["unmapped_reads_total"]) else [] for x, y in zip(sum_stats["unmapped_reads_other"], sum_stats["unmapped_reads_total"])],
    "groups": [x / y * 100 if all(y != 0 for y in sum_stats["mapped_reads_total"]) else 0 for x, y in zip(sum_stats["groups"], sum_stats["mapped_reads_total"])],
    "clusters": [x / y * 100 if all(y != 0 for y in sum_stats["groups"]) else 0 for x, y in zip(sum_stats["clusters"], sum_stats["groups"])],
    "annotated_clusters": [x / y * 100 if all(y != 0 for y in sum_stats["clusters"]) else [] for x, y in zip(sum_stats["annotated_clusters"], sum_stats["clusters"])]
    }
        
    if configuration['rem_rep'] in [False, 'False', '0', 0, 'no']:
        freq_dict["mapped_reads_total"] = [x / y * 100 if all(y != 0 for y in sum_stats["collapsed_reads"]) else 0 for x, y in zip(sum_stats["mapped_reads_total"], sum_stats["collapsed_reads"])]
        freq_dict["unmapped_reads_total"] = [x / y * 100 if all(y != 0 for y in sum_stats["collapsed_reads"]) else 0 for x, y in zip(sum_stats["unmapped_reads_total"], sum_stats["collapsed_reads"])]
        

    keys_to_delete = []
    for x in sum_stats.keys():
        if sum_stats[x] == []:
            keys_to_delete.append(x)

    for key in keys_to_delete:
        del sum_stats[key]
        del freq_dict[key]
    
    sum_stats_tab = pd.DataFrame.from_dict(sum_stats, orient = 'index', columns = samplename)
    freq_tab = pd.DataFrame.from_dict(freq_dict, orient = 'index', columns = ['{0}_freq'.format(sample) for sample in samplename])

    final_df = pd.merge(sum_stats_tab, freq_tab, left_index=True, right_index=True)
    if final_df.shape[0] == 11:
        analysis_step = ['reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'alignment', 'alignment', 'alignment', 'clusters_calling', 'clusters_calling', 'clusters_calling']
        denominator = ['raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'collapsed_reads', 'collapsed_reads', 'unmapped_reads_total', 'mapped_reads_total', 'groups', 'clusters']
    elif final_df.shape[0] == 12:
        analysis_step = ['reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'alignment', 'alignment', 'alignment', 'clusters_calling', 'clusters_calling', 'clusters_calling']
        denominator = ['raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'rm_rep_reads', 'rm_rep_reads', 'unmapped_reads_total', 'mapped_reads_total', 'groups', 'clusters']
    elif final_df.shape[0] == 16:
        analysis_step = ['reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'clusters_calling', 'clusters_calling', 'clusters_calling']
        denominator = ['raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'collampsed_reads', 'mapped_reads_total', 'mapped_reads_total', 'collapsed_reads', 'unmapped_reads_total', 'unmapped_reads_total', 'unmapped_reads_total', 'unmapped_reads_total', 'mapped_reads_total', 'groups', 'clusters']
    elif final_df.shape[0] == 17:
        analysis_step = ['reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'reads_processing', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'alignment', 'clusters_calling', 'clusters_calling', 'clusters_calling']
        denominator = ['raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'raw_reads', 'rm_rep_reads', 'mapped_reads_total', 'mapped_reads_total', 'rm_reps_reads', 'unmapped_reads_total', 'unmapped_reads_total', 'unmapped_reads_total', 'unmapped_reads_total', 'mapped_reads_total', 'groups', 'clusters']
    
    final_df["analysis_step"] = analysis_step
    final_df["denominator"] = denominator
    final_df.to_csv(out_file, sep = "\t", index=True)
    print(final_df)
    
    
    return out_file   
    
def prepare_plots(configuration, sum_stats):
    
    sum_stats = sum_stats
    outdir = os.path.join(cwd, 'output')
    pp2_dir = configuration['parpipe2_dir']
    generate_plots_R = os.path.join(pp2_dir, 'bin', 'generate_plots.R')
    logfile =  os.path.join(outdir, 'summary_plots.{0}.log'.format(ts))
    errfile =  os.path.join(outdir, 'summary_plots.{0}.err'.format(ts))
    
    command = ['Rscript', generate_plots_R, sum_stats, outdir, ts, '>', logfile, '2>', errfile]
    
    print("\nPreparing summary plots...")
    outcode = os.WEXITSTATUS(os.system(" ".join(command)))    
    if outcode == 0:     
        print("Summary plots prepared")
    else:
        print("an ERROR occured during summary plots preparation :(\ncheck {0} for details".format(errfile))

def main(args):
    config_file = os.path.abspath(args.config_file)
    cfg = read_configfile(config_file)
    threads = cfg['threads']
    aligner = cfg['aligner']
    reads = os.path.abspath(args.reads)    
    samplename = args.samplename
    if args.umi:
        umi = args.umi
    elif not args.umi:
        umi = None
        
    genome, annotation, genome_2bit, main_annot_rds  = prepare_genome_and_annotation(cfg)
    if aligner in ['S', 's', 'STAR', 'star', 'Star']:
        star_idx_outdir = prepare_star_index(cfg, genome, annotation, threads)
    else :
        bowtie_idx_outdir = prepare_bowtie_index(cfg, genome, threads)
    trim_reads1, cutadapt_log = trim_adapters_round1(cfg, reads, samplename)
    trim_reads2 = trim_adapters_round2(cfg, samplename, trim_reads1)
    collapsed_reads = collapse_reads(samplename, trim_reads2, umi)
    rm_rep_reads = star_rm_reps(cfg, samplename, collapsed_reads, threads)
    if aligner in ['S', 's', 'STAR', 'star', 'Star']:
        out_sam, ali_logs = star_genome_map(cfg, samplename, star_idx_outdir, rm_rep_reads, threads)
    else:
        out_sam, ali_logs = bowtie_genome_map(cfg, samplename, bowtie_idx_outdir, rm_rep_reads, threads)
    #merged_sam, merged_bed, genome_annot = annotate_mapped_reads(cfg, samplename, out_sam, main_annot_rds)
    #genome_conv_stats(samplename, merged_sam, merged_bed, genome_annot, genome)
    out_groups, out_clust, = config_and_run_PARalyzer(cfg, samplename, out_sam, genome_2bit)
    out_stats = run_clusters_stats(samplename, out_sam, out_clust, genome)
    annot = run_annot(cfg, samplename, out_stats, main_annot_rds)
    sum_stats = get_summary_stats(cfg, samplename, reads, cutadapt_log, collapsed_reads, rm_rep_reads, ali_logs, out_groups, out_clust, annot)
    prepare_plots(cfg, sum_stats)
    
    print('\nPAR-PIPE2 FINISHED REQUESTED ANALYSIS SUCCESFULLY FOR JOB {0}. CHECK YOUR OUTPUT FILES'.format(ts))
    
   
               
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Analyze PAR-CLIP')
    parser.add_argument('--config_file', type = str, help = 'Yaml file with cutadapt and STAR parameters', required = True)
    parser.add_argument('--reads', type = str, help = 'Reads containing RNA sequence (dependnent on library protocol may be either read1 or read2, if > 1, comma separated)). Gzipped fastq.', required = True)
    parser.add_argument('--samplename', type = str, help = 'Sample name (if > 1, comma separated). Will be appended to output files.', required = True)
    parser.add_argument('--umi', type = int, help = "UMI location, if 5' UMI - positive number, 3' UMI negative number,  default None", required = False)
    args = parser.parse_args()   
    main(args)

import os
import sys
import wget
import time
import yaml
import json
import argparse

cwd = os.getcwd()

def read_configfile(config_file):

    print("Loading config file..." )

    with open(config_file) as f_cfg:
        cfg = yaml.load(f_cfg, Loader=yaml.FullLoader)
    print("Config file succesfully loaded")
   
    return cfg
    
def prepare_genome_and_annotation(configuration):
    cwd = os.getcwd()
    if not os.path.exists('genome'):
        os.mkdir('genome')
    outdir = os.path.join(cwd, 'genome')
    if not os.path.exists('logs/genome_files'):
        os.makedirs('logs/genome_files')
    logdir = os.path.join(cwd, 'logs', 'genome_files')
  
    if configuration["genome_fasta"] and configuration["genome_fasta"] != "":
        genome = configuration["genome_fasta"]
        annotation = configuration["gtf"]
        print("Genome and annotation files are provided by user")
    
    else:

        outdir = os.path.join(cwd, 'genome')
        genome_file = "GRCh38.primary_assembly.genome.fa"
        genome_ver = configuration.get('genome_version', '45')
        if not genome_ver:
            genome_ver = '45'
        if int(genome_ver) > 41:
            gtf_file = f"gencode.v{genome_ver}.primary_assembly.basic.annotation.gtf"
        else:
            gtf_file = f"gencode.v{genome_ver}.primary_assembly.annotation.gtf"
        
        print(f"Selected genome version is {genome_ver}")
        genome = os.path.join(outdir, genome_file)
        annotation = os.path.join(outdir, gtf_file)
        
        if not os.path.isfile(annotation) or not os.path.isfile(genome):
            genome_gz = os.path.join(outdir, f"{genome_file}.gz")
            gtf_gz = os.path.join(outdir, f"{gtf_file}.gz")
            logfile = os.path.join(logdir, f'genome_and_gtf_unpacking.log')
            errfile = os.path.join(logdir, f'genome_and_gtf_unpacking.err')
            url_genome = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{genome_ver}/{genome_file}.gz"
            url_gtf = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{genome_ver}/{gtf_file}.gz"
            print(f"Downloading GRCh38 v.{genome_ver} genome from Gencode...\n")
            wget.download(url_genome, out=outdir)
            print("\nDownloading genome annotation file...\n")
            wget.download(url_gtf, out=outdir)
            print("\nUnpacking genome and gtf files...")
            gunzip_genome_cmd = f"gunzip {genome_gz} > {logfile} 2> {errfile}"
            gunzip_gtf_cmd = f"gunzip {gtf_gz} >> {logfile} 2>> {errfile}"
            outcode_genome = os.WEXITSTATUS(os.system(gunzip_genome_cmd))
            outcode_gtf = os.WEXITSTATUS(os.system(gunzip_gtf_cmd))
            if outcode_genome == 0:     
                print("Genome unpacked successfully")
            else:
                print(f"An ERROR occurred during genome unpacking :( check {errfile} for details")
            if outcode_gtf == 0:     
                print("Gtf unpacked successfully")
            else:
                print(f"An ERROR occurred during gtf file unpacking :( check {errfile} for details")
        else:
            print('Genome and annotation files already exist')
                      
    pp2_path = configuration['parpipe2_dir']
    main_annot_rds = os.path.join(outdir, 'main_annotation.rds')
    
    if os.path.isfile(main_annot_rds) and os.path.getmtime(main_annot_rds) > os.path.getmtime(annotation):
        print('main_annotation.rds file from current gtf already exists')
    else:
        print("Preparing main_annotation.rds file")
        logfile = os.path.join(logdir, f'prep_main_rds_from_gtf.log')
        errfile = os.path.join(logdir, f'prep_main_rds_from_gtf.err')
        ann_script = os.path.join(pp2_path, "bin", "prepareAnnotFilesFromGtf.R")
        ann_cmd = f"Rscript {ann_script} {outdir} {annotation} main > {logfile} 2> {errfile}"
        outcode = os.WEXITSTATUS(os.system(ann_cmd))
    
        if outcode == 0:     
            print(f"main_annotation.rds file created successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during creating main_annotation.rds file :( check {errfile} for details")
    
    genome_2bit = os.path.normpath(str(genome).replace('fasta', '2bit').replace('fa', '2bit'))
    
    if os.path.isfile(genome_2bit) and os.path.isfile(f'{genome}.fai'):
        print('Index and 2bit files exist. Moving to the next step')
    else:
        logfile = os.path.join(logdir, f'genome_indexing_and_2bit_preparation.log')
        errfile = os.path.join(logdir, f'genome_indexing_and_2bit_preparation.err')
        index_genome_cmd = f"samtools faidx {genome} > {logfile} 2> {errfile}"
        command_2bit = f"faToTwoBit {genome} {genome_2bit} >> {logfile} 2>> {errfile}"
        outcode = os.WEXITSTATUS(os.system(index_genome_cmd))
        if outcode == 0:     
            print(f"Genome file indexed successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during genome indexing :( check {errfile} for details")
        print("Preparing genome 2bit file...")
        outcode = os.WEXITSTATUS(os.system(command_2bit))
        if outcode == 0:     
            print(f"Genome 2bit file prepared successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during genome 2bit preparation :( check {errfile} for details")

    # Write the paths to a metadata file
    metadata = {
        "genome": genome,
        "annotation": annotation,
        "genome_2bit": genome_2bit,
        "main_annot_rds": main_annot_rds
    }
    with open(os.path.join("metadata.json"), "w") as metafile:
        json.dump(metadata, metafile, indent=4)
    
    return genome, annotation, genome_2bit, main_annot_rds


def prepare_star_index(configuration, genome, annotation, threads):

    if not os.path.exists('logs/star_idx'):
        os.makedirs('logs/star_idx')
    logdir = os.path.join(cwd, 'logs', 'star_idx')

    if not configuration.get('star_index_dir'):
        if not os.path.exists('star_index'):
            os.mkdir('star_index')
        star_idx_outdir = os.path.join(cwd, 'star_index')
        logfile = os.path.join(logdir, f'star_index.log')
        errfile = os.path.join(logdir, f'star_index.err')
        if not configuration.get('star_index_params'):
            command = f"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {star_idx_outdir} --genomeFastaFiles {genome} --sjdbGTFfile {annotation} --sjdbOverhang 100 > {logfile} 2> {errfile}"
        else:
            params = " ".join(configuration['star_index_params'].split(","))
            command = f"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {star_idx_outdir} --genomeFastaFiles {genome} --sjdbGTFfile {annotation} {params} > {logfile} 2> {errfile}"
        
        print("Generating STAR index...")
        outcode = os.WEXITSTATUS(os.system(command))
        if outcode == 0:     
            print(f"STAR index generated successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during STAR index generation :(\ncheck {errfile} for details")
    else:
        print("STAR index provided by user")
        star_idx_outdir = configuration['star_index_dir']
    
    with open("metadata.json", "r") as metafile:
        metadata = json.load(metafile)
    
    metadata['star_index'] = star_idx_outdir
    
    with open("metadata.json", "w") as metafile:
        json.dump(metadata, metafile, indent=4)
    
    return star_idx_outdir

def prepare_bowtie_index(configuration, genome, threads):

    if not os.path.exists('logs/bowtie_idx'):
        os.makedirs('logs/bowtie_idx')
    logdir = os.path.join(cwd, 'logs', 'bowtie_idx')

    if not configuration.get('bowtie_index_dir'):
        if not os.path.exists('bowtie_index'):
            os.mkdir('bowtie_index')
        bowtie_idx_outdir = os.path.join(cwd, 'bowtie_index')
        logfile = os.path.join(logdir, f'bowtie_index.log')
        errfile = os.path.join(logdir, f'bowtie_index.err')
        if not configuration.get('bowtie_index_params'):
            command = f"bowtie-build --threads {threads} -f {genome} {bowtie_idx_outdir}/bwt > {logfile} 2> {errfile}"
        else:
            params = " ".join(configuration['bowtie_index_params'].split(","))
            command = f"bowtie-build {params} -f {genome} {bowtie_idx_outdir}/bwt > {logfile} 2> {errfile}"
        
        print("Generating bowtie index...")
        outcode = os.WEXITSTATUS(os.system(command))
        if outcode == 0:     
            print(f"Bowtie index generated successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during bowtie index generation :(\ncheck {errfile} for details")
    else:
        print("Bowtie index provided by user")
        bowtie_idx_outdir = configuration['bowtie_index_dir']
    
    with open("metadata.json", "r") as metafile:
        metadata = json.load(metafile)
    
    metadata['bowtie_index'] = bowtie_idx_outdir
    
    with open("metadata.json", "w") as metafile:
        json.dump(metadata, metafile, indent=4)
    
    return bowtie_idx_outdir


def main(args):
    config_file = os.path.abspath(args.config_file)
    cfg = read_configfile(config_file)
    threads = cfg['threads']
    aligner = cfg['aligner']
    genome, annotation, genome_2bit, main_annot_rds  = prepare_genome_and_annotation(cfg)
    if aligner in ['S', 's', 'STAR', 'star', 'Star']:
        prepare_star_index(cfg, genome, annotation, threads)
    else :
        prepare_bowtie_index(cfg, genome, threads)
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Analyze PAR-CLIP')
    parser.add_argument('--config_file', type = str, help = 'Yaml file with required parameters', required = True)
    args = parser.parse_args()   
    main(args)

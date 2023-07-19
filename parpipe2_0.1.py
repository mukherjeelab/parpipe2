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
        command = cfg['cutadapt_params1'].split(",") +  ["-o", output, read]
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
            command = cfg['cutadapt_params2'].split(',') +  ['-o', output, infile]
            print("Running cutadapt round 2 for {0}...".format(sample))
            print(command)
            
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

#TODO                    
#star_rm_reps
#star_genome_map 
#config_and_run_PARalyzer  
#run perbase
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

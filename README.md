# PAR-PIPE2
PAR-PIPE2 is the python wrapper to analyze PAR-CLIP data. 

## Overview

[Introduction](#introduction)

[Prerequisites](#prerequisites)

[Dependencies](#dependencies)

[Preparing configuration file](#preparing_configuration_file)

[Usage](#usage)

[Contributors](#contributors)

[License](#license)

## Introduction
PAR-PIPE2 runs whole analysis starting from gzipped fastq files end ending with annotation tables using single command and yaml configuration file. It is user friendly and highly customizable tool.

## Prerequisites
Before PAR-PIPE2 run you have to prepare some additional files.
1. Indexed genome sequence in FASTA format. We reccomend primary assembly files from GENCODE. To index genome simply run:
```
samtools faidx genome_file.fa
```
2. Genome sequence in .2bit format.
3. [STAR](https://github.com/alexdobin/STAR) index for your genome (optimally generated using STAR 2.7.10b). Star manual can be found [here](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf).
4. STAR index for repetitive elements (optional).
5. Gtf annotation file for your genome assembly.
6. Gtf annotation file for tRNAs in your genome version (optional).
7. Gtf annotation file for repetitive elemants in the genome.

## Dependencies
Required libraries are stored in the `environment.yml` file. All required dependencies can be installed using *conda* by executing the following command:
```
$ conda env create -f environment.yml
```

The environment then needs to be activated in order to run PAR-PIPE2:
```
$ conda activate parpipe2
```
We also plan to provide Docker image soon.

## Preparing configuration file
Be described soon.

## Usage
First clone the repository using: 
```
$ gh repo clone mukherjeelab/parpipe2
```
or:
```
$ github clone mukherjeelab/parpipe2
```
Second, create directory for your project, e.g.:
```
$ mkdir my_parclip_dir
```
Then you have to prepare annotation files. There are three types of annotation files: **main** - which is created from standard genome annotation in gtf format, **trna** - which contains tRNA coordinates and **rep**, which contains repetitive elemants coordinates. Main annotation file is required, remaining two are optional. To prepare these files for PAR-PIPE2 simply run:
```
$ Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/main_annotation.gtf main
$ Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/trna.gtf trna
$ Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/rep_elements.gtf rep
```
You can also run:
```
$ Rscript prepareAnnotFilesFromGtf.R -h
```
or:
```
$ Rscript prepareAnnotFilesFromGtf.R --help
```
to print help message.
Now you can go to your project directory:
```
$ cd my_parclip_dir
```
**NOTE: PAR-PIPE2 must be run from your project directory.**

And run PAR-PIPE2:
```
$ python /path/to/parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/reads.fastq.gz --samplename prefix_added_to_output_files
```
NOTE: If you want to run PAR-PIPE2 for more than one sample, specify more than one reads (comma separated, no spaces) and samplename (comma separated, no spaces). Number of argumants in --reads and --samplename have to be the same, e.g.:
```
$ python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B
```
NOTE: If your reads contain UMIs add _--umi_ followed by positive or negative number to your command. Positive number means that UMIs are located at 5' end of the read, negative - 3' end (in both cases after removing the adapters), e.g.:
```
$ python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B --umi 6 # 6 first nucleotides are UMIs
$ python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B --umi -4 # 4 last nucleotides are UMIs
```
After running PAR-PIPE2 your project folder should have the following structure:
```
my_parclip_dir
│
└───annotation_files
|      | main_annotation.rds
|      | repeats.rds
|      | trna.rds
│
└───reads
│      │ A_trimmed1.fq.gz
│      │ A_trimmed_final.fq.gz
│      │ A_trimmed.fa
│      │ A_trimmed_collapsed.fa
│      │ A_rmreps.fa
│      │ B_trimmed1.fq.gz
│      │ B_trimmed_final.fq.gz
│      │ B_trimmed.fa
│      │ B_trimmed_collapsed.fa
│      │ B_rmreps.fa
│
└───ali
|      | A.reps.bam
|      | A.aligned.bam
|      | A.aligned.mapped.sam
|      | A.aligned.mapped.sorted.sam
|      | B.reps.bam
|      | B.aligned.bam
|      | B.aligned.mapped.sam
|      | B.aligned.mapped.sorted.sam
│
└───PARalyzer
|      | A.ini
|      | A_PARalyzer_Utilized.sam
|      | A_PARalyzer_Utilized.sorted.bam
|      | A_PARalyzer_Utilized.sorted.bam.bai
|      | A.distribution
|      | A.groups
|      | A.clusters
|      | A_clusters_conv_stats.tsv
|      | A.ini
|      | B_PARalyzer_Utilized.sam
|      | B_PARalyzer_Utilized.sorted.bam
|      | B_PARalyzer_Utilized.sorted.bam.bai
|      | B.distribution
|      | B.groups
|      | B.clusters
|      | B_clusters_conv_stats.tsv
│
└───output
|      | A_annotation.tsv
|      | B_annotation.tsv
```
NOTE: Directories will also contain log, error and some other temporary files not listed above.

## Contributors
* Marcin Sajek
* Neelanjan Mukherjee
* Samantha Lisy
* Manuel Ascano Jr

## License
GNU GPL license

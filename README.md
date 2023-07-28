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
* Indexed genome sequence in FASTA format. We reccomend primary assembly files from [GENCODE](https://www.gencodegenes.org/). Samtools manual can be found [here](http://www.htslib.org/doc/samtools.html). To index genome simply run:
```
samtools faidx genome_file.fa
```
* Genome sequence in .2bit format. You can use faToTwoBit file from this repo:
```
faToTwoBit genome.fa genome.2bit
```
* [STAR](https://github.com/alexdobin/STAR) index for your genome (optimally generated using STAR 2.7.10b). Star manual can be found [here](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf).
* [STAR](https://github.com/alexdobin/STAR) index for repetitive elements (optimally generated using STAR 2.7.10b) (optional).
* Gtf annotation file for your genome assembly, available in [GENCODE](https://www.gencodegenes.org/).
* Gtf annotation file for tRNAs in your genome version, available in [GENCODE](https://www.gencodegenes.org/) (optional).
* Gtf annotation file for repetitive elemants in the genome (optional).
* Annotation rank file.
Annotation rank file contain ranked list of features from your gtf(s) files. When PAR-CLIP peak is mapping to more than one genomic feature, top feature from annotation rank file is selected to annotate this peak. ```annot_categories.txt``` is an example of annotation rank file.

## Dependencies
Required libraries are stored in the `environment.yml` file. All required dependencies can be installed using *conda* by executing the following command:
```
conda env create -f environment.yml
```

The environment then needs to be activated in order to run PAR-PIPE2:
```
conda activate parpipe2
```
We also plan to provide Docker image soon.

## Preparing configuration file
Most of the parameters for PAR-PIPE2 run have to be set up in ```config.yaml`` file. The configuration file contains five main sections:
* general parameters
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) parameters
* [STAR](https://github.com/alexdobin/STAR) parameters
* [PARalyzer](https://ohlerlab.mdc-berlin.de/files/duke/PARalyzer/README.txt) parameters
* annotation parameters
General parameters contain two fields: ```parpipe2_dir``` and ```threads```. The first one contains absolute path to your PAR-PIPE2 directory (cloned from this repo). The second, number of threads that should be used for the analysis. Example below:
```
{parpipe2_dir : "/media/marcin/data2/par_pipe2_test/",
threads : 10,
```
Cutadapt section contain adapter trimming parameters. Please familiarize yourself with [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) before filling this section. For some analyses it may be better to split adapter trimming for two rounds. If you want to do it at one step fill ```cutadapt_params1``` field and leave ```cutadapt_params2``` field with empty quotes ```cutadapt_params2 : ""```. To fill these fields properly provide a cutadapt parameter name and its value separated by coma withount any spaces. If there is a boolean parameter (without any value passed to it), put it in a file without any additional values. Do not put -o (--output) and input files in ```cutadapt_params1``` and ```cutadapt_params2```. They will be determined automatically by PAR-PIPE2 based on provided sample names. An example of filled ```cutadapt_params1``` and ```cutadapt_params2``` fields is showed below:
```
#cutadapt parameters
cutadapt_params1 : "-g,GTTCAGAGTTCTACAGTCCGACGATC,-b,CGTACGCGGGTTTAAACGA,-b,CGTACGCGGAATAGTTTAAACTGT,-n,3,-j,0,--nextseq-trim=15",
cutadapt_params2 : "-a,TGGAATTCTCG,-n,3,-j,0,-m,23,-M,50,--discard-untrimmed,--max-n,1",
```
TIP: ```-j 0 ``` in cutadapt determined automatically numbers of available threads (should be written in config file as ```,-j,0,```.
TIP: If your sequencing was performed on Novaseq platform add ```--nextseq-trim``` parameter, which automatically remove empty runs.
TIP: To long reads can cause further issues with PARalyzer. Consider setting max read length parameter, e.g. ```-M,50```.

The next section of configuration file contain parameters list passed to  STAR. Detailed description of  STAR parameters may be find in this [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf). STAR alignment is divided into two steps. First, optional step is mapping your reads to repetitive elements in the genome. To perfom this step you have to fill ```rep_idx``` field with absolute path to the STAR index for repetitive elements and ```star_rem_reps_params``` with STAR parameters that will be passed to PAR-PIPE2. If you leave ```star_rem_reps_params``` blank ```star_rem_reps_params : ""``` this step will be omitted. At the next step reads that were not aligned to repetitive elments (or all reads if first step was omitted) are mapped to the genome. ```genome_idx``` required absolute path to the STAR index for your genome and ```star_map_params``` list of STAR parameters. ```star_rem_reps_params``` and ``star_map_params``` must be filled in exactly the same way as ```cutadapt_params1``` and ```cutadapt_params2```. Do not include ```--genomeDir``` parameter in ```star_rem_reps_params``` and ``star_map_params```. It will be determined based on ```rep_idx``` and|or ```genome_idx``` fields. Also do not include ```--readFilesIn``` and ```--outFileNamePrefix``` parameters. They will be determined automatically by PAR-PIPE2. The example of filled STAR parameters section is showed below:
```
#STAR parameters
rep_idx : "/media/marcin/data2/par_pipe2_test/star_test/rep_idx",
star_rem_reps_params : "--runMode,alignReads,--genomeLoad,NoSharedMemory,--alignEndsType,EndToEnd,--outSAMunmapped,Within,--outFilterMultimapNmax,30,--outFilterMultimapScoreRange,1,--outSAMtype,BAM,Unsorted,--outFilterType,BySJout,--outBAMcompression,10,--outReadsUnmapped,Fastx,--outFilterScoreMin,10,--outSAMattrRGline,ID:foo,--outSAMattributes,All,--outSAMmode,Full",
genome_idx : "/media/marcin/data2/par_pipe2_test/star_test/genome_idx",
star_map_params : "--runMode,alignReads,--genomeLoad,NoSharedMemory,--alignEndsType,EndToEnd,--outSAMunmapped,Within,--outFilterMultimapNmax,1,--outFilterMultimapScoreRange,1,--outSAMattributes,All,--outSAMtype,BAM,Unsorted,--outFilterType,BySJout,--outFilterMismatchNoverReadLmax,0.2,--outSAMattrRGline,ID:foo,--outStd,Log,--outBAMcompression,10,--outSAMmode,Full",
```
TIP: If parmeter takes more than one value, simply separate both values by coma, e.g. ```--outSAMtype BAM Unsorted``` should be passed to config file as ```--outSAMtype,BAM,Unsorted```.

PARalyzer parameters section is directly use to create .ini file for PARalyzer run and setup memory limit for PARalyzer (```paralyzer_memory``` field). To properly fill all fields please read PARalyzer [manual](https://ohlerlab.mdc-berlin.de/files/duke/PARalyzer/README.txt). Field ```ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL``` is required only if ```signal_extension : ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL```. In other cases will be skipped. Example of PARalyzer parameters section below:
```
#PARalyzer parameters
paralyzer_memory : "4G",
BANDWIDTH : 3,
CONVERSION : "T>C",
MINIMUM_READ_COUNT_PER_GROUP : 5,
MINIMUM_READ_COUNT_PER_CLUSTER : 2,
MINIMUM_READ_COUNT_FOR_KDE : 3,
MINIMUM_CLUSTER_SIZE : 11,
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER : 2,
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER : 2,
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION : 1,
MINIMUM_READ_LENGTH : 15,
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES : 1,
signal_extension : "EXTEND_BY_READ", #one of 'EXTEND_BY_READ', 'HAFNER_APPROACH', 'ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL'
ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL : 2, #only valid when signal_extension : 'ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL'
GENOME_2BIT_FILE: "/media/marcin/data2/par_pipe2_test/ZFP36/PARalyzer/GRCh38_no_alt_analysis_set_GCA_000001405.15.2bit",
```
Annotation section contains ```genome_fasta``` field - absolute path to your genome sequence and ```annot_rank``` field - absolute path to your annotation rank file. Example below:
```
#ANNOTATION
genome_fasta : "/media/marcin/data2/par_pipe2_test/paralyzer_test/PARalyzer/GRCh38_no_alt_analysis_set_GCA_000001405.15.fa",
annot_rank : "/media/marcin/data2/par_pipe2_test/annot_files/annot_categories.txt"}
```

## Usage
First clone the repository using: 
```
gh repo clone mukherjeelab/parpipe2
```
or:
```
github clone mukherjeelab/parpipe2
```
Second, create directory for your project, e.g.:
```
mkdir my_parclip_dir
```
Then you have to prepare annotation files. There are three types of annotation files: **main** - which is created from standard genome annotation in gtf format, **trna** - which contains tRNA coordinates and **rep**, which contains repetitive elemants coordinates. Main annotation file is required, remaining two are optional. To prepare these files for PAR-PIPE2 simply run:
```
Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/main_annotation.gtf main
Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/trna.gtf trna
Rscript prepareAnnotFilesFromGtf.R /absolute/path/to/my_parclip_dir /absolute/path/to/rep_elements.gtf rep
```
You can also run:
```
Rscript prepareAnnotFilesFromGtf.R -h
```
or:
```
$ Rscript prepareAnnotFilesFromGtf.R --help
```
to print help message.
Now you can go to your project directory:
```
cd my_parclip_dir
```
**NOTE: PAR-PIPE2 must be run from your project directory.**

And run PAR-PIPE2:
```
python /path/to/parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/reads.fastq.gz --samplename prefix_added_to_output_files
```
NOTE: If you want to run PAR-PIPE2 for more than one sample, specify more than one reads (comma separated, no spaces) and samplename (comma separated, no spaces). Number of argumants in --reads and --samplename have to be the same, e.g.:
```
python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B
```
NOTE: If your reads contain UMIs add _--umi_ followed by positive or negative number to your command. Positive number means that UMIs are located at 5' end of the read, negative - 3' end (in both cases after removing the adapters), e.g.:
```
python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B --umi 6 # 6 first nucleotides are UMIs
python parpipe2_0.1.py --config_file /path/to/yaml/configuration/file --reads /path/to/A.fastq.gz,/path/to/B.fastq.gz --samplename A,B --umi -4 # 4 last nucleotides are UMIs
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

# PAR-PIPE2
PAR-PIPE2 is the python wrapper to analyze PAR-CLIP data. 

## Overview

[Introduction](#introduction)

[Dependencies](#dependencies)

[Usage](#usage)

[Contributors](#contributors)

[License](#license)

## Introduction
PAR-PIPE2 runs whole analysis starting from gzipped fastq files end ending with annotation tables using single command and yaml configuration file. It is user friendly and highly customizable tool.

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

## Usage
First clone the repository using: 
```
$ gh repo clone mukherjeelab/parpipe2
```
or:
```
$ github clone mukherjeelab/parpipe2
```
Second, create directory for your project, eg:
```
$ mkdir my_parclip_dir
```
Then you have to prepare annotation files. There are three types of annotation files: **main** - wchich is created from standard genome annotation in gtf format, **trna** - which contains tRNA coordinates and **rep**, which contains repetitive elemants coordinates. Main annotation file is required, remaining two are optional. To prepare these files for PAR-PIPE2 simply run:
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
to get help message
## Contributors
* Marcin Sajek
* Neelanjan Mukherjee
* Samantha Lisy
* Manuel Ascano Jr

## License
GNU GPL license

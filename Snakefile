import os
import sys
import json

def get_input(wildcards):
    ext = sample_names[wildcards.sample]
    return f"reads/{wildcards.sample}.{ext}"
    
def get_index_path(aligner):
    with open("metadata.json") as f:
        data = json.load(f)
    return data[f"{aligner}_index"]

def get_2bit_path():
    with open("metadata.json") as f:
        data = json.load(f)
    return data["genome_2bit"]

def get_genome_path():
    with open("metadata.json") as f:
        data = json.load(f)
    return data["genome"]

# Ensure required config parameters are provided
if not config.get("input_files"):
    raise ValueError("Config parameter 'input_files' must be provided.")

# Extract input file names from the config
config_file_path = None
if "--configfile" in sys.argv:
    # Get the index of the --configfile argument
    configfile_index = sys.argv.index("--configfile")
    
    # Retrieve the path to the config file from the next argument
    config_file_path = sys.argv[configfile_index + 1]
    
    print("Config file path:", config_file_path)
else:
    raise ValueError("No config file path specified via command line arguments.")

input_filenames = config["input_files"].split(",")
num_cores = int(config.get("threads", 1))

# List of possible extensions
sample_names = {}
extensions = ['fq.gz', 'fq', 'fastq.gz', 'fastq']
for filename in input_filenames:
    sample_name = os.path.basename(filename).split(".")[0]
    for ext in extensions:
        if filename.endswith(ext):
            sample_names[sample_name] = ext
            break
print(sample_names)

if config.get('umi', '') == '' or int(config['umi']) == 0:
    umi_param = ""
elif int(config['umi']) > 0:
    umi_param = f"| fastx_trimmer -f {int(config['umi']) + 1}"
else:
    umi_param = f"| fastx_trimmer -t {int(config['umi'])}"

ruleorder: create_symlinks > trim_adapters_r1

rule all:
    input:
        expand("qc/fastqc/{sample}_fastqc.html", sample=sample_names.keys()),
        expand("qc/fastqc/{sample}_fastqc.zip", sample=sample_names.keys()),
        "metadata.json",
        #expand("reads/{sample}_trimmed1.fq.gz", sample=sample_names.keys()),
        #expand("reads/{sample}_trimmed_final.fq.gz", sample=sample_names.keys())
        #expand("reads/{sample}_trimmed_collapsed.fa", sample=sample_names.keys())
        #expand("reads/{sample}_rm_rep.fa", sample=sample_names.keys())
        #expand("reads/{sample}_rm_rep.fa", sample=sample_names.keys())
        #expand("ali/{sample}.aligned.mapped.sorted.sam", sample=sample_names.keys())
        #expand("PARalyzer/{sample}.ini", sample=sample_names.keys())
        #expand("PARalyzer/{sample}_PARalyzer_Utilized.bam.bai", sample=sample_names.keys())
        expand("stats/{sample}_clusters_conv_stats.tsv", sample=sample_names.keys())

rule create_symlinks:
    input:
        input_filenames
    output:
        expand("reads/{sample}.{ext}", zip, sample=sample_names.keys(), ext=sample_names.values())
    run:
        # Create the 'reads' directory if it doesn't exist
        if not os.path.exists("reads"):
            print("Creating reads directory")
            os.makedirs("reads")
            
        for filename in input_filenames:
            sample_name = os.path.basename(filename).split(".")[0]
            ext = sample_names[sample_name]
            symlink_path = os.path.join("reads", f"{sample_name}.{ext}")
            if not os.path.lexists(symlink_path):
                os.symlink(filename, symlink_path)

rule fastqc:
    input:
        lambda wildcards: f"reads/{wildcards.sample}.{sample_names[wildcards.sample]}"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    params:
        num_cores=num_cores
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p qc/fastqc
        echo "Running FastQC on {input}" > {log}
        fastqc --quiet -t {params.num_cores} {input} --outdir qc/fastqc >> {log} 2>&1
        """

rule get_genomic_files_and_index:
    input:
        script=config["parpipe2_dir"] + "/bin/genomic_files_and_indexes.py"
    output:
        "metadata.json"
    params:
        config_file=config_file_path
    shell:
        """
        python {input.script} --config_file {params.config_file}
        """

rule trim_adapters_r1:
    input:
        read=get_input
    output:
        trim_reads1 = "reads/{sample}_trimmed1.fq.gz"
    params:
        three_prime_adapter=config['three_prime_adapter'],
        five_prime_adapter=config['five_prime_adapter'],
        cutadapt_params=config['cutadapt_params1'].split(","),
        nextseq=config['nextseq']
    log:
        out_log = "logs/cutadapt_r1/{sample}.log",
        err_log = "logs/cutadapt_r1/{sample}.err"
    shell:
        """
        # Check if input file exists
        if [ ! -f {input.read} ]; then
            echo "Error: Input file {input.read} does not exist!" >> {log.out_log}
            exit 1
        fi
        
        # Check if cutadapt is available
        if ! command -v cutadapt &> /dev/null; then
            echo "Error: cutadapt could not be found" >> {log.out_log}
            exit 1
        fi

        nextseq_trim=""
        if [ "{params.nextseq}" = "True" ]; then
            nextseq_trim="--nextseq-trim=15"
        fi

        adapters=""
        if [ "{params.three_prime_adapter}" != "" ] && [ "{params.five_prime_adapter}" != "" ]; then
            adapters="-a {params.three_prime_adapter} -g {params.five_prime_adapter}"
        elif [ "{params.three_prime_adapter}" != "" ]; then
            adapters="-a {params.three_prime_adapter}"
        elif [ "{params.five_prime_adapter}" != "" ]; then
            adapters="-g {params.five_prime_adapter}"
        fi

        cutadapt_params="-n 3 -m 15 -M 50 -j 0 --max-n 1"
        if [ ! -z "{params.cutadapt_params}" ] && [ "{params.cutadapt_params}" != "['']" ]; then
            cutadapt_params="{params.cutadapt_params}"
        fi

        echo "Command: cutadapt $adapters $nextseq_trim $cutadapt_params -o {output.trim_reads1} {input.read}" >> {log.out_log}
        cutadapt $adapters $nextseq_trim $cutadapt_params -o {output.trim_reads1} {input.read} >> {log.out_log} 2>> {log.err_log}
        status=$?
        if [ $status -ne 0 ]; then
            echo "cutadapt failed with status $status" >> {log.out_log}
            cat {log.err_log} >> {log.out_log}
            exit $status
        fi
        """
        
if config["cutadapt_params2"] in ["", None]:
    rule move_trimmed_files:
        input:
            trim_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed1.fq.gz"
        output:
            trim_final="reads/{sample}_trimmed_final.fq.gz"
        shell:
            """
            echo "Skipping cutadapt round 2"
            mv {input.trim_read} {output.trim_final}
            """
else:
    rule trim_adapters_r2:
        input:
            trim_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed1.fq.gz"
        output:
            trim_final="reads/{sample}_trimmed_final.fq.gz"
        params:
            cutadapt_params=config['cutadapt_params2'].split(",")
        log:
            out_log = "logs/cutadapt_r2/{sample}.log",
            err_log = "logs/cutadapt_r2/{sample}.err"
        shell:
            """
            # Check if input file exists
            if [ ! -f {input.trim_read} ]; then
                echo "Error: Input file {input.trim_read} does not exist!" >> {log.out_log}
                exit 1
            fi
            cutadapt_params="{params.cutadapt_params}"
            echo "Command: cutadapt $cutadapt_params -o {output.trim_final} {input.trim_read}" >> {log.out_log}
            cutadapt $cutadapt_params -o {output.trim_final} {input.trim_read} >> {log.out_log} 2>> {log.err_log}
            status=$?
            if [ $status -ne 0 ]; then
                echo "cutadapt failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
            
rule collapse_reads:
    input:
        trim_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed_final.fq.gz"
    output:
        collapsed="reads/{sample}_trimmed_collapsed.fa"
    params:
        umi_p=umi_param
    log:
        err_log = "logs/collapsing_reads/{sample}.err"
    shell:
        """
        # Check if input file exists
        if [ ! -f {input.trim_read} ]; then
            echo "Error: Input file {input.trim_read} does not exist!" >> {log.err_log}
            exit 1
        fi
        seqtk seq -a {input.trim_read} | fastx_collapser -i - {params.umi_p} | sed 's/>/>Read/g' > {output.collapsed} 2>> {log.err_log}
        status=$?
        if [ $status -ne 0 ]; then
            echo "Collapsing command failed with status $status" >> {log.err_log}
            exit $status
        fi
        """
        
if config["rem_rep"] in ["True", True, 1, "1"]:
    rule remove_rep:
        input:
            collapsed_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed_collapsed.fa"
        output:
            rm_rep_read="reads/{sample}_rm_rep.fa",
            unaligned_bam="ali/{sample}_reps.bam"
        params:
            star_prefix="ali/{sample}_rm_rep",
            index_path=config["rep_idx"],
            params=config["star_rem_reps_params"].split(","),
            pp2_dir=config["parpipe2_dir"],
            t=num_cores
        log:
            out_log = "logs/rm_reps/{sample}.log",
            err_log = "logs/rm_reps/{sample}.err"
        shell:
            """
            threads="{params.t}"
            prefix="{params.star_prefix}"
            bam="{params.star_prefix}Aligned.out.bam"
            fa="{params.star_prefix}Unmapped.out.mate1"

            if [ ! -f {input.collapsed_read} ]; then
                echo "Error: Input file {input.collapsed_read} does not exist!" >> {log.out_log}
                exit 1
            fi

            rmrep_params="--runMode alignReads --genomeLoad NoSharedMemory --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outSAMtype BAM Unsorted --outFilterType BySJout --outBAMcompression 10 --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --outSAMattributes All --outSAMmode Full"
            if [ ! -z "{params.params}" ] && [ "{params.params}" != "['']" ]; then
                rmrep_params="{params.params}"
            fi
            
            idx_path="{params.pp2_dir}/annotation/rep_idx"
            
            if [ ! -z "{params.index_path}" ] && [ "{params.index_path}" != "['']" ]; then
                rmrep_params="{params.index_path}"
            fi
            
            STAR $rmrep_params --runThreadN $threads --genomeDir $idx_path --outFileNamePrefix $prefix --readFilesIn {input.collapsed_read} >> {log.out_log} 2>> {log.err_log}

            mv $fa {output.rm_rep_read}
            mv $bam {output.unaligned_bam}
            
            status=$?
            if [ $status -ne 0 ]; then
                echo "Removing repetitive elaments failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
else:
    rule copy_collapsed_files:
        input:
            collapsed_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed_collapsed.fa"
        output:
            rm_rep_read="reads/{sample}_rm_rep.fa"
        shell:
            """
            cp {input.collapsed_read} {output.rm_rep_read}
            """

            
if config["aligner"] in ['S', 's', 'star', 'Star', 'STAR']:
    rule genome_alignment_star:
        input:
            rm_rep_read=lambda wildcards: f"reads/{wildcards.sample}_rm_rep.fa"
        output:
            aligned_sam= "ali/{sample}.aligned.mapped.sorted.sam"
            #add all other star output files here (not to rule all)
        params:
            star_prefix="ali/{sample}",
            index_path=get_index_path('star'),
            params=config["star_map_params"].split(","),
            t=num_cores
        log:
            out_log = "logs/star_ali/{sample}.log",
            err_log = "logs/star_ali/{sample}.err"
        shell:
            """
            threads="{params.t}"
            prefix="{params.star_prefix}"
            bam="{params.star_prefix}Aligned.out.bam"
            sam="{params.star_prefix}.aligned.sam"
            
            if [ ! -f {input.rm_rep_read} ]; then
                echo "Error: Input file {input.rm_rep_read} does not exist!" >> {log.out_log}
                exit 1
            fi
            
            star_params="--runMode alignReads --genomeLoad NoSharedMemory --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outSAMattributes All --outSAMtype BAM Unsorted --outFilterType BySJout --outFilterMismatchNoverReadLmax 0.05 --outSAMattrRGline ID:foo --outStd Log --outBAMcompression 10 --outSAMmode Full"
            if [ ! -z "{params.params}" ] && [ "{params.params}" != "['']" ]; then
                star_params="{params.params}"
            fi
            
            STAR $star_params --runThreadN $threads --genomeDir {params.index_path} --outFileNamePrefix $prefix --readFilesIn {input.rm_rep_read} >> {log.out_log} 2>> {log.err_log}
            
            samtools view -F 4 -h -o $sam $bam >> {log.out_log} 2>> {log.err_log}
            samtools sort -o {output.aligned.sam} $sam >> {log.out_log} 2>> {log.err_log}
            
            status=$?
            if [ $status -ne 0 ]; then
                echo "Mapping reads with STAR failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
else:
    rule genome_alignment_bowtie:
        input:
            rm_rep_read=lambda wildcards: f"reads/{wildcards.sample}_rm_rep.fa"
        output:
            aligned_sam= "ali/{sample}.aligned.mapped.sorted.sam"
        params:
            index_path=get_index_path('bowtie'),
            params=config["bowtie_map_params"].split(","),
            t=num_cores
        log:
            out_log = "logs/bowtie_ali/{sample}.log",
            err_log = "logs/bowtie_ali/{sample}.err"
        shell:
            """
            threads="{params.t}"
            
            if [ ! -f {input.rm_rep_read} ]; then
                echo "Error: Input file {input.rm_rep_read} does not exist!" >> {log.out_log}
                exit 1
            fi
            
            bowtie_params="-v 1 -m 10 --best --strata -S --quiet"   
            if [ ! -z "{params.params}" ] && [ "{params.params}" != "['']" ]; then
                bowtie_params="{params.params}"
            fi
            
            bowtie {params.index_path}/bwt $bowtie_params -p $threads -f {input.rm_rep_read} 2>> {log.err_log} | samtools view -hS -F 4 - | samtools sort -O SAM - -o {output.aligned_sam} >> {log.out_log} 2>> {log.err_log}
            
            status=$?
            if [ $status -ne 0 ]; then
                echo "Mapping reads with bowtie failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
            
rule prepare_paralyzer_ini:
    input:
        sam=lambda wildcards: f"ali/{wildcards.sample}.aligned.mapped.sorted.sam"
    output:
        paralyzer_ini="PARalyzer/{sample}.ini"
    params:
        genome_2bit=get_2bit_path(),
        bandwidth=config["BANDWIDTH"],
        conv=config["CONVERSION"],
        min_read_cnt_group=config["MINIMUM_READ_COUNT_PER_GROUP"],
        min_read_cnt_clust=config["MINIMUM_READ_COUNT_PER_CLUSTER"],
        min_read_cnt_kde=config["MINIMUM_READ_COUNT_FOR_KDE"],
        min_clust_size=config["MINIMUM_CLUSTER_SIZE"],
        min_conv_loc_clust=config["MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER"],
        min_conv_cnt_clust=config["MINIMUM_CONVERSION_COUNT_FOR_CLUSTER"],
        min_read_cnt_clust_incl=config["MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION"],
        min_read_l=config["MINIMUM_READ_LENGTH"],
        max_non_conv_mism=config["MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES"],
        signal_ext=config["signal_extension"],
        add_nt_signal=config["ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL"],
        pcr_dups=config["use_pcr_duplicates"]
    shell:
        """
        in_sam="$(pwd)/{input.sam}"
        if [ "{params.pcr_dups}" = "True" ]; then
            in_sam="$(pwd)/{input.sam}=COLLAPSED"
        fi        
        bandwidth="3"
        if [ ! -z "{params.bandwidth}" ] && [ "{params.bandwidth}" != "['']" ]; then
            bandwidth="{params.bandwidth}"
        fi
        conv="T>C"
        if [ ! -z "{params.conv}" ] && [ "{params.conv}" != "['']" ]; then
            conv="{params.conv}"
        fi
        min_read_cnt_group="5"
        if [ ! -z "{params.min_read_cnt_group}" ] && [ "{params.min_read_cnt_group}" != "['']" ]; then
            min_read_cnt_group="{params.min_read_cnt_group}"
        fi
        min_read_cnt_clust="2"
        if [ ! -z "{params.min_read_cnt_clust}" ] && [ "{params.min_read_cnt_clust}" != "['']" ]; then
            min_read_cnt_clust="{params.min_read_cnt_clust}"
        fi
        min_read_cnt_kde="3"
        if [ ! -z "{params.min_read_cnt_kde}" ] && [ "{params.min_read_cnt_kde}" != "['']" ]; then
            min_read_cnt_kde="{params.min_read_cnt_kde}"
        fi
        min_clust_size="11"
        if [ ! -z "{params.min_clust_size}" ] && [ "{params.min_clust_size}" != "['']" ]; then
            min_clust_size="{params.min_clust_size}"
        fi
        min_conv_loc_clust="2"
        if [ ! -z "{params.min_conv_loc_clust}" ] && [ "{params.min_conv_loc_clust}" != "['']" ]; then
            min_conv_loc_clust="{params.min_conv_loc_clust}"
        fi
        min_conv_cnt_clust="2"
        if [ ! -z "{params.min_conv_cnt_clust}" ] && [ "{params.min_conv_cnt_clust}" != "['']" ]; then
            min_conv_cnt_clust="{params.min_conv_cnt_clust}"
        fi
        min_read_cnt_clust_incl="1"
        if [ ! -z "{params.min_read_cnt_clust_incl}" ] && [ "{params.min_read_cnt_clust_incl}" != "['']" ]; then
            min_read_cnt_clust_incl="{params.min_read_cnt_clust_incl}"
        fi
        min_read_l="15"
        if [ ! -z "{params.min_read_l}" ] && [ "{params.min_read_l}" != "['']" ]; then
            min_read_l="{params.min_read_l}"
        fi
        max_non_conv_mism="1"
        if [ ! -z "{params.max_non_conv_mism}" ] && [ "{params.max_non_conv_mism}" != "['']" ]; then
            max_non_conv_mism="{params.max_non_conv_mism}"
        fi
        signal_ext="EXTEND_BY_READ"
        if [ ! -z "{params.signal_ext}" ] && [ "{params.signal_ext}" != "['']" ]; then
            signal_ext="{params.signal_ext}"
        fi
        add_nt_signal=""
        if [ "{params.signal_ext}" = "ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL" ]; then
            add_nt_signal="=2"
            if [ ! -z "{params.add_nt_signal}" ] && [ "{params.add_nt_signal}" != "['']" ]; then
                signal_ext="={params.add_nt_signal}"
            fi
        fi
        
        cat <<EOL > {output.paralyzer_ini}
#
BANDWIDTH=$bandwidth
CONVERSION=$conv
MINIMUM_READ_COUNT_PER_GROUP=$min_read_cnt_group
MINIMUM_READ_COUNT_PER_CLUSTER=$min_read_cnt_clust
MINIMUM_READ_COUNT_FOR_KDE=$min_read_cnt_kde
MINIMUM_CLUSTER_SIZE=$min_clust_size
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=$min_conv_loc_clust
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=$min_conv_cnt_clust
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=$min_read_cnt_clust_incl
MINIMUM_READ_LENGTH=$min_read_l
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=$max_non_conv_mism

$signal_ext$add_nt_signal

GENOME_2BIT_FILE={params.genome_2bit}
SAM_FILE=$in_sam

OUTPUT_DISTRIBUTIONS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.distribution
OUTPUT_GROUPS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.groups
OUTPUT_CLUSTERS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.clusters
OUTPUT_READS_FILE=$(pwd)/PARalyzer/{wildcards.sample}_PARalyzer_Utilized.sam

EOL
        """
        
rule run_paralyzer: 
    input:        
        ini_file=lambda wildcards: f"PARalyzer/{wildcards.sample}.ini"
    output:
        bam="PARalyzer/{sample}_PARalyzer_Utilized.bam",
        bam_idx="PARalyzer/{sample}_PARalyzer_Utilized.bam.bai",
        distr="PARalyzer/{sample}.distribution",
        groups="PARalyzer/{sample}.groups",
        clusters="PARalyzer/{sample}.clusters"
    params:
        mem=config["paralyzer_memory"]
    log:
        out_log = "logs/PARalyzer/{sample}.log",
        err_log = "logs/PARalyzer/{sample}.err"
    shell:
        """
        mem="4G"
        if [ ! -z "{params.mem}" ] && [ "{params.mem}" != "['']" ]; then
            mem="{params.mem}"
        fi
        
        PARalyzer $mem {input.ini_file} >> {log.out_log} 2>> {log.err_log}
        samtools sort -o {output.bam} PARalyzer/{wildcards.sample}_PARalyzer_Utilized.sam  >> {log.out_log} 2>> {log.err_log}
        samtools index {output.bam} >> {log.out_log} 2>> {log.err_log}
        
        status=$?
            if [ $status -ne 0 ]; then
                echo "Running PARalyzer failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
        """

rule calc_conv_stats:
    input:
        bam=lambda wildcards: f"PARalyzer/{wildcards.sample}_PARalyzer_Utilized.bam",
        groups=lambda wildcards: f"PARalyzer/{wildcards.sample}.groups",
        clusters=lambda wildcards: f"PARalyzer/{wildcards.sample}.clusters", 
    output:
        groups_conv="stats/{sample}_groups_conv_stats.tsv",
        clusters_conv="stats/{sample}_clusters_conv_stats.tsv"
    params:
        script=config["parpipe2_dir"] + "/bin/add_conv_stats.py",
        genome=get_genome_path(),
        t=num_cores  
    log:
        out_log = "logs/conv_stats/{sample}.log",
        err_log = "logs/conv_stats/{sample}.err"
    shell:
        """
        python {params.script} --bamfile {input.bam} --groups {input.groups} --clusters {input.clusters} --genome {params.genome} --threads {params.t} 2>> {log.err_log}
        """


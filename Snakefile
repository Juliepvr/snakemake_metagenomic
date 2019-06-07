################################################################################
# Config file and variable declaration.
################################################################################
configfile: "config.yaml"

SAMPLES, = glob_wildcards("data/{sample}_1.fastq.gz")
PATH_DB=config['path_DB']
PATH_ANNOT=config['path_annot']
PATH_RESF=config['path_resf']
SEQ_LEN=config['seq_len']

#create adaptersfile? fastqc?

################################################################################
# Rule "all" executes the entire pipeline when no argument is given.
################################################################################

rule all:
    input:
        expand("10_annotation/{sample}/{sample}_counts_AMR.tsv", sample=SAMPLES),
        expand("06_checkm_result/{sample}/qa_result_checkm.tsv", sample=SAMPLES),
        "01_fastqc/multiqc_report.html"

################################################################################
# 
# Rules to execute the pipeline:
#
################################################################################

################################################################################        
# 1) Run FastQC
################################################################################
rule fastQC:
    input:
        A="data/{sample}_1.fastq.gz",
        B="data/{sample}_2.fastq.gz"
    output:
        "01_fastqc/{sample}_1_fastqc.zip"
    params:
        dir="01_fastqc"
    threads: 4
    shell:
        "fastqc -t {threads} -o {params.dir} {input}"

################################################################################        
# 2) MultiQC to summarize fastQC
################################################################################
rule multiQC:
    input:
        expand("01_fastqc/{sample}_1_fastqc.zip", sample=SAMPLES)
    output:
        "01_fastqc/multiqc_report.html"
    params:
        dir="01_fastqc/"
    shell:
        "multiqc {params.dir} -o {params.dir}"

################################################################################        
# 3) Trim adapters using trimmomatic.
################################################################################
rule trim:
    input:
        A="data/{sample}_1.fastq.gz",
        B="data/{sample}_2.fastq.gz",
        adp="adapters.fa"
    output:
        "02_trimmed/{sample}/{sample}_1_paired.fastq.gz",
        "02_trimmed/{sample}/{sample}_1_unpaired.fastq.gz",
        "02_trimmed/{sample}/{sample}_2_paired.fastq.gz",
        "02_trimmed/{sample}/{sample}_2_unpaired.fastq.gz"
    threads: 8
    shell:
        '''
        trimmomatic PE -threads {threads} -phred33 {input.A} {input.B} {output} \
        ILLUMINACLIP:{input.adp}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
        '''

################################################################################
# 4) Create contigs using megahit, with specified minimum contig length 1500 bp.
################################################################################
rule assemble:
    input:
        A="02_trimmed/{sample}/{sample}_1_paired.fastq.gz",
        B="02_trimmed/{sample}/{sample}_2_paired.fastq.gz"
    output:
        "03_megahit_assembly/{sample}/final.contigs.fa"
    params:
        dir="03_megahit_assembly/{sample}"
    threads: 16
    shell:
# -fo : force outputfolder; without -f megahit requires a non-existing directory. 
# -f is needed because snakemake creates te outputfolder before executing the command
        "megahit -1 {input.A} -2 {input.B} --min-contig-len 1500 -fo {params.dir} -t {threads}"

################################################################################
# 5) Using BWA (Burrows-Wheeler Aligner),
#   align trimmed reads to contigs created by megahit.
################################################################################
rule align:
    input:
        fa="03_megahit_assembly/{sample}/final.contigs.fa",
        A="02_trimmed/{sample}/{sample}_1_paired.fastq.gz",
        B="02_trimmed/{sample}/{sample}_2_paired.fastq.gz"
    output:
        bam="04_align/{sample}/{sample}_results.bwa.bam",
        flag="04_align/{sample}/{sample}_results.flagstat"
    params: 
        tmp="04_align/{sample}/temp.bwa"
    threads: 8
    shell:
        '''
        bwa index {input.fa}
        bwa mem -t {threads} {input.fa} {input.A} {input.B} | samtools view -T \
        {input.fa} -bS - | samtools sort -T {params.tmp} -o {output.bam} -
        samtools flagstat {output.bam} > {output.flag}
        '''

################################################################################
# 6) Calculate coverage & depth of contigs, the jgi script is part of MetaBAT.
################################################################################
rule coverage:
    input:
        "04_align/{sample}/{sample}_results.bwa.bam"
    output:
        "04_align/{sample}/{sample}_results_depth.txt"
    shell: 
        '''
        jgi_summarize_bam_contig_depths --outputDepth {output} --minContigLength \
        2000 --minContigDepth 2 {input}
        '''

################################################################################
# 7) Using MetaBAT, create bins with contigs having similar depth. 
#   Based on contig length and frequency of seq.
################################################################################
rule binning:
    input:
        asm="03_megahit_assembly/{sample}/final.contigs.fa",
        depth="04_align/{sample}/{sample}_results_depth.txt"
    output:
        "05_metabat_bins/{sample}/{sample}_bin.1.fa"
    params:
        dir="05_metabat_bins/{sample}/{sample}_bin"
    threads: 8
    shell: 
        # -m = min contig length, -o = prefix for output
        "metabat2 -m 2000 -i {input.asm} -a {input.depth} -o {params.dir} -t {threads}"

################################################################################
# 8) Using checkm, check bins for marker genes and give taxonomic rank.
# checkm needs +/- 35 G RAM !
# execute snakemake --use-conda to read conda file for checkm rule
# add --cores <int> to make use of threads with <int> >= threads
# checkm uses python 2.7 , other software python 3 --> need seperate conda file 
# with checkm environment.
################################################################################

################################################################################
# 8a) Download the checkm database and set root for checkm.
################################################################################
rule checkm_root:
    output:
        "00_checkm_data/done_tar.log" 
    threads: 4
    params:
        dir="00_checkm_data/"
    conda: "envs/checkm.yaml"
    shell:
        # root needs to be set for every new checkm environment
        # wget -N and tar --skip-old-files prevent redownloading and unpacking same thing
        '''
        wget -N https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P {params.dir}
        tar --skip-old-files -xzvf {params.dir}*.tar.gz -C {params.dir} --index-file={output}
        echo {params.dir} | checkm data setRoot {params.dir}
        '''

################################################################################
# 8b) Checkm lineage_wf: look for taxonomic markers in bins in order to classify 
#     these bins.
################################################################################
rule checkm_lineage:
    input:
        "00_checkm_data/done_tar.log",
        "05_metabat_bins/{sample}/{sample}_bin.1.fa"
    output:
        "06_checkm_result/{sample}/lineage.ms"
    threads: 16
    params:
        dir="06_checkm_result/{sample}/",
        bins="05_metabat_bins/{sample}/"
    conda: "envs/checkm.yaml"
    shell:
    # add -f to save output in default outputtype 1
        "checkm lineage_wf -t {threads} -x fa {params.bins} {params.dir}"

################################################################################
# 8c) Checkm qa: this step is incorporated by the previous command, but is 
#     repeated to generate an extended summary of bin quality, specified by -o 2.
################################################################################
rule checkm_qa:
    input:
        "06_checkm_result/{sample}/lineage.ms"
    output:
        "06_checkm_result/{sample}/qa_result_checkm.tsv"
    threads: 4
    params:
        dir="06_checkm_result/{sample}/"
    conda: "envs/checkm.yaml"
    shell:
    # -o 2 for tsv version of output
        "checkm qa -t {threads} {input} {params.dir} -f {output} --tab_table -o 2"

################################################################################
# 9) Prodigal to predict proteins from contigs of step 2.
################################################################################
rule predict:
    input:
        "03_megahit_assembly/{sample}/final.contigs.fa"
    output:
        meta="07_prodigal/{sample}/{sample}_metadata_genes.txt",
        pred="07_prodigal/{sample}/{sample}_predicted_genes.fa"
    #threads: 4 # no option for threads
    shell:
        "prodigal -i {input} -o {output.meta} -a {output.pred}"

################################################################################
# 10) Run diamond results against database
################################################################################
rule diamond_kegg:
    input:
        db= PATH_DB,
        pred="07_prodigal/{sample}/{sample}_predicted_genes.fa"
    output:
        "08_diamond_blastP/{sample}/{sample}_blastp_matches_kegg.m8"
    threads: 8
    shell:
        "diamond blastp --max-target-seqs 10 -d {input.db} -q {input.pred} -o {output} --threads {threads}"

################################################################################
# 11) Run diamond results against original reads
################################################################################
########################
# 11a) make protein DB #
########################
rule makedb:
    input:
        "07_prodigal/{sample}/{sample}_predicted_genes.fa"
    output: 
        "09_diamond_blastX/{sample}/{sample}_predicted.dmnd"
    threads: 16
    shell:
        "diamond makedb --in {input} -d {output} --threads {threads}"

####################
# 11b) run blastX  #
####################
rule blastx:
    input:
        db="09_diamond_blastX/{sample}/{sample}_predicted.dmnd",
        reads="02_trimmed/{sample}/{sample}_1_paired.fastq.gz"
    output:
        "09_diamond_blastX/{sample}/{sample}_blastx_matches_reads.m8"
    threads: 16
    shell:
        "diamond blastx --max-target-seqs 1 -d {input.db} -q {input.reads} -o {output} --threads {threads}"

################################################################################
# 12) Annotate blastP results from step 10
################################################################################
rule annotate_blastP:
    input:
        m8="08_diamond_blastP/{sample}/{sample}_blastp_matches_kegg.m8"
    output:
        "10_annotation/{sample}/{sample}_summary.tsv"
    params:
        dir="10_annotation/{sample}/{sample}_",
        annot=PATH_ANNOT
    shell:
        "scripts/kegg_info.py {input.m8} {params.annot} {params.dir}"

################################################################################
# 13) Count and annotate blastX results from step 11
################################################################################
rule count_blastX:
    input:
        m8="09_diamond_blastX/{sample}/{sample}_blastx_matches_reads.m8",
        sl=SEQ_LEN,
        need="10_annotation/{sample}/{sample}_summary.tsv"
    output:
        "10_annotation/{sample}/{sample}_counts_protein"
    params:
        dir="10_annotation/{sample}/{sample}_"
    shell:
        "scripts/count_annotate.py {input.m8} {params.dir} {input.sl} "

################################################################################
# 14) Look for AMR genes using ABRicate
################################################################################
rule abricate:
    input:
        fa="03_megahit_assembly/{sample}/final.contigs.fa"
    output:
        tsv="11_AMR/{sample}/{sample}_AMR.tsv"
    shell:
        "abricate {input.fa} > {output.tsv}"

################################################################################
# 15) Annotate AMR genes
################################################################################
rule amr:
    input:
        cts="10_annotation/{sample}/{sample}_counts_protein",
        amr="11_AMR/{sample}/{sample}_AMR.tsv",
        resf=PATH_RESF
    output:
        tsv="10_annotation/{sample}/{sample}_counts_AMR.tsv"
    shell:
        "scripts/count_amr.py {input.cts} {input.amr} {input.resf} {output.tsv}"


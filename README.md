# Reproducible Snakemake pipelines for metagenomic analysis in farm animals

Metagenomic analysis is generally time consuming and requires the implementation of multiple tools. A pipeline stitching the different steps together and allowing the jobs to run in parallel would be time-saving and guarantee the exact same process for every sample.

Using a subset of an Illumina HiSeq 3000 dataset, generated from 359 European pig herds’ and poultry flocks’ gut microbiomes, a bioinformatics pipeline for metagenomic analysis was developed. Using conda environments and Snakemake, all steps from quality control to genome binning and annotation were incorporated. Because of their significance in the food chain, detection of antimicrobial resistance genes was also included. By inspecting pathways, enzymes and CAZymes, a clearer understanding of the gut microbiome is intended. 

User friendliness was guaranteed by using conda. Apart from anaconda, no installations are required.
The pipeline was designed for use with metagenomic, paired-end, short read sequences.

An internet connection and at least 35Gb of RAM memory are necessary.

## STEP 1 : Setup directory to run snakemake. 

This directory should contain: 
- **"data"** : either a link called "data", to the fastq files ; or a directory called "data" , containing the fastq files; format {sample}_1.fastq.gz and {sample}_2.fastq.gz
- **"envs"** : a directory containing the conda environments env.yaml and checkm.yaml
- **"scripts"** : a directory containing the python scripts to be used by snakemake
- *"config.yaml"* : a file containing paths to a database.dmnd (created with diamond makedb, see https://github.com/bbuchfink/diamond), the directory containing annotation files and the sequence length, adjust this to your own data
- *"Snakefile"* : the snakefile itself, containing the tasks and their execution order
- *"cluster.json"* : a file containing the specifications per task, such as time allowed to run, nr of cores, project, ...
- *"custombash.sh"* : a bash script to setup the environment on a cluster
- *"runjob.sh"* : the script containing the command to send the snakemake jobs to a cluster via qsub, check if the commands match these for the server you are working on
- *"adapters.fa"* : a fasta file containing the adapters to be trimmed

## STEP 2: create the base conda environment

requires anaconda (https://www.anaconda.com/)

`$ conda env create -f envs/env.yaml `

## STEP 3 : Start a terminal multiplexer 

e.g. `$ tmux new -s snakemake `

## STEP 4 : Start the job 

make sure you are in the right directory and 

- run the pipeline as a whole:

	`$ conda activate meta-assembly `

	`$ snakemake --use-conda --cores <int> `

- or  submit jobs to the server (advisable with lots of samples):

	`$ ./runjob.sh `

## STEP 5 : Sit back and relax! 
Maybe book a cruise to the carribean, you deserve it

      _                   .-=-.          .-==-.
     { }      __        .' O o '.       /  -<' )
     { }    .' O'.     / o .-. O \     /  .--v`
     { }   / .-. o\   /O  /   \  o\   /O /
      \ `-` /   \ O`-'o  /     \  O`-`o /
       `-.-`     '.____.'       `.____.'

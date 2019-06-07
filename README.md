# STEP 1 : Setup directory to run snakemake. 

This directory should contain: 
- **"data"** : either a link called "data", to the fastq files ; or a directory called "data" , containing the fastq files
- **"envs"** : a directory containing the conda environments env.yaml and checkm.yaml
- **"scripts"** : a directoy containing the python scripts to be used by snakemake
- *"config.yaml"* : a file containing paths to a db, annotation files and sequence length, adjust this to your own data
- *"Snakefile"* : the snakefile itself, a seperate one in this case with specified threads per task
- *"cluster.json"* : a file containing the specifications per task, such as time allowed to run, nr of cores, project, ...
- *"custombash.sh"* : a bash script to setup the environment on the cluster
- *"runjob.sh"* : the script containing the command to send the snakemake jobs to a cluster via qsub, check if the commands match these for the server you are working on

# STEP 2: create the base conda environment

requires anaconda (https://www.anaconda.com/)

`$ conda env create -f envs/env.yaml `

# STEP 3 : Start a terminal multiplexer 

e.g. `$ tmux new -s snakemake `

# STEP 4 : Start the job 

make sure you are in the right directory and 

run the pipeline as a whole:

`$ conda activate meta-assembly `

`$ snakemake --use-conda --cores <int> `

or  submit jobs to the server (advisable with lots of samples):

`$ ./runjob.sh `

# STEP 5 : Sit back and relax! 
Maybe book a cruise to the carribean, you deserve it

      _                   .-=-.          .-==-.
     { }      __        .' O o '.       /  -<' )
     { }    .' O'.     / o .-. O \     /  .--v`
     { }   / .-. o\   /O  /   \  o\   /O /
      \ `-` /   \ O`-'o  /     \  O`-`o /
       `-.-`     '.____.'       `.____.'

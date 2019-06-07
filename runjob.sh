#!/bin/bash
source $HOME/.bashrc

# load sge and anaconda
module load sge
module load anaconda/5.0.1

# activate the conda environment
source activate meta-assembly

snakemake --keep-going --rerun-incomplete --jobscript custombash.sh \
--use-conda --cluster-config cluster.json --cluster "qsub -R yes \
-V -S /bin/bash -cwd -pe sharedmem {cluster.core} -l h_rt={cluster.time} \
-l h_vmem={cluster.vmem} -P {cluster.proj}" --jobs 5000 -s Snakefile


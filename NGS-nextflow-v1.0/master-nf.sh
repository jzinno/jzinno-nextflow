#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=nextflow
#SBATCH --mail-type=FAIL,END
#SBATCH --output=nextflow-log_%j.out


module purge

module load singularity-ce/3.9.2

module load nextflow/21.10.6

nextflow main.nf -with-report report-nextflow-log.html -with-dag flowchart.html -with-timeline timeline.html -resume

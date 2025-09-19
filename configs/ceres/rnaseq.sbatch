#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p long
#SBATCH -A musca
#SBATCH -o "rnaseq.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "rnaseq.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=edwardbird@ksu.edu
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load apptainer
module load nextflow


nextflow run -resume /project/culicoides/bird_projects/pipelines/rnaseq -c /project/culicoides/bird_projects/pipelines/rnaseq/conf/ceres.cfg --input ./trim_samplesheet.csv --outdir exapx_retrim_data --remove_ribo_rna true -N edwardbird@ksu.edu --gtf ./cson_egapx.gtf --fasta ./cson_egapx.fna -with-report nextflow.report.html -with-timeline nextflow.timeline.html
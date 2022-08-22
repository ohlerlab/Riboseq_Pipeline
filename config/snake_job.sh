#! /bin/bash
#$ -cwd #start from current directory
#$ -l h_vmem=20G -l h_rt=48:00:00 -l os=centos7 -pe smp 2
#$ -V #export all the environmental variables into the context of the job
#$ -j yes #merge the stderr with the stdout
#$ -o logs/ #stdout, job log
#$ -m ea # send email beginning, end, and suspension
#$ -M gabriel.villamil@mdc-berlin.de
#$ -N 'myjob'


source ~/.bashrc
mkdir -p sge_log
conda activate snakemake
snakemake -j 10 -k -p --restart-times 1 --max-jobs-per-second 5 -s Snakefile \
	 --cluster-config ../config/config_pipeline.json  --rerun-incomplete \
	 --use-singularity --singularity-args "-B /fast/AG_Ohler/:/fast/AG_Ohler"  \
	 --cluster="qsub -cwd -V -l m_mem_free={cluster.m_mem_free} -l h_rt={cluster.h_rt} -pe {cluster.pe} -j yes -o sge_log" "all"

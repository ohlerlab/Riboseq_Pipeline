# Riboseq_Pipeline
This is the Lab's standard Ribo-seq processing pipeline, a snakemake workflow and some R scripts.


Running the pipeline consist of _ main steps

1. Set up
	- Clone this repo

	- (optionally) copy in the input fastq files to this directory, if they don't already have a stable home elsewhere.
	- modify 'read_files.csv'  to point at the files you're using.
		- it has four columns but the last two - 'mate' and 'pair_id' can be left blank if you're using single end data
		- the first two are mandatory - sample_id groups files that will be grouped together for QC, alignment etc.
		- 'sample_id' links the files to the lines of 'sample paramter'
    - sample parameter should be modified to include any metadata and covariates you want to model, as well as fragment length (or just use defaults) for RNA samples
	- sample_parameter.csv also contains data on whether things are strand specific or not, and if they are paired.
	- Now modify config.yaml to point at the annotation etc. you want to use, as well as the locations of RiboseQC and ORFquant (they can be cloned if they are missing)
	- Sink the snakemake file to the pipeline folder with ln -s ../src/pipeline.smk Snakefile. Now you can do snakemake -n instead of snakemake -s ../src/pipeline.smk -n.
	- Create a conda environment from the conda ../ribopipebase.yml file using `conda env create -f ../ribopipebase.yml`, this may take a while because conda. and activate the envoriment before running the pipeline manually `conda activate ribopipe`. This will be done automaticall when running the snake_job.
	- Install ORFquant and RiboseQC
		- when not in the ribobase conda enviroment `conda activate ribopipe`
		- `R -e 'library("devtools");install_github(repo = "ohlerlab/ORFquant")'`
		- `conda install -c bioconda riboseqc`

2. Running the pipeline
    - first do a 'dry run' and see if the snakefile works, without runnning commands - `snakemake -n`. Debug
    - The pipeline is run by doing ```bash ../src/snake_job.sh```
3. Inevitable debugging
    - The bugs you will face (and you will face them) can be divided into several main catagories.
    1. Bugs at the level of Snakemake
        - often the snakemake won't run the way you think it should, this is often due to irregularities in the input files (did you use a commma where a space should be, typos in sample ids, etc etc) test run the snakemake by doing snakemake -n
    2. Bugs in the bash code for a rule - for these it's often useful to do snakemake -p -j2 problem_file - this will show you the code being run, you can copy paste it to run line by line and find out what's wrong.
    

#copy in files (so nodes can see them)
	find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | xargs -I{} rsync -avshP {} input/
#make a read_files.csv - note that the sed call there does some redex surgery on the file names, modify it, and the rest of the code, as needed to deal with your files and directory structure
    (head -n1 src/read_files.csv.example ; find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | while read f ; do  f=$(pwd)/input/$(basename $f);echo $(basename $f | sed 's/_S[0-9]_R[0-9]_.*fastq.gz//' )","$f ; done ) | tee /dev/stderr  > src/read_files.csv



#code to check cutaddapt workedd
    ```grep Summary -A7 pipeline/cutadapt_reads/*/*fastq.gz.cutadaptstats.txt```
#code to see how many reads lost to collaps_reads
```Sys.glob('pipeline/collapse_reads/*/*.fastq.gz.collreadstats.txt')%>%setNames(.,basename(dirname(.)))%>%map(readLines)%>%map(head,4)%>%map(tail,2)%>%map(str_extract,'\\d+')%>%simplify2array%>%t%>%set_colnames(c('input','uniq'))%>%as.data.frame(stringsAsFactors=F)%>%rownames_to_column('sample')%>%mutate(unique = round(as.numeric(uniq)/as.numeric(input),3))```


HOLY SHIT- I Have RIboseq data! What do I do with it???

1) get yourself a node with say 10 cores and 10 megs each
qrsh -V -now no -pe smp 8 -l  m_mem_free=10G

2) Clone the Riboseq pipeline into a new folder.
cd /fast/AG_Ohler/dharnet
git clone https://github.com/ohlerlab/Riboseq_Pipeline/ SRP_Riboseq
cd SRP_Riboseq

2) install the base conda environment if it's not already installed, and then activate it
conda create -f ribopipebase.yaml -n ribopipe; conda activate ribopipe

3)create a branch for this project
git branch SRP_Riboseq
git checkout SRP_Riboseq

4) Edit src/read_files.csv - you can do this by hand, or just use a bash loop like..

echo "sample_id,file_id,mate,pair_id" > src/read_files.csv; for fastq in /fast/AG_Landthaler/emanuel/Jonas_ER_Riboseq/*fastq ; do echo $(basename ${fastq%.fastq}),$fastq ; done >> src/read_files.csv

It has columns sample_id,file_id,mate,pair_id - the first groups technical replicates, the second is the file path, the 3rd identifys which end for paired end reads, and the last one groups paired end samples.

5) Edit sample_parameter.csv. This has a few columns

sample_id - the id which links fastq files to a biological sample. Paired-end mates, as well as technical replicates (resequencing) will have the same sample_id and e.g. fastqc will combine them.

libtype - see here, this tells salmon (and some other rules in the snakemake file) what kind of library it is
this will vary from library to library, but most frequently you'll have SF for Riboseq (single end forward strand) and SR or SF for the matching RNAseq.
https://salmon.readthedocs.io/en/latest/library_type.html

group - this column groups your biological replicates together.
isriboseq - this column should be either True or False

6) Edit config.yaml - see the comments in file, this is wehre you put in e.g. the path to annotation, genome sequence, etc.

7) make a pipeline directory - mkdir pipeline;cd pipeline
8) make a link to the snakefile - ln -s ../src/pipeline.smk Snakefile
9) do snakemake -n  - this will run the snakemake code without executing the rules. Most likely there'll be misspecified paths etc the first time, you may need to look at src/pipeline.smk to figure out what's wrong, or even drop a set_trace() somewhere inside it and then e.g. look at sampledf or seqfilesdf to make sure they look right.


....


10) When snakemake -n runs set the whole thing to run on the cluster with bash ../src/snake_job.sh

#when this fails you can see the errors for specific jobs by looking at the log files in pipeline/sge_logs/
#it's often useful to just rerun a specific file on your node, without submitting ot hte cluster, by doing snakemake -p -j2 problem_file - this will show you the command that's being fun and the error message.
#Cluster logs etc.
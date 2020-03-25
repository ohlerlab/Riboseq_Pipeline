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
    

#copy in files (so nodes can see them)
	find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | xargs -I{} rsync -avshP {} input/
#make a read_files.csv - note that the sed call there does some redex surgery on the file names, modify it, and the rest of the code, as needed to deal with your files and directory structure
    (head -n1 src/read_files.csv.example ; find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | while read f ; do  f=$(pwd)/input/$(basename $f);echo $(basename $f | sed 's/_S[0-9]_R[0-9]_.*fastq.gz//' )","$f ; done ) | tee /dev/stderr  > src/read_files.csv



#code to check cutaddapt workedd
    ```grep Summary -A7 pipeline/cutadapt_reads/*/*fastq.gz.cutadaptstats.txt```
#code to see how many reads lost to collaps_reads
```Sys.glob('pipeline/collapse_reads/*/*.fastq.gz.collreadstats.txt')%>%setNames(.,basename(dirname(.)))%>%map(readLines)%>%map(head,4)%>%map(tail,2)%>%map(str_extract,'\\d+')%>%simplify2array%>%t%>%set_colnames(c('input','uniq'))%>%as.data.frame(stringsAsFactors=F)%>%rownames_to_column('sample')%>%mutate(unique = round(as.numeric(uniq)/as.numeric(input),3))```

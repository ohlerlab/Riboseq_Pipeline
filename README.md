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
	- 

#copy in files (so nodes can see them)
	find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | xargs -I{} rsync -avshP {} input/
#make a read_files.csv - note that the sed call there does some redex surgery on the file names, modify it, and the rest of the code, as needed to deal with your files and directory structure
    (head -n1 src/read_files.csv.example ; find /data/ohler/seqData/819-836/RiboFP_rRNA_depletion_test/ -iname '*fastq.gz' | grep R1  | while read f ; do  f=$(pwd)/input/$(basename $f);echo $(basename $f | sed 's/_S[0-9]_R[0-9]_.*fastq.gz//' )","$f ; done ) | tee /dev/stderr  > src/read_files.csv

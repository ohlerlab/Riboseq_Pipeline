# Riboseq_Pipeline
This is the Lab's standard Ribo-seq processing pipeline, a snakemake workflow and some R scripts, including ORFquant and RiboseQC.

- On the max cluster, get yourself a node with say 10 cores and 10 gigs each
``qrsh -V -now no -pe smp 8 -l  m_mem_free=10G``

- Clone the Riboseq pipeline into a new folder.
``cd /fast/AG_Ohler/{$USER}``
``git clone https://github.com/ohlerlab/Riboseq_Pipeline/ Pipeline_Demo``
`cd Pipeline_Demo`

- install the base conda environment if it's not already installed, and then activate it
``conda create -f ribopipebase.yaml -n ribopipe;
conda activate ribopipe``

- create a branch for this project
`git branch Pipeline_Demo 
git checkout Pipeline_Demo`

- install RiboseQC and ORFquant
``mkdir Applications``
``git clone https://github.com/ohlerlab/RiboseQC.git Applications/RiboseQC``
``git clone https://github.com/ohlerlab/ORFquant.git Applications/ORFquant``

- Edit ``src/read_files.csv`` to point it at your fastq files (zipped or not, either is fine)- you can do this by hand, or just use a bash loop like..

``echo "sample_id,file_id,mate,pair_id" > src/read_files.csv; for fastq in /fast/dharnet/Ebp1_Riboseq/input/*fastq.gz ; do echo $(basename ${fastq%.fastq.gz}),$fastq ; done >> src/read_files.csv``

It has columns: ``sample_id,file_id,mate,pair_id`` - the first groups technical replicates, the second is the file path, the 3rd identifys which end for paired end reads (1 or 2) , and the last one groups paired end samples.

- Edit sample_parameter.csv. This has a few columns

sample_id - the id which links fastq files to a biological sample. Paired-end mates, as well as technical replicates (resequencing) will have the same sample_id and e.g. fastqc will combine them.

libtype - see here, this tells salmon (and some other rules in the snakemake file) what kind of library it is
this will vary from library to library, but most frequently you'll have SF for Riboseq (single end forward strand) and SR or SF for the matching RNAseq.
https://salmon.readthedocs.io/en/latest/library_type.html

group - this column groups your biological replicates together.
isriboseq - this column should be either True or False

- Edit config.yaml - see the comments in file, this is wehre you put in e.g. the path to annotation, genome sequence, etc.

- make and enter a pipeline directory - ``mkdir pipeline;cd pipeline``
- make a link to the snakefile - ``ln -s ../src/pipeline.smk Snakefile``
9) do ``snakemake -n``  - this will run the snakemake code without executing the rules.  somewhere inside it and then e.g. look at sampledf or seqfilesdf to make sure they look right.
- Running the pipeline
    - first do a 'dry run' and see if the snakefile works, without runnning commands - `snakemake -n`.
    - Most likely there'll be misspecified paths etc the first time, you may need to look at src/pipeline.smk to figure out what's wrong, or even drop a set_trace()
    - When that works, The pipeline is run by doing ```bash ../src/snake_job.sh```
- Inevitable debugging
    - The bugs you will face (and you will face them) can be divided into several main catagories.
    - Bugs at the level of Snakemake
        - often the snakemake won't run the way you think it should, this is often due to irregularities in the input files (did you use a commma where a space should be, typos in sample ids, etc etc) test run the snakemake by doing snakemake -n
    - Bugs in the bash code for a rule - for these it's often useful to do snakemake -p -j2 problem_file - this will show you the code being run, you can copy paste it to run line by line and find out what's wrong.
    

.... a few iterations later


- When snakemake -n runs set the whole thing to run on the cluster with ``bash ../src/snake_job.sh`` (I generally run interactively with screen, so I get colored feedback - you can also qsub this script)

	- when this fails you can see the errors for specific jobs by looking at the log files in ``sge_logs/``
	- it's often useful to just rerun a specific file on your node, without submitting ot hte cluster, by doing ``snakemake - p -j2 problem_file`` - this will show you the command that's being fun and the error message.

#some useful code snippets

#code to check cutaddapt workedd
    ```grep Summary -A7 pipeline/cutadapt_reads/*/*fastq.gz.cutadaptstats.txt```
#code to see how many reads lost to collaps_reads
```Sys.glob('pipeline/collapse_reads/*/*.fastq.gz.collreadstats.txt')%>%setNames(.,basename(dirname(.)))%>%map(readLines)%>%map(head,4)%>%map(tail,2)%>%map(str_extract,'\\d+')%>%simplify2array%>%t%>%set_colnames(c('input','uniq'))%>%as.data.frame(stringsAsFactors=F)%>%rownames_to_column('sample')%>%mutate(unique = round(as.numeric(uniq)/as.numeric(input),3))```


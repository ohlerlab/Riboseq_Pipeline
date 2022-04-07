# Ribo-Seq Pipeline

This is the lab's standard Ribo-Seq processing pipeline. It consists of a [docker](7)ised [snakemake](6) workflow that ties together R scripts, as well as [ORFquant](3), [Ribostan](4) and [RiboseQC](5) and can be run in [Singularity](8).  

## Installation

1. **Install Singularity:**  

    [Link to guide.](1)   
    *Note: If you have conda installed, it might interfer with Singularity.  
    Add this to the end of your path, e.g. `.bashrc` on linux, to avoid this behavior:*
    ```
    export PATH=$PATH:/usr/bin
    export PATH=$PATH:/usr/sbin
    if [ $SINGULARITY_NAME ]; then
        echo "in a singularity container, miniconda and guix entries are removed from path";
        export PATH=$(echo $PATH | tr ':' '\n' | grep  -v '/guix' | grep -v '/miniconda' | tr '\n' ':' )
        echo $PATH
    fi
    ```

2. **Install our lab's RiboseQC, ORFquant, and Ribostan packages:**  

    By default, the pipeline will look for a folder above the project folder called Applications, so create this folder (in e.g. `/fast/AG_Ohler/user/Applications`) and populate it like so:  
    ```
    mkdir Applications #create a folder in the current directory
    git clone https://github.com/ohlerlab/RiboseQC.git Applications/RiboseQC
    git clone https://github.com/ohlerlab/ORFquant.git Applications/ORFquant
    git clone https://github.com/zslastman/Ribostan.git Applications/Ribostan
    ```

3. **Accessing the container for running locally:**  
    ```
    singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://dermotharnett/riboseq_pipeline
    ```
    *Note: the `-B` entry mounts file paths, so that Docker can see the folders.*   

    It will take a while to download all the necessary programs and libraries the first time.


## Usage

### Initial configuration

1. Edit [sample_config.tsv](/README_new.md#sampleconfigtsv) in the folder `/src/` to point it to your fastq files (zipped or unzipped) with the appropriate parameters.

2. Edit [config.yaml.](/README_new.md#configyaml) This is where you put in for example, the path to annotation, genome sequence, etc.

3. Dry run of the pipeline  
    Make and enter a pipeline directory: 
    ```
    mkdir pipeline;cd pipeline
    ```
    Make a link to the snakefile:  
    ```
    ln -s ../src/pipeline.smk Snakefile
    ```
    Run the snakemake code without executing the rules:
    ```
    snakemake -n
    ```  
    Check whether your output files (sampledf or seqfilesdf) look right.  
    A common reason for errors are misspecified paths or misplaced symbols like **'. ,'**.   

### Debugging

- You can look at the individual rules (code run for a specific file) in the snakemake file `/src/pipeline.smk`. You can also rerun a specific file without submitting it to the cluster. Using this approach will show you the command that's being run and the error message:
    ```
    snakemake -p -j2 problem_file
    ```
    *Note: Errors for specific jobs are saved in the log files in `/sge_logs/`*

- Command to check *cutadapt* worked:
    ```
    grep Summary -A7 pipeline/cutadapt_reads/*/*fastq.gz.cutadaptstats.txt
    ```

- Command to see how many reads were lost to *collaps_reads*:
    ```
    Sys.glob('pipeline/collapse_reads/*/*.fastq.gz.collreadstats.txt')%>%setNames(.,basename(dirname(.)))%>%map(readLines)%>%map(head,4)%>%map(tail,2)%>%map(str_extract,'\\d+')%>%simplify2array%>%t%>%set_colnames(c('input','uniq'))%>%as.data.frame(stringsAsFactors=F)%>%rownames_to_column('sample')%>%mutate(unique = round(as.numeric(uniq)/as.numeric(input),3))
    ```

### Running the pipeline
- After the initial configuration and dry run use:
    ```
    bash ../src/snake_job.sh
    ```

    *Note: This command can be run interactively with screen for color-coded feedback or with qsub.*

## Development

>**What are Docker and Singularity?**  
[Docker and Singularity](9) are both container management tools. A container is a layer you run programs with that allows you to virtualize everything below the operating system - i.e. you can run programs as if you have installed a totally different set of software. Docker needs root access for installation and usage which is a problem for server-sided use where permissions are often restricted. Singularity is able to use Docker containers but without root access.  

1. Build a container on a system where you have root access.
2. Make this container available to the cluster (e.g. by uploading it to Docker hub). 
3. Use singularity on the cluster to download and use this container, thus obviating the need to install software locally.  

### Making changes to the container  

1. Install Docker
2. Create or login to a profile on Docker hub.
3. Create a folder and put the ‘DOCKERFILE’ from the repo inside.
4. Edit it as needed.  
    *Note: When editing, add lines AFTER the existing ones or it will have to rerun everything above (it will run everything on the first time regardless).*  
    E.g., to install gplot2, add 
    ```
    RUN R -e ‘BiocManager::install(c(“ggplot2”))’ 
    ```
5. Build the container with: 
    ```
    docker build -t YOUR_DOCKERHUB_NAME/riboseq_pipeline
    ```
6. Push it to Docker hub with: 
    ```
    docker push YOUR_DOCKERHUB_NAME/riboseq_pipeline
    ```
7. To make snakemake refresh the container, delete `myproject/pipeline/.snakemake/singularity`.
8. Either   
    - step into the container to run snakemake from there:
        ```
        singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://YOUR_DOCKERHUB_NAME/riboseq_pipeline
        ```  
    - or use the snake_job script:
        ```
        bash ../src/snake_job.sh
        ```
        >*Note: The script passes the flags `—use-singularity  --singularity-args "-B /fast/AG_Ohler/:/fast/AG_Ohler/"` to snakemake so that each cluster node uses it.*

### Tips for docker environments

- It’s easy to unintentionally push the incorrect Docker environment as you build.  
Make sure the commands are correct, especially the tag names. You can go on Docker hub and go to tags `>` latest, and see what commands were run to make the container. Click on each line to see the complete command.
- Often, problems can occur when other entries in your `.bashrc` change your path.


## Features
TODO

## Configuration

### `sample_config.tsv`

|COLUMNS |Description|
|:---:|:---:|
|``sample_id``|The id which links fastq files to a biological sample. Paired-end mates,as well as technical replicates (resequencing) will have the same sample_id and fastqc will combine them.|
|``file_id``| The file path|
|``mate``|Identifier for paired end reads (1 or 2)|
|``pair_id``|Groups paired end samples (the pair get the idential number)|
|`libtype`|Tell salmon (and some other rules in the snakemake file) what kind of library it is, [see here](2). This will vary from library to library, but most frequently you'll have SF for Riboseq (single end forward strand) and SR or SF for the matching RNAseq.|
|`group`|Groups the biological replicates together.|
|`isriboseq`|Should be either True or False when the sample is Riboseq or RNAseq, respectively.|  

### `config.yaml`
TODO add some options for tuning the pipeline like  
- `FILT_RIBO_FOLDER`  
- `FILT_RNA_FOLDER`

## Contributing
TODO

## Links
https://sylabs.io/guides/3.0/user-guide/installation.html 'installation guide for singularity'  
https://salmon.readthedocs.io/en/latest/library_type.html 'salmon libtype documentation'  
https://github.com/ohlerlab/ORFquant.git 'ORFquant Git'  
https://github.com/zslastman/Ribostan.git 'Ribostan Git'  
https://github.com/ohlerlab/RiboseQC.git 'RiboseQC Git'  
https://snakemake.readthedocs.io/en/stable/ 'Snakemake docs'  
https://docs.docker.com/get-started/overview/ 'What is Docker'  
https://sylabs.io/guides/2.6/user-guide/introduction.html 'Singularity intro'  
https://sylabs.io/guides/2.6/user-guide/singularity_and_docker.html 'Singularity with Docker'

[1]: https://sylabs.io/guides/3.0/user-guide/installation.html 'installation guide for singularity'  
[2]: https://salmon.readthedocs.io/en/latest/library_type.html 'salmon libtype documentation'  
[3]: https://github.com/ohlerlab/ORFquant.git 'ORFquant Git'  
[4]: https://github.com/zslastman/Ribostan.git 'Ribostan Git'  
[5]: https://github.com/ohlerlab/RiboseQC.git 'RiboseQC Git'  
[6]: https://snakemake.readthedocs.io/en/stable/ 'Snakemake docs'  
[7]: https://docs.docker.com/get-started/overview/ 'What is Docker'  
[8]: https://sylabs.io/guides/2.6/user-guide/introduction.html 'Singularity intro'  
[9]: https://sylabs.io/guides/2.6/user-guide/singularity_and_docker.html 'Singularity with Docker'

## Licensing
TODO
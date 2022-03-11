# Riboseq_Pipeline

This is the lab's standard Ribo-seq processing pipeline. It consists of a dockerised snakemake workflow that ties together R scripts, as well as ORFquant, Ribostan and RiboseQC.


## Installation

1. Install Singularity: [guide][1] 
    - *If you have conda installed, add this to the end of your path (.bashrc) for conda to not interfere with singularity*
    ```
    export PATH=$PATH:/usr/bin
    export PATH=$PATH:/usr/sbin
    if [ $SINGULARITY_NAME ]; then
        echo "in a sing container,removing miniconda and guix entries from path";
        export PATH=$(echo $PATH | tr ':' '\n' | grep  -v '/guix' | grep -v '/miniconda' | tr '\n' ':' )
        echo $PATH
    fi
    ```

2. Install our lab's packages RiboseQC, ORFquant, and Ribostan. 
    - By default, the pipeline will look for a folder above the project folder called Applications, so create this folder (in e.g. `/fast/AG_Ohler/user/Applications` ) and populate it like so:  
    ```
    mkdir Applications # this creates a folder in the current directory
    git clone https://github.com/ohlerlab/RiboseQC.git Applications/RiboseQC
    git clone https://github.com/ohlerlab/ORFquant.git Applications/ORFquant
    git clone https://github.com/zslastman/Ribostan.git Applications/Ribostan
    ```

3. To run locally on one node you can just step inside the container:
    `singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://dermotharnett/riboseq_pipeline`
    - Note: the ‘-B’ entry mounts file paths, so Docker can see the folders.
    - It will take a while to download all the necessary programs and libraries the first time.


## Usage

### Initial configuration

this section is under construction


### Docker and singularity - a brief overview

Docker and singularity are both container management tools. A container is a layer you run programs with that allows you to virtualize everything below the operating system - i.e. you can run programs as if you have installed a totally different set of software. Docker needs root access for installation and usage which is a problem for server-sided use where permissions are often restricted. Singularity is able to use Docker containers but without root access.

#### The principle for working with this pipeline

1. Build a container on some system where you have root access - your MacBook, a workstation.
2. Make this container available to the cluster - e.g. by uploading it to Docker hub. 
3. Use singularity on the cluster to download and use this container - thus obviating the need to install software locally and prevent crashes.


### Making changes to the container

This is a solution for those using a MAC - on a linux machine this might be easier.
1. Install Docker
2. Create a profile on Docker hub, or use a shared one if you want.
3. Create a folder, and put the ‘DOCKERFILE’ from the repo in there.
4. Edit it as needed. Note that when editing you should add lines AFTER the existing ones or it will have to rerun everything above, (which it will have to do anyway the first time). Format should be self explanatory. It first loads an existing one from bioconductor, with R etc, then installs R libraries, then uses conda to install snakemake, star, samtools etc etc.
    e.g. add ```RUN R -e ‘BiocManager::install(c(“ggplot2”))’ ``` to install ggplot2.
5. Build the container with: `docker build -t YOUR_DOCKERHUB_NAME/riboseq_pipeline`.
6. Push it to Docker hub with: `docker push YOUR_DOCKERHUB_NAME/riboseq_pipeline`.
7. If necessary delete myproject/pipeline/.snakemake/singularity to make snakemake refresh the container
8. either 
    a) step into the container with `singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://YOUR_DOCKERHUB_NAME/riboseq_pipeline` and run snakemake
    b) use the snake_job script, which passes the flags `—use-singularity  --singularity-args "-B /fast/AG_Ohler/:/fast/AG_Ohler/"` to snakemake so that each cluster node uses it.

#### Useful tips for docker environments

- It’s easy to unintentionally push the incorrect Docker environment as you build. Make sure you get the commands right, especially the tag names. You can go on Docker hub and go to tags>latest, and see what commands were run to make the container. Click on each line to see the complete command.
- Often, problems can occur when other entries in your .bashrc change your path. The code in •Installation above should stop this, but check your path if issues persist.


## Features


## Configuration


## Contributing


## Links

[1] https://sylabs.io/guides/3.0/user-guide/installation.html 'installation guide for singularity'


## Licensing
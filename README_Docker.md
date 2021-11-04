# Docker and singularity - a brief overview.

A container is a layer you run programs without that allows you to virtualize everything below the operating system - i.e. you can run programs as if you have installed a totally different set of software using sudo apt-get, conda, whatever (sound complex, but the principle is obvious once you use it). Docker and singularity are both container management tools. Has the major flaw that you need root access to install and use it, but is easy to install on Mac. Singularity doesn’t require root access, and it able to use docker containers. The idea is therefore:

Build a container on some system where you have root access - your MacBook, a workstation.
Make this container available to the cluster - e.g. by uploading it to docker hub. 
Use singularity on the cluster to download and use this container - thus obviating the need to install software locally and prevent crashes.

## To use. 

1) Put the code at the bottom of this file in your .bashrc so that you can run singularity, and not interfere with it by putting conda etc on your path.
2)To run locally on one node you can just step inside the container:

singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://dermotharnett/riboseq_pipeline

(Note the ‘-B’ entry that mounts file paths, so docker can see the folders we need)

This will take a while to download the first time. Now you have snakemake, STAR, the necessary R libraries, etc etc, to run the pipeline. This is also good for debugging, making sure you have the software you want in the container, etc.

## To change the container:
This is a solution for those using a MAC - on a linux machine this might be easier.
1) Install docker (google it)
2) Create a profile on docker hub, or use a shared one if you want.
3) Create a folder, and put the ‘DOCKERFILE’ from the repo in there.
4) Edit it as needed. Note that when editing you should add lines AFTER the existing ones or it will have to rerun everything above, (which it will have to do anyway the first time). Format should be self explanatory. It first loads an existing one from bioconductor, with R etc, then installs R libraries, then uses conda to install snakemake, star, samtools etc etc.
    e.g. add ```RUN R -e ‘BiocManager::install(c(“ggplot2”))’ ``` to install ggplot2.
5) Build the container with : docker build -t YOUR_DOCKERHUB_NAME/riboseq_pipeline .
6) Push it to docker hub with  : docker push YOUR_DOCKERHUB_NAME/riboseq_pipeline
7) If necessary delete myproject/pipeline/.snakemake/singularity to make snakemake refresh the container
8) either 
    a) step into the container with `singularity run -B /fast/AG_Ohler/:/fast/AG_Ohler/ docker://YOUR_DOCKERHUB_NAME/riboseq_pipeline` and run snakemake
    b) use the snake_job script, which passes the flags ```—use-singularity  --singularity-args "-B /fast/AG_Ohler/:/fast/AG_Ohler/" ``` to snakemake so that each cluster node uses it.



### these should be in your path for singularity to work
```
export PATH=$PATH:/usr/bin
export PATH=$PATH:/usr/sbin
if [ $SINGULARITY_NAME ]; then
    echo "in a sing container,removing miniconda and guix entries from path";
    export PATH=$(echo $PATH | tr ':' '\n' | grep  -v '/guix' | grep -v '/miniconda' | tr '\n' ':' )
    echo $PATH
fi
```

Possible problems:
1) It’s kind of easy to push the incorrect docker environment as you build - make sure you get the commands right, especially the tag names. You can go on docker hub and go to tags>latest , and see what commands were run to make the container. Click on each line to see the whole thing.
2) It’s also easy to mess with singularity by having things in your bashrc that change your path, so that e.g. even in the container you are trying to use guix or miniconda somewhere else. The code above should stop this, but check your path if youre having issues.



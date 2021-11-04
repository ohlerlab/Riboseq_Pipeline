FROM bioconductor/bioconductor_docker:RELEASE_3_13


# Installing R packages
RUN   R -e "BiocManager::install(c('BSgenome','GenomicAlignments','GenomicFeatures','GenomicFiles','Rtracklayer','Rsamtools','data.table','DT','foreach','domc','testthat','ggpubr','gridextra','multitaper','rstan','Matrix','magrittr','R.utils','DESeq2','assertthat','biomaRt','limma','ggrepel','dbscan','zeallot','clusterProfiler','tidyverse','devtools'))" 

RUN   R -e "BiocManager::install(c('doMC','knitr','rmarkdown','cowplot'))"

# Installing miniconda
RUN    touch /.condarc \
&&     wget -O Miniconda_installer.sh -c https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh \
&&     /bin/bash Miniconda_installer.sh -bfp /usr/local \
&&     conda config --file /.condarc --add channels defaults \
&&     conda config --file /.condarc --add channels bioconda \
&&     conda config --file /.condarc --add channels conda-forge

# Installing conda dependancies
RUN  conda install -c conda-forge mamba \ 
&&   mamba install -c conda-forge -c bioconda snakemake=6.2.1
RUN mamba install -c conda-forge -c bioconda cutadapt gffread 
RUN mamba install -c bioconda star
RUN mamba install -c bioconda bedtools
RUN mamba install -c bioconda samtools
RUN conda clean --all

RUN apt-get update && apt-get install less
RUN mamba install perl
RUN mamba install -y -c conda-forge -c bioconda tbb=2020.2 salmon=1.4.0
RUN mamba install -y -c conda-forge ipython ipdb pandas pathlib

CMD /bin/bash "$@"

SHELL ["/bin/bash"]


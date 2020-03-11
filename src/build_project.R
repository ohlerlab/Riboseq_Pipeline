# --------------------------------
# Section 1: Setup
# --------------------------------

#load necessary libraries

library(quietly=TRUE, assertthat)
library(quietly=TRUE, htmltools)
library(quietly=TRUE, tools)
library(quietly=TRUE, rmarkdown)
library(quietly=TRUE, knitr)
library(quietly=TRUE, stringr)
library(quietly=TRUE, magrittr)
library(quietly=TRUE, tidyverse)
library(quietly=TRUE, stringr)
library(quietly=TRUE, GenomicRanges)
library(quietly=TRUE, zeallot)
library(quietly=TRUE, pathview)
library(quietly=TRUE, splines)
# library(quietly=TRUE, rseqdata)
library(quietly=TRUE, pathview)
library(quietly=TRUE, clusterProfiler)
library(quietly=TRUE, memoise)
library(quietly=TRUE, data.table)

matches<-dplyr::matches
select<-dplyr::select
addfileinf <-  .%>%{attr(.,'info')=file.info(.);.}
my_cache <- cache_filesystem(here::here('R_cache')%T>%dir.create(showWarnings=F))

library(here)

# library(quietly=TRUE, org.Mm.eg.db)
# library(quietly=TRUE, org.Hs.eg.db)

#this is a useful library if you're running in the console - but it isn't required.
try({library(colorout)}, silent=TRUE)

##read our design file and load necessary parameters
argv=here::here('src/rseq_design.yaml')
design_file <- commandArgs(trail=TRUE)[1]
if((length(design_file)==0) | is.na(design_file) ) design_file<-normalizePath(file.path('src/rseq_design.yaml'),mustWork=TRUE) 
message(paste0("Reading the design file at ",design_file,'\n\n'))
#now load the and meta data
design_list <- yaml::yaml.load_file(design_file)
#and put the parameters from our design file in the global environment
parameters <- design_list$parameters
assert_that(file.exists(parameters$root))
for(param in names(parameters)) {
  message(paste0('Assigning ',param,' = ',parameters[[param]]))
  assign(param,parameters[[param]])
}
message('\n\n')

message('Building project in')
message(root)


message('Caching results in:')


# source(file.path(rmdfold,'build_project_functions.R'))
devtools::load_all(rseq_path)

subanalysesfolder<-normalizePath(root,'subanalyses',mustWork=TRUE)

metadata_file <- file.path(root,'metadata_file.txt')
covariates_file <- file.path(root,'covariates.txt')
genes_of_interest_file <- file.path(parameters$genes_of_interest_file)
#Read Covariates File



covar_from_samparam <- function(samparam,mcols) {
  out <- read_sample_table(samparam)
  stopifnot(all(mcols %in% colnames(out)))

  out%>%
    dplyr::select(sample_id,group,one_of(mcols)) %>%
    group_by(group) %>%
    mutate(rep=seq_len(n())) %>%
    mutate(sample_name=paste0(group,'_',rep)) %>%
    identity

}


# if(!file.exists(covariates_file)) {
  covar_from_samparam(samparam=parameters$sample_param_file,parameters$mcols) %>%

    write_tsv(file.path(root,'covariates.txt'))
# }

if (file.exists(genes_of_interest_file)) {
  genes_of_interest <- read.csv(genes_of_interest_file,header=T,sep='\t')
}

#Read Covariates File
covariates <- read_covariates(covariates_file)
covariates$sample_id <- as.character(covariates$sample_id)


#get the design (subsets/models/contrasts)
subsets <- get_sample_subsets(design_list, covariates)

#relevel the covariates with factor names in our design file
covariates <- relevel_covariates(covariates, design_list)

#create the folder the subanalyses will go in
subanalysisfolder <- file.path(root,'subanalyses')
message(paste0('creating subanalysis folder ',subanalysisfolder))
# system('rm -rf subanalyses')
dir.create(subanalysisfolder,showWarnings = FALSE)

#create subfolders, including copy of the rmd
devtools::load_all(rseq_path)
modeldirs <- write_parameters(subsets, subanalysisfolder, parameters)

# modeldirs %>% map(convert_to_group_model, covariates)
if(testcontrasts){
  #create test dds objects
  ddsfiles <- modeldirs %>%
    map(create_dds_rld, 
    covariates,
    test=T, 
    corenum=8
  )

  #and see if the contrasts work
  contrastfiles <- lapply(modeldirs, calculate_contrasts)
}
if(testrunonly) stop('unset the testrunonly parameter to do the full analysis')

# --------------------------------
# Section 2: Run full analysis
# --------------------------------

source(file.path(rmdfold,'project_prep.R'))

sample_list <- covariates$sample_id

#TODO, fold teh ensembl entrez map into teh feature annotationp
#TODO - do something about the need for entrez mapping at all

if(!exists('feature_annot')) c(ensembl_to_entrez_map,rseqdata_feature_annot,feature_annot) %<-% projmemoise(get_feature_annot)(my_org_eg_db,annotation_gtf=addfileinf(annotation_gtf))

# aln_stats <-
  # read_aln_stats(file.path(pipeline_results_path, 'qc', 'data'), sample_list)

c(feature_counts_data,feature_counts_summary) %<-% get_feature_counts(
  feature_counts_path=file.path(pipeline_results_path, 'feature_counts', 'data'),
  sample_list)

c(counts_list_all,counts_list_prot_coding) %<-% get_counts_lists(
  feature_annot,
  feature_counts_data)

options('mc.cores'=4, knitr.duplicate.label = 'allow')#need this if generating rmd iteratively, as for each contrast

qcrmd<-normalizePath(file.path(rmdfold,'project_qc.Rmd'))
# source(knitr::purl(qcrmd))
make_qc_html <-function(qcrmd,parameters,covariates,counts_list_prot_coding,useBetaPrior){
  rmarkdown::render(qcrmd,output_file = file.path(parameters$root,'project_qc.html'))
}
make_qc_html <- projmemoise(make_qc_html)
make_qc_html(qcrmd,parameters,covariates,counts_list_prot_coding,useBetaPrior)



ExpressionSet(assayData = counts,phenoData=covariates,featureData=features)


# purely(exclude=c('qcrmd','parameters','covariates','counts_list_prot_coding','useBetaPrior'),
#   {)()

# test_that(" lce1f is right",{
#   gtf_gr<-rtracklayer::import(annotation_gtf)


#   testgid = feature_annot%>%filter(gene_name=='LCE1F')%>%.$gene_id
#   testgr<-gtf_gr%>%subset(type=='gene')%>%subset(gene_id==testgid)  


#   bams <- Sys.glob('pipeline/star/data/*/*bam')%>%setNames(.,basename(dirname(.)))


#   testgr<-gtf_gr%>%subset(type=='gene')%>%subset(gene_id==sample(gene_id,1))  

#   bams%>%map(bamsignals::bamCount,testgr)

#   #What. 
#   expect_true((bams%>%map(bamsignals::bamCount,testgr)%>%.['rnaseq_E2'])==0)
#   expect_true((bams%>%map(bamsignals::bamCount,testgr)%>%.['riboseq_E2'])!=0)

#   #let's do this on the sequence level
#   commonriboseq <- readGAlignments(bams[2],param=ScanBamParam(which=testgr))%>%
#     subset(njunc(.)==0)%>%GRanges%>%table%>%
#     sort%>%tail(1)%>%names%>%GRanges%>%
#     getSeq(x=FaFile('pipeline/hg38.fa'),.)

      
#   readswithsameseq<-str_interp('samtools view ${bams[8]} | grep -e  ${as.character(commonriboseq)}')%>%fread
#   readswithsameseq%>%set_colnames(c('name','seqnames','start',''))



# })


#make our fitted deseq2 objects
create_dds_rld <- projmemoise(create_dds_rld)
forget(create_dds_rld)
ddsfiles <- modeldirs %>% map(.%>%{attr(.,'info')=file.info(.);.})%>%
  map(create_dds_rld, 
    counts_list_prot_coding_sub=counts_list_prot_coding,
    test=FALSE,
    corenum=8)

options('mc.cores'=4, knitr.duplicate.label = 'allow')#need this if generating rmd iteratively, as for each contrast

#perform contrasts 
calculate_contrasts<-projmemoise(calculate_contrasts)
contrastfiles <- mclapply(modeldirs%>%map(addfileinf), calculate_contrasts )

contrastfiles %<>% map(flatten_chr)%>%flatten_chr

#now, produce tables
contrasttables <- contrastfiles %>% map(create_contrast_tables)
create_counts_tables<-projmemoise(create_counts_tables)
counttables <- ddsfiles %>% map(addfileinf)%>%
  map(create_counts_tables, counts_list_prot_coding, rfun=DESeq2::vst)

# here a function for rendering single submodel.Rmd's, not sure it's necessary
render_report <- function(modeldir, recreate=FALSE, is_markdowntest=FALSE) {
  file.copy(file.path(rmdfold,"model.Rmd"),
                                            file.path(modeldir,'submodel.Rmd'),
                                            overwrite=TRUE)

  outputfile <- file.path(modeldir,'submodel.html')
  
  #don't render if the html is already there
  if (file.exists(outputfile) & (recreate==FALSE)) return(outputfile)
 
  # copy in the rmd if not already there
  file.copy(file.path(rmdfold,"model.Rmd"),
            file.path(modeldir,'submodel.Rmd'),
            overwrite=TRUE)
  
  modelrmd <- file.path(modeldir,"submodel.Rmd")

  if(is_markdowntest){

    message('testing the markdown with purl')

    params = list(modeldir = modeldir)
    
    
    modeltestscript <- modelrmd%>%{paste0(dirname(.),'/',basename(file_path_sans_ext(.)),'.R')}
    tmpfile=tempfile()
    tmpfile%>%dirname%>%dir.create
    modelrmd%>%purl(out=tmpfile)
    #trim out the top (so params isn't set) and then move to appropriate location
    readLines(tmpfile)%>%
      {.[cumsum(str_detect(.,'^$'))>1]}%>%
      cat(file=modeltestscript,sep='\n')
    base::source(modeltestscript,local=TRUE)
  }else{
    #and render it
    rmarkdown::render(
      modelrmd,
      output_file = outputfile,
      knit_root_dir = modeldir,
      params = list(
        modeldir = modeldir
        )
    )
  }
}

#this currently fails at random for a few of the contrasts, not sure why
#this set up allows you to just re run this last line until all have worked
#I've uncommented this part because the reports are re-created at the end anyway
#modelreports <- lapply(modeldirs,render_report,recreate=TRUE)

if(!exists('is_markdowntest')) is_markdowntest<-FALSE
# is_markdowntest=FALSE
do_DE=T
modelreports <- lapply(modeldirs,render_report,recreate=TRUE,is_markdowntest=is_markdowntest)
traceback()
header <- '---
title: "Differential expression analysis - all data"
author: "Dermot Harnett"
date: "`r format(Sys.time(), \'%d %B, %Y\')`"
output:
  html_document:
    toc: true
    toc_float: true
---
'

rmd <-  unlist(map(modeldirs, function(md) normalizePath(file.path(md,'submodel.Rmd'))))
# rmd <- c(qcrmd,rmd)
# rmd <- c(normalizePath(file.path(rmdfold,'project_qc.Rmd')))
chunks <- paste0("```{r child = '", rmd, "'}\n```\n")
stitchfile <- file.path(root,'all_reports.Rmd')
cat(header,chunks, sep = '\n', file = stitchfile)


rmarkdown::render(stitchfile)

message(normalizePath(stitchfile))

save.image('data/build_project.Rdata')
#rmarkdown::render(stitchfile)


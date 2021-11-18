##This script is for running deseq - but it's probably going to need modifying every time you use it.
##See the comments with NB in them below to see what you might want to change

DO_XTAIL=F
library(tximport)
library(DESeq2)

library(assertthat)
source(here::here('src/Rprofile.R'))
source(here::here('src/functions.R'))

if(!exists('gtf_gr')) gtf_gr<- readRDS(here("data/gtf_gr.rds"))
if(!exists('tx_countdata')) tx_countdata<- readRDS(here("data/tx_countdata.rds"))


#table with gene 2 gname transformation
gid2gname = gtf_gr%>%mcols%>%
  as.data.frame%>%
  distinct(gene_id,gene_name)%>%
  select(gene_id,gene_name)

#Below we load up the sample parameter df and make sure the factor levels
#are in the order they appear in sample_parameter file.
#read in the sample_parameter file from the piepline
get_sampdf <- function(sampfile){
  sampdf = here('src/sample_parameter.csv')%>%
    read_csv
  stopifnot(colnames(tx_countdata$counts)%in%sampdf$sample_id)
  #use the factor ordering in the sample_parameters file - by default it gets
  sampdf <- sampdf%>%
    as.data.frame%>%
    set_rownames(.$sample_id)%>%
    mutate(group = sample_id%>%str_replace('(_)?\\d+$','')%>%str_replace('\\+','pls')%>%str_replace('\\-','neg'))%>%
    mutate_if(is.character,as_factor)
  if(is.null(sampdf$sample_name)) sampdf$sample_name <- sampdf$sample_id
  for(i in 1:ncol(sampdf)){
    colnm = colnames(sampdf)[[i]]
    if(colnm%in%c('sample_id','sample_name')) next
    if(is.factor(sampdf[[i]])){
     message(paste0('factor levels for ',colnm,':'))
     message(paste(collapse=',',levels(sampdf[[i]])) )
    }
  }
  sampdf
}
sampdf <- get_sampdf(here('src/sample_parameter.csv')) 


#read the counts data - this takes a while sometimes
# file.remove(here('data/dds.rds'))
if(!file.exists(here('data/dds.rds'))){
  samples <- colnames(tx_countdata$counts[,])
  #create a granges object with the coordinates in the first fcounts file
  dds <-   DESeqDataSetFromTximport(tx_countdata,
                               colData= sampdf[match(samples,sampdf$sample_id),],
                               design = ~group
  )
  rownames(dds)<-tx_countdata$counts%>%rownames
  rowRanges(dds) <- gtf_gr%>%subset(type=='exon')%>%
    split(.,trimids(.$gene_id))%>%
    .[rownames(dds)]
  rowData(dds) <- rowRanges(dds)%>%{.@unlistData[.@partitioning@end]}%>%
    mcols
  #save deseq object 
  rownames(dds) <- rowData(dds)$gene_id
  rownames(colData(dds)) <- colData(dds)$sample_id
  saveRDS(dds,here('data/dds.rds'))
}else{
  dds <-readRDS(here('data/dds.rds'))
}
cat('successfully loaded  countdata and annotation')
stopifnot(exists('sampdf'))

#get normalized counts.
#you might want to use vst if there's a lot of data, for speed
dir.create(here('data'))
# file.remove('data/normcounts.rds')
if(!file.exists(here('data/normcounts.rds'))){
  # 
  normfunc <- if(ncol(dds)>20) DESeq2::vst else DESeq2::rlog
  # 
  normcounts <- normfunc(dds)
  rownames(normcounts) <- rownames(dds)
  saveRDS(normcounts,here('data/normcounts.rds'))
}else{
  normcounts<-readRDS(here('data/normcounts.rds'))
}
cat('successfully normalized  countdata')

stopifnot(exists('dds'))
stopifnot(exists('normcounts'))
stopifnot(exists('sampdf'))

{
  #do a trial run of DESeq to get the contrast names
  #NB - modify this formula to match you covariates - genotype, timepoint, treatment, whatever.
  design <- as.formula('~ fraction*genotype')
  design(dds) = design
  testdds = DESeq(head(dds,2*ncol(dds)))
  resnames = resultsNames(testdds)
  resnames
}


################################################################################
########COnstruct Model
################################################################################

  #NB - after doing your formula, look at the value of 'resnames' and choose
#which contrasts you want (or specify a combination of them as a binary vector)
resnames = resultsNames(testdds)
#
contrasts <- list("fraction_membrane_vs_cyto", "genotype_KO_vs_WT",
"fractionmembrane.genotypeKO")%>%setNames(.,.)
#
#check the contrasts
for(contrastname in names(contrasts)){
  contrast = contrasts[[contrastname]]
  if(is.character(contrast) & length(contrast)==1){
    assert_that(contrast %in% resnames)
    contrast <- as.numeric(resnames==contrast)
    contrasts[[contrastname]] <- contrast
  }
  message(contrastname)
  try({
    results(testdds,contrast)
  })
}

cat('tested contrasts')



dir.create(here('data'))
# file.remove('data/resultslist.rds')
if(!file.exists(here('data/resultslist.rds'))){
  design(dds) <- design
 
  dds <- DESeq(dds,betaPrior = F)
  #and get the contrasts we need
  resultslist <- lapply(contrasts,results,object = dds)
  names(resultslist) <- names(contrasts)
  for(i in seq_along(resultslist))resultslist[[i]]$gene_id = rowData(dds)$gene_id
  saveRDS(resultslist,here('data/resultslist.rds'))

  ddstodf<- function(ddsdf)results(ddsLRT)%>%as.data.frame%>%rownames_to_column('gene_id')
  ddsLRT = DESeq(dds, test="LRT", full=design(dds),reduced = ~ 1)
  changegenes <- results(ddsLRT)%>%ddstodf %>%subset%>%filter(padj < 0.05)%>%.$gene_id
  
  saveRDS(changegenes,here('data/changegenes.rds'))
  
  if(DO_XTAIL){
    run_xtail = function()Sys.glob(here('pipeline','xtail/*'))%>%setNames(.,basename(.)%>%str_extract(regex('(?<=xtail_).*(?=\\.tsv)')))
#    xtailcontrasts = 'condition'
 #   xtailfiles = run_xtail(dds)  
    xtailfiles <- run_xtail()
    xtailres = map_df(.id='contrast',xtailfiles,fread)
    xtailres$isxtail = TRUE
    xtailres$contrast = paste0('xtailTE_',xtailres$contrast)
    
    ddsreslist<-resultslist[!str_detect(names(resultslist),'xtail')]
    # xtailres%<>%rename(log2FoldChange := log2fc)
    xtailres%<>%rename(gene_id := feature_id)
    xtailres %<>% rename(baseMean := base_mean)
    resultslist <- c(ddsreslist,split(xtailres,xtailres$contrast))  
  }
}else{
  resultslist<-readRDS(here('data/resultslist.rds'))
  changegenes<-readRDS(here('data/changegenes.rds'))
}

cat('successfully calculated contrasts:')
message(names(resultslist)%>%paste(collapse=' '))

columns <- c("constrast", "baseMean", "log2FoldChange", "lfcSE", "stat",
"pvalue", "padj", "gene_id")

dir.create(here('tables'))
resultslist%>%map_df(.id='constrast',.%>%as.data.frame%>%set_rownames(NULL))%>%
  select(one_of(colvals))%>%
  write_tsv(here('tables/all_log2fc.tsv'))



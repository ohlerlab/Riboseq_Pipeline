#' # Clustering - Zea Mays Meoisis Riboseq data

#+ setup, include=FALSE, echo=FALSE, eval=T
{
knitr::opts_chunk$set(root.dir = here::here(),eval=FALSE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=FALSE)

library(rmarkdown)
library(knitr)
library(here)
library(magrittr)
library(stringr)
library(ggplot2)
library(dbscan)

source('../cortexomics/src/R/Rprofile.R')

isknitr<-isTRUE(getOption('knitr.in.progress'))
# is_sourcecall <- 
if(!isknitr){
	dir.create(here('Reports'),showWarnings=F)
	projdir <- here()
	sourcecall <- sys.calls()[1]%>%as.character
	filename = sourcecall%>%str_extract(str_interp('(?<=${projdir}/).*\\.[rR](?=["\'])'))
	filehtml <- filename%>%str_replace('\\.[rR]','\\.html')%>%basename%>%here('Reports',.)
	message(paste0('got filename from call stack, spinning: ',filename))
	rmarkdown::render(knitr::spin(here(filename),knit=F),output_file=filehtml,output_dir=here('Reports'),knit_root_dir=here())
	stop(paste0('stopping, knit spin complete:',normalizePath(filehtml)))
}

}

#+ setup2
library(rseq)
library(memoise)
library(here)
library(stringr)
library(purrr)
library(tidyverse)
library(data.table)
library(magrittr)
devtools::load_all('/fast/work/groups/ag_ohler/dharnet_m/senescence_ribo/Applications/rseq/')

gtf <- 'pipeline/genes.chrnm.gtf'
anno <- rtracklayer::import(gtf)

counts <- Sys.glob('pipeline/feature_counts/data/*/feature_counts')%>%map(fread)
counts <- counts%>%map(~.[,c(1,7)])%>%reduce(left_join,by='Geneid')%>%{colnames(.)%<>%dirname%>%basename%>%{c('gene_id',.[-1])};.}
counts <- set_rownames(counts[,-1],counts[[1]])

anno <- mcols(anno)[match(rownames(counts),anno$gene_id),c('gene_id','gene_name')]%>%as.data.frame

sampledata <- fread('src/sample_parameter.csv')%>%select(-(library_layout:fragment_length_sd))%>%as.data.frame%>%set_rownames(.,.$sample_id)

colnames(counts)==rownames(sampledata)

countdata <- ExpressionSet(assaydata=AnnotatedDataFrame(counts),phenoData=AnnotatedDataFrame(sampledata)
	#,featureData=AnnotatedDataFrame(anno)
)

DESeqDataSet

library(GenomicRanges)


{
if(!exists(here('data/ddsall.rds'))){
	dds <-  create_dds(counts_list_prot_coding,
             num_cores = 8,
             design = ~ assay ,
             useBetaPrior = useBetaPrior)
	saveRDS(dds,here('data/ddsall.rds'))
}else{
	dds<-readRDS(here('data/ddsall.rds'))
}


if(!exists(here('data/vstall.rds'))){
	vst <- DESeq2::vst(dds, blind = TRUE)
	saveRDS(vst,here('data/vstall.rds'))
}else{
	vst<-readRDS(here('data/vstall.rds'))
}

assay(vst)%>%rowMeans%>%{
	qplot(data=data.frame(x=.),x=x,geom='density',fill=I('blue'))+theme_bw()+
	ggtitle('Normalized Expression Density plot')
}

clustsizes <- c(5,10,20,50)

dbscans <- mclapply(clustsizes,function(clustsize){ dbscan::hdbscan(assay(vst),
	minPts = clustsize )})

if(!exists(here('data/dbscans.rds'))){
	clustsizes <- c(5,10,20,50)
	dbscans <- mclapply(clustsizes,function(clustsize){ dbscan::hdbscan(assay(vst),
		minPts = clustsize )})
	dbscans%<>%setNames(paste0('mp',clustsizes))
	saveRDS(dbscans,here('data/dbscans.rds'))
}else{
	dbscans<-readRDS(here('data/dbscans.rds'))
}

dbscanmaps <- dbscans%>%map(.%>%.$cluster%>%safe_hashmap(keys=rownames(dds),values=.) )

clustermap <- dbscanmaps[[4]]
clustervals <- clustermap$values()%>%unique
clusterval_i = clustervals[1]
results.tab <- rownames(dds)%>%setNames(is_in(clustermap[[.]],clusterval_i),.)%>%rungo(GTOGO,'BP')
myontology <- 'BP'
clustgores<-clustervals%>%mclapply(function(clusterval_i){
	results.tab <- rownames(dds)%>%setNames(is_in(clustermap[[.]],clusterval_i),.)%>%rungo(GTOGO,myontology)
})

#For this plot we'll take the top N terms in each cluster, then plot cluster on the x axis
#term on the y axis, 
# plot_go_enrich(results.tab,sort_var = 'elimFisher',"Gabalike_Gliaresp")

mapname = 'dbscan_mp50'
names(clustgores)<-paste0('cluster_',seq_along(clustgores))

go_comparison_plot <- clustgores%>%map_df(.id='clust',.%>%arrange(elimFisher)%>%head(10)%>%select(Term,elimFisher))%>%
	mutate(Term = as_factor(Term))%>%
	mutate(cluster = as_factor(cluster))%>%
	ggplot(.,aes(x=cluster,y=Term,color=-log10(elimFisher),size=-log10(elimFisher)))+
	geom_point()+
	ggtitle(mapname)+
	theme_bw()
}

#' ## subheading
#' 
#' 
#' 





stopifnot(project_cache=='/fast/work/groups/ag_ohler/dharnet_m/Mitosis_Riboseq/R_cache')
create_counts_tables<-projmemoise(create_counts_tables)
counttables <- ddsfiles %>% map(addfileinf)%>%
  map(create_counts_tables, counts_list_prot_coding, rfun=DESeq2::vst)



#+ aplot, fig.width =12,fig.height=4,out.width=1200,out.height=450,dev='svg',include=TRUE
qplot(1)
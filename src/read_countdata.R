# if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
library(tximport)
source(here::here('src/Rprofile.R'))
source(here::here('src/functions.R'))


gtf <- here(paste0('pipeline/',basename(yaml::yaml.load_file(here('config/config.yaml'))$GTF_orig)))
stopifnot(file.exists(gtf))
if(!exists('gtf_gr')) gtf_gr <- rtracklayer::import(gtf)


sampleinfo <- read_csv(here('config/sample_parameter.csv'))
ribosamples <- sampleinfo%>%filter(isriboseq)%>%.$sample_id
rnasamples <- sampleinfo%>%filter(!isriboseq)%>%.$sample_id
salmonfiles <- ''[0]
if(length(rnasamples)>0) salmonfiles <- paste0('salmon/data/',rnasamples,'/quant.sf')
if(length(salmonfiles)>0) stopifnot(all(file.exists(salmonfiles)))
ribostanfiles <- here(paste0('pipeline/ribostan/',ribosamples,'/',ribosamples,'.ribostan.tsv'))
if(length(ribostanfiles)>0) stopifnot(all(file.exists(ribostanfiles)))
countdatafiles <- c(salmonfiles,ribostanfiles)
#

################################################################################
########Now load annotation data
################################################################################
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
tx2genetbl = gtf_gr%>%mcols%>%as.data.frame%>%filter(!is.na(gene_id),!is.na(transcript_id))%>%
	distinct(transcript_id,gene_id)%>%select(transcript_id,gene_id)
trid2gid=tx2genetbl%>%{setNames(.$gene_id,.$transcript_id)}
#
ribostantrs = ribostanfiles[1]%>%read_tsv%>%.$Name
trs = ribostantrs
if(length(salmonfiles)>0){
	salmontrs =  salmonfiles[1]%>%read_tsv%>%.$Name%>%str_extract('[^|]+')
	salmoncds = salmontrs%>%intersect(gtf_gr%>%subset(type=='CDS')%>%.$transcript_id)
	# cdswidth = gtf_gr%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%width%>%sum
	# cdsrangedf <- tibble(Name=salmoncds,annolength = cdswidth[salmoncds])
	# cdsrangedf$Name%<>%trimids
	# trs = cdsrangedf$Name
	trs = intersect(ribostantrs,salmontrs)
}
trid2gid = trid2gid[trs] 
#we have to trim the ids for orfquant
tx2genetbl$transcript_id%<>%trimids
tx2genetbl$gene_id%<>%trimids

################################################################################
########Collect transcript level info for top transcripts
################################################################################

countdatafiles%<>%setNames(.,basename(dirname(.)))
tx_countdata = tximport(files=countdatafiles,
	ignoreTxVersion=TRUE,
	tx2gene=tx2genetbl,
	type='salmon',
	# countsFromAbundance='scaledTPM',
	ignoreAfterBar=TRUE,
	# importer=filereadfunc
)

iso_tx_countdata = 	tximport(files=countdatafiles,
	txOut=TRUE,
	ignoreTxVersion=TRUE,
	tx2gene=tx2genetbl,
	type='salmon',
	countsFromAbundance='scaledTPM'
)

stopifnot(iso_tx_countdata$abundance%>%rownames%>%setequal(trs))

dir.create(here('data'))
tx_countdata%>%saveRDS(here('data/tx_countdata.rds'))
iso_tx_countdata%>%saveRDS(here('data/iso_tx_countdata.rds'))

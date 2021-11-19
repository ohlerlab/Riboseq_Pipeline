
################################################################################
########Metagene plots and positions specific analysis
################################################################################
HIGHCOUNTLIM=32

{
# source('src/0_load_annotation.R')
library(here)
library(GenomicAlignments)
library(GenomicFeatures)
library(abind)
library(Biostrings)
library(assertthat)
library(magrittr)
library(data.table)
library(tidyverse)
select <- dplyr::select
slice <- dplyr::slice
#
fafile <- here(paste0('pipeline/',basename(yaml::yaml.load_file(here('config/config.yaml'))$REF)))
fafile%<>%str_replace('.gz$','')
fafileob <- FaFile(fafile)

gtf <- here(paste0('pipeline/',basename(yaml::yaml.load_file(here('config/config.yaml'))$GTF)))
if(!exists('gtf_gr')) gtf_gr <- rtracklayer::import(gtf)

source(here('src/functions.R'))
#
roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}


STARTCDSSIZE = 60
STOPCDSSIZE = 60
FPEXT = 36
STARTWINDSIZE = STARTCDSSIZE+FPEXT
TPUTREXT = 36
STOPWINDSIZE = STOPCDSSIZE+TPUTREXT

MINCDSSIZE = STARTCDSSIZE+STOPCDSSIZE+1

TOTBINS = (STARTWINDSIZE)+(STOPWINDSIZE)+1
}


# displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
# stageconv = names(displaystagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
# stagecols <- displaystagecols%>%setNames(names(stageconv))

CODOONSOFINTEREST<-DNAStringSet(c('TCA','TCG','TCC','TCT'))
codons2scan <- c(CODOONSOFINTEREST,reverseComplement(CODOONSOFINTEREST))
# codons2scan <- c(CODOONSOFINTEREST,(CODOONSOFINTEREST))
ctrlcodons<-c('TTT','GGG','AAA','CCC')
codons2scan <- c(CODOONSOFINTEREST,ctrlcodons)
#codons2scan <- CODOONSOFINTEREST

FLANKCODS<-15

gn2tr = mcols(gtf_gr)%>%as.data.frame%>%
	filter(!is.na(transcript_id), !is.na(gene_id))%>%
	distinct(transcript_id,gene_id)
cdsgrl = gtf_gr%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%sort_grl_st
gn2tr = gn2tr%>%left_join(cdsgrl%>%width%>%sum%>%enframe('transcript_id','length'))
splocs = read_tsv(here("ext_data/sp_list.txt"))
splocs$transcript_id = splocs$Transcript.stable.ID.version
tmlocs = read_tsv(here("ext_data/tm_list.txt"))
tmlocs$transcript_id = tmlocs$Transcript.stable.ID.version

#
if(!file.exists(here('data/longesttrs.rds'))){

	#code to get the highest abundance tr per gene
	# longesttrs = iso_tx_countdata$abundance[,TRUE]%>%
	# 	as.data.frame%>%
	# 	rownames_to_column('transcript_id')%>%
	# 	left_join(tx2genemap%>%set_colnames(c('transcript_id','gene_id')))%>%
	# 	pivot_longer(-one_of('gene_id','transcript_id'))%>%
	# 	# separate(name,c('time','assay','rep'))%>%
	# 	group_by(gene_id,transcript_id)%>%summarise(value=median(value))%>%
	# 	group_by(gene_id)%>%
	# 	filter(gene_id %in% highcountgenes)%>%
	# 	slice(which.max(value))%>%.$transcript_id
	longesttrs = gn2tr%>%group_by(gene_id)%>%
		mutate(hastm = transcript_id%in%tmlocs$transcript_id)%>%
		mutate(hassp = transcript_id%in%splocs$transcript_id)%>%
		arrange(-hastm, -hassp, -length)%>%
		slice(1)%>%
		filter(!is.na(length))%>%
		.$transcript_id
	saveRDS(longesttrs,here('data/longesttrs.rds'))
}else{
	longesttrs<-readRDS(here('data/longesttrs.rds'))
}

trid2gid = gn2tr$gene_id%>%setNames(gn2tr$transcript_id)
gid2trid = gn2tr$transcript_id%>%setNames(gn2tr$gene_id)

cds2use <- cdsgrl[longesttrs]

cdsseq <- cds2use%>%{GenomicFeatures::extractTranscriptSeqs(.,x=fafileob)}
allcodons=getGeneticCode()

#get genes with reasonably high counts
counts <- readRDS(here('data/tx_countdata.rds'))
counts<-counts$counts%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	pivot_longer(-gene_id,names_to='sample',values_to='count')
highcountgenes = counts%>%group_by(gene_id)%>%filter(any(count>HIGHCOUNTLIM))%>%.$gene_id%>%unique
highcountcovtrs = longesttrs[trid2gid[longesttrs]%>%str_replace('\\.\\d+$','')%>%is_in(highcountgenes)]

exonsgrl<-gtf_gr%>%subset(type=='exon')%>%split(.,.$transcript_id)
trlens = exonsgrl%>%width%>%sum
trlensgr = trlens[longesttrs]%>%enframe('seqnames','end')%>%mutate(start=1)%>%GRanges

# offsets <- read_tsv('ext_data/offsets_manual.tsv')

trseqinfo = Seqinfo(seqnames=longesttrs,seqlengths=trlens[longesttrs])

cdsstarts = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1)%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
cdsends = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1,'end')%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
trcds = GRanges(names(cdsstarts),IRanges(cdsstarts,cdsends))


# tx_countdata$counts%>%colSums%>%divide_by(1e6)
#riboslevels <- seqlevels(ribogr)

ribobams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(neg=T,'rnaseq')
trimids <- .%>% str_replace('\\.\\d+$','')
# ribobams <- ribobams%>%str_subset('4E')
names(ribobams) <- ribobams%>%basename%>%str_replace('_\\d+.bam','')%>%str_replace('rep','')%>%str_replace('riboseq','ribo')%>%str_replace('\\+','pos')%>%str_replace('\\-','neg')

#filter for only long enough ones, if we're doing windows
trspacecds = pmapToTranscripts(cdsgrl[highcountcovtrs],exonsgrl[highcountcovtrs])
trspacecds%<>%unlist
ltrspacecds = trspacecds
longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
ltrspacecds = ltrspacecds[longcdstrs]

ribobam <-ribobams[1]
if(!file.exists(here('data/fpcovlist.rds'))){
	fpcovlist = ribobams%>%mclapply(mc.cores=10,function(ribobam){
		ribogr <- GenomicAlignments::readGAlignments(ribobam)
		mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
		ribogr%<>%as("GenomicRanges")
		cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl[highcountcovtrs])
		cov$readlen = mcols(ribogr)$readlen[cov$xHits]
		cov%<>%subset(between(readlen,25,35))
		message('warning shifting fps')
		strand(cov)='+'
		cov%<>%resize(1,'start')%>%GenomicRanges::shift(12)
		# cdscov <- cov%>%mapToTranscripts(trspacecds[highcountcovtrs])
		# c
		# cdscov$readlen <- cov$readlen[cdscov$xHits]
		# cdscov$xHits <- NULL
		# cdscov$transcriptsHits<-NULL
		cdscov <- cov	
		split(cdscov,cdscov$readlen)%>%
			lapply(coverage)
	})
	names(fpcovlist)<-names(ribobams)
	saveRDS(fpcovlist,here('data/fpcovlist.rds'))
}else{
	fpcovlist<-readRDS(here('data/fpcovlist.rds'))
}

sum(runLength(fpcovlist[[1]][[1]]))[ttrs]==seqlengths(metasrpwindows)[ttrs]
sum(runLength(fpcovlist[[1]][[1]]))[ttrs]==sum(width(cdsgrl))[ttrs]
sum(runLength(fpcovlist[[1]][[1]]))[ttrs]==sum(width(exonsgrl))[ttrs]

# #see how many we eliminate for length reasons
# message(length(ltrspacecds))
# #spand the cds as necessary

# is_offchr<-function(gr,si){
#   if(is(gr,'GenomicRangesList')){
#    (end(gr) > split(seqlengths(gr)[as.character(unlist(seqnames(gr)))],gr@partitioning) ) %in% TRUE
#   }else{
#     seqinfo(gr)<-si
#     end(gr) > seqlengths(gr)[as.character(seqnames(gr))]

#   }
# }
# is_out_of_bounds <- function(gr,si = seqinfo(gr)){
#   start(gr)<1 | is_offchr(gr,si) 
# }


# regScoreSums<-function(srle,gr){
# 	scoregr=srle%>%as("GRanges")%>%subset(score!=0)%>%width1grs
# 	ov = findOverlaps(scoregr,gr)
# 	sumscores = tibble(reg=ov@to,score=scoregr$score[ov@from])%>%
# 		group_by(reg)%>%summarise_at('score',sum)
# 	tibble(reg=1:length(gr))%>%left_join(sumscores, by='reg')%>%
# 		mutate(score=replace_na(score,0))%>%
# 		pluck('score')%>%
# 		setNames(names(gr))
# }
# width1grs <- function(gr){
# 	stopifnot(Negate(is.unsorted)(gr))
# 	isw1 <- width(gr)==1
# 	broad <- gr[!isw1]
# 	#vector of integers - 1,2,3 for range 1-3
# 	narrowstarts <- unlist(as(broad@ranges,'IntegerList'))
# 	narrow <- {GRanges(
# 			rep(seqnames(broad),width(broad)),
# 			IRanges(narrowstarts,w=1)
# 		)}
# 	mcols(narrow) <- mcols(broad)[rep(seq_along(broad),width(broad)),,drop=F]
# 	sort(c(gr[isw1],narrow))
# }

# sampnames<-names(fpcovlist)
# rls<-names(fpcovlist[[1]])


# fputrs <- gaps(ltrspacecds)%>%subset(start==1)%>%subset(strand=='+')%>%setNames(seqnames(.))
# tputrs <- gaps(ltrspacecds)%>%subset(start!=1)%>%subset(strand=='+')%>%setNames(seqnames(.))


# tr_regions<-list(fputr=fputrs,cds=ltrspacecds,tputrs=tputrs)
# scores <- map_df(.id='region',tr_regions,function(tr_region){
# 	lengths<-tr_region%>%width%>%setNames(.,names(tr_region))%>%enframe('tr_id','length')
# 	out = map_df(.id='sample',sampnames%>%setNames(.,.),function(sampname){
# 		map_df(.id='readlen',rls%>%setNames(.,.),function(rl){
# 		scoregr <- fpcovlist[[sampname]][[rl]]%>%
# 			GRanges%>%
# 			subset(score!=0)%>%
# 			width1grs%>%
# 			GenomicRanges::shift(12)%>%
# 			identity
# 		ov = findOverlaps(scoregr,tr_region)
# 		sumscores = tibble(tr_region=ov@to,score=scoregr$score[ov@from])%>%
# 				group_by(tr_region)%>%summarise_at('score',sum)
# 		tibble(tr_region=1:length(tr_region))%>%left_join(sumscores, by='tr_region')%>%
# 				mutate(score=replace_na(score,0))%>%
# 				pluck('score')%>%
# 				setNames(names(tr_region))%>%
# 				enframe('tr_id','score')
# 		})
# 	})
# 	out = out%>%left_join(lengths)
# 	out
# })
# rlsumscores <- scores%>%group_by(region,sample,tr_id,length)%>%summarise(score=sum(score))%>%
# 	mutate(density = score/length)

# dir.create('tables')
# rlsumscores%>%write_tsv('tables/region_scores.tsv')

# rlsumscores$sample%>%unique
# l
# stop()

# #create our Drimseq objects
# iso_tx_countdata <- readRDS('data/iso_tx_countdata.rds')
# cts <- iso_tx_countdata$counts
# cts <- cts[rowSums(cts) > 0,]
# #
# # tx2genemap%<>%set_colnames(c('tr_id','g_id'))
# #
# #
# rlsumscores$gene_id <-  trid2gid[rlsumscores$tr_id]
# rlsumscores <- rlsumscores%>%mutate(feature_id = paste0(gene_id,'_',region))
# counts <- rlsumscores%>%ungroup%>%select(gene_id,feature_id,sample,score)%>%
# 	pivot_wider(names_from='sample',values_from='score')
# library(DRIMSeq)
# #
# library(DRIMSeq)
# subunit = '4E'
# subunits <- colnames(iso_tx_countdata$counts)%>%str_split('_')%>%map_chr(1)%>%unique%>%setdiff('rnaseq')
# #
# drimseqlist <- lapply(subunits,function(subunit){
# 	subunit_counts <- counts%>%select(gene_id,feature_id,matches(str_interp('^${subunit}')))
# 	allcountdesign = subunit_counts%>%colnames%>%tail(-2)%>%unique%>%data.frame(sample=.)%>%separate(sample,into=c('subunit','induced','rep'),remove=F)
# 	allcountdesign%<>%mutate(sample_id=sample)
# 	d <- dmDSdata(counts=subunit_counts%>%as.data.frame, samples=allcountdesign%>%as.data.frame)
# 	n.small = allcountdesign%>%group_by(induced)%>%tally%>%.$n%>%min
# 	n = nrow(allcountdesign)
# 	d <- dmFilter(d,
# 	              min_samps_feature_expr=n.small, min_feature_expr=10,
# 	              min_samps_feature_prop=n.small, min_feature_prop=0.1,
# 	              min_samps_gene_expr=n, min_gene_expr=10)
# 	design_full <- model.matrix(~induced, data=DRIMSeq::samples(d))
# 	#Run drimseq
# 	d = d
# 	d <- dmPrecision(d, design=design_full)
# 	d <- dmFit(d, design=design_full)
# 	d <- dmTest(d)
# 	d
# })

# drimseqlist%>%saveRDS(here('data/drimseqlist.rds'))



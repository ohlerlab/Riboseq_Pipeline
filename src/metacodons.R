source('src/load_pos_data.R')

################################################################################
########Verify offsets with metacodon plots
################################################################################
STARTBUFF=60
ENDBUFF=60
exonseq = exonsgrl[highcountcovtrs]%>%extractTranscriptSeqs(x=fafileob)
i=1
#get allcodlist granges object descxribing codon positions in the transcripts
if(!file.exists(here('data/allcodlist.rds'))){
	allcodlist <- lapply(seq_along(allcodons)%>%setNames(names(allcodons)),function(i){
		#	
		codon=names(allcodons)[[i]]
		message(codon)
		codmatches<-vmatchPattern(pattern=codon,exonseq[highcountcovtrs])#exclude the start ccodon
		#
		nmatches = 	codmatches%>%elementNROWS 
		#
		matchgr<-codmatches%>%unlist%>%GRanges(names(.),.)
		matchgr$cdspos = start(matchgr) - start(trspacecds[as.vector(seqnames(matchgr))])
		matchgr%<>%subset(cdspos %%3 == 0)
		seqlengths(matchgr) = exonsgrl%>%width%>%.[seqlevels(matchgr)]%>%sum
		innercds = trspacecds%>%subset(width>(3+STARTBUFF+ENDBUFF))%>%
			resize(width(.)-STARTBUFF,'end')%>%
			resize(width(.)-ENDBUFF,'start')
		matchgr = matchgr%>%subsetByOverlaps(innercds)
		codmatchwindows<-matchgr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')
		codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
		codmatchwindows%<>%subsetByOverlaps(innercds)
		codmatchwindows
	})
	allcodlist=allcodlist%>%GRangesList%>%unlist
	saveRDS(allcodlist,here('data/allcodlist.rds'))
}else{
	allcodlist<-readRDS(here('data/allcodlist.rds'))
	stopifnot(allcodlist@seqinfo@seqnames%>%setequal(highcountcovtrs))
	stopifnot(allcodlist@seqinfo@seqlengths%>%setequal(sum(width(exonsgrl[highcountcovtrs]))))
}



# fpcovlist <- fpcovlist[names(allbamtbls)]
# if(!file.exists(here('data/fpprofilelist.rds'))){
# 	fpprofilelist <-imap(fpcovlist[mainsamps],function(sampfpcov,sampname){
# 		trsums = sampfpcov%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
# 		sampfpcov%>%lapply(function(rlfpcov){
# 			rlfpcov = rlfpcov/(trsums)
# 			allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
# 			cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
# 			message('.')
# 			out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
# 			out
# 		})
# 	})
# 	saveRDS(fpprofilelist,here('data/fpprofilelist.rds'))
# }else{
# 	fpprofilelist<-readRDS(here('data/fpprofilelist.rds'))
# }

# codonprofiledat = fpprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
# codonprofiledat%<>%mutate(position = position - 1 - (FLANKCODS*3))
# codonprofiledat%<>%group_by(sample,readlen,codon)%>%mutate(count= count / median(count))
# codonprofiledat%<>%filter(!codon %in% c('TAG','TAA','TGA'))

################################################################################
########Now get the profiles of the codons
################################################################################
	
if(!file.exists(here('data/fprustprofilelist.rds'))){
	fprustprofilelist <-mapply(fpcovlist,names(fpcovlist),F=function(sampfpcov,sampname){
		trsums = sampfpcov%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
		sampfpcov%>%lapply(function(rlfpcov){
			rlfpcov = rlfpcov > mean(rlfpcov)
			allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
			cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
			message('.')
			out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
			out
		})
	})
	saveRDS(fprustprofilelist,here('data/fprustprofilelist.rds'))
}else{
	fprustprofilelist<-readRDS(here('data/fprustprofilelist.rds'))
}

rustprofiledat = fprustprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
rustprofiledat%<>%mutate(position = position - 1 - (FLANKCODS*3))
rustprofiledat%<>%group_by(sample,readlen,codon)%>%mutate(count= count / median(count))
rustprofiledat%<>%filter(!codon %in% c('TAG','TAA','TGA'))


################################################################################
########testing the pca based a-site calls
################################################################################
offsets = tibble(readlen=paste0(28:32),offset=12)

{
plotfile <- 'plots/fppos_vs_codon_variance.pdf'
plotfile%>%dirname%>%dir.create(showWarn=F)
# offsets%<>%mutate(readlen=paste0(length))
pdf(plotfile,w=12,h=12)
#plotting variance amongst codons at each point.
# codonprofiledat%>%
rustprofiledat%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	# filter(count==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	# filter(sample%>%str_detect('ribo_1'))%>%
	filter(!is.nan(count))%>%
	summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
	# separate(sample,c('time','assay','rep'))%>%
	group_by(readlen,sample,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	filter(position> -numreadlen+6,position < -6)%>%
	filter(numreadlen>=27,numreadlen<=32)%>%
	filter(!sample%>%str_detect('ribo_4E_posIAA_3'))%>%
	arrange(position)%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~sample,scale='free')+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}


newcodonoccs = rustprofiledat%>%left_join(offsets)%>%filter(position == -offset-3)%>%
	# group_by(sample)%>%group_slice(1)
	identity
newcodonoccs%<>%separate(sample,c('time','assay','rep'))
newcodonoccs%<>%group_by(codon,time)%>%summarise(count=mean(count))
tRNAstats <- read_tsv('tables/tRNA_stat_df')
tRNAstats%<>%filter(fraction=='total')

tRNAstats%>%select(codon,time,dwell_time)%>%left_join(newcodonoccs%>%select(codon,time,count))%>%
	{quicktest(.$dwell_time,.$count)}


codonprofiledat = fpprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
#
windowoffsets = fpprofilelist[[1]][[1]][[1]]%>%length%>%seq(1,.)%>%subtract(1+(FLANKCODS*3))%>%multiply_by(-1)
windowpos = fpprofilelist[[1]][[1]][[1]]%>%length%>%seq(1,.)%>%subtract(1+(FLANKCODS*3))
# codonprofiledat$position = windowpos[codonprofiledat$position]

#
profvarpca = codonprofiledat%>%
	split(.,.$sample)%>%
	map_df(.id='sample',.%>%
		split(.,list(.$readlen))%>%
		# .[[1]]%>%
		map_df( .id='readlen',.%>%
			mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
			# filter(position> -numreadlen-1,position < 1)%>%
			filter(position> -numreadlen+6,position < -6)%>%
			# filter(sample=='E13_ribo_1')%>%
			ungroup%>%
			select(-numreadlen,-readlen,-sample)%>%
			# group_by(codon)%>%mutate(count/median(count))%>%
			spread(position,count)%>%
			{set_rownames(.[,-1],.$codon)}%>%
			princomp%>%{.$loadings[,1]}%>%{./.[which.max(abs(.))]}%>%enframe('position','pca1')
		)
	)
#
profvarpca%<>%select(sample,readlen,position,pca1)
profvarpca$readlen = paste0('rl',profvarpca$readlen)
library(rlang)
offsets <- read_tsv('ext_data/offsets_manual.tsv')
#

plotfile<-'plots/figures/figure2/trna_codons/fppos_vs_codon_pcascore.pdf'
offsets%<>%mutate(readlen=paste0('rl',length))
pdf(plotfile,w=12,h=12)
profvarpca%>%
	# slice_by(sample,c(1,2,3,4,5,6))%>%
	filter(sample%>%str_detect('ribo_1'))%>%
ggplot(data=.,aes(y=pca1,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)




pcaderrivedpos = profvarpca%>%
	# mutate(pca13wind = pca1+lead(pca1)+lead(pca1,2))%>%
	mutate(pca13wind = pca1)%>%
	group_by(readlen,position)%>%summarise(pca13wind=mean(pca13wind))%>%
	group_by(readlen)%>%slice(which.max(pca13wind))
# pcaderrivedpos=pcaderrivedpos%>%mutate(position = map(as.numeric(position),~ (.:(.+2))))%>%unnest

codonprofiles_pcawind = codonprofiles%>%semi_join(pcaderrivedpos%>%mutate(position=as.numeric(position)))
codonprofiles_pcawind = codonprofiles_pcawind%>%group_by(sample,codon,fraction)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(codon,time)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))
	
codonprofiles_pcawind%>%inner_join(timecodon_tAb)%>%
	filter(.,time=='E13')%>%{
	quicktest(.$abundance,.$occ_nonorm)
}













offsets = tibble(readlen=(25:35),offset=12)

################################################################################
########
################################################################################

if(!file.exists(here('data/fpcovlist.rds'))){
	psitecovlist = ribobams%>%mclapply(mc.cores=4,function(ribobam){
		ribogr <- GenomicAlignments::readGAlignments(ribobam)
		mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
		ribogr%<>%as("GenomicRanges")
		ribogr%<>%as.data.table 
		ribogr$seqnames%<>%str_extract('[^|]+')
		ribogr%<>%select(transcript_id=seqnames,start,readlen)
		ribogr%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%highcountcovtrs)%>%
			filter(between(readlen,25,35))%>%
			inner_join(offsets%>%select(offset,readlen))%>%
			mutate(start = start+offset)%>%
			{GRanges(.$transcript_id,IRanges(.$start,w=1),readlen=.$readlen)}%>%
			{seqlevels(.) = highcountcovtrs ;.}%>%
			{seqlengths(.) = trlens[highcountcovtrs] ;.}%>%
			coverage
	})
	names(psitecovlist)<-names(ribobams)
	saveRDS(psitecovlist,here('data/psitecovlist.rds'))
}else{
	psitecovlist<-readRDS(here('data/psitecovlist.rds'))
	stopifnot(names(psitecovlist)==names(ribobam))
}

if(!file.exists(here('data/psitecovnorm.rds'))){
	psitecovnorm = lapply(psitecovlist,function(psitecov){psitecov/(psitecov%>%sum%>%pmax(.,1))})
	saveRDS(psitecovnorm,here('data/psitecovnorm.rds'))
}else{
	psitecovnorm<-readRDS(here('data/psitecovnorm.rds'))
}

get_utr_exts <- function(cov,trspacecds,FPEXT,TPUTREXT){
	trlens = cov%>%runLength%>%sum
	seqnms = names(trlens)
	fputrlens = start(trspacecds[seqnms])-1
	fputrext = pmax(0,FPEXT - fputrlens)%>%setNames(seqnms)

	tputrlens = trlens - end(trspacecds[seqnms])
	tputrext = pmax(0,TPUTREXT - tputrlens)%>%setNames(seqnms)

	list(fputrext,tputrext)
}

ext_cov <- function(cov,fputrext,tputrext){
	gr = as(cov,"GRanges")
	seqnms = names(fputrext)
	stopifnot(identical(seqnms,names(tputrext)))
	stopifnot(identical(seqnms,names(cov)))
	seqlengths(gr)[seqnms]%<>%add(fputrext+tputrext)
	gr %<>% shift(fputrext)
	coverage(gr,weight='score')
}
{
#filter for only long enough ones, if we're doing windows
trspacecds = pmapToTranscripts(cdsgrl[highcountcovtrs],exonsgrl[highcountcovtrs])
trspacecds%<>%unlist

ltrspacecds = trspacecds
longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
ltrspacecds = ltrspacecds[longcdstrs]
#see how many we eliminate for length reasons
message(psitecovlist[[1]]%>%length)
message(length(ltrspacecds))
#spand the cds as necessary

cov = psitecovlist[[1]][longcdstrs]

c(fputrext,tputrext) %<-% (psitecovlist[[1]][longcdstrs]%>%get_utr_exts(ltrspacecds,FPEXT,TPUTREXT))
stopifnot(all(names(FPEXT)==longcdstrs))
stopifnot(all(names(tputrext)==longcdstrs))

if(!file.exists(here('data/epsitecovlist.rds'))){
	epsitecovlist  <- psitecovlist%>% 
		lapply(.%>%.[longcdstrs])%>%
		lapply(ext_cov,fputrext,tputrext)
	saveRDS(epsitecovlist,here('data/epsitecovlist.rds'))
}else{
	epsitecovlist<-readRDS(here('data/epsitecovlist.rds'))

}
stopifnot(all(names(epsitecovlist[[1]])==longcdstrs))
# stopifnot(seqnames(epsitecovlist[[1]])==longcdstrs)
#also lift cds to this new space
eltrspacecds <- ltrspacecds[longcdstrs]
seqlengths(eltrspacecds)[longcdstrs] %<>% add(fputrext+tputrext)
eltrspacecds %<>% shift(fputrext)
eltrspacecds%<>%resize(width(.)+FPEXT,'end')%>%resize(width(.)+TPUTREXT,'start')
#also the cds in this space

# highcovtrs = psitecovlist[[1]]%>%sum%>%sort%>%tail(500)%>%names%>%intersect(names(ltrspacecds))
# highcovtrs = names(ltrspacecds)

cov = epsitecovlist[[1]]
}
# file.remove(here('data/bindatamats.rds'))
if(!file.exists(here('data/bindatamats.rds'))){
	bindatamats = lapply(epsitecovlist,function(cov){
		midwind = eltrspacecds%>%
			resize(width(.)-(STOPWINDSIZE),'end')%>%
			resize(width(.)-(STARTWINDSIZE),'start')
		midmat = cov[midwind]%>%
			sum%>%#compress these variable length Rles down to 1 value per gene
			matrix#make a 1 column matrix
		# midmat = midmat%>%mean
		stwind = eltrspacecds%>%resize(STARTWINDSIZE,ignore.strand=TRUE)
		startmat = cov[stwind]%>%as.matrix
	#	%>%colMeans(na.rm=T)
		endwind = eltrspacecds%>%resize(STOPWINDSIZE,fix='end',ignore.strand=TRUE)
		endmat = cov[endwind]%>%as.matrix
		#%>%colMeans(na.rm=T)
		#sum codons
		# startmat%<>%matrix(nrow=3)%>%colSums
		# endmat%<>%matrix(nrow=3)%>%colSums
		message('.')
		#output as a data frame
		out = cbind(startmat,midmat,endmat)
		# c(startmat,midmat,endmat)%>%enframe('start','signal')%>%mutate(section=c(rep('AUG',150/3),'middle',rep('end',150/3)))
	})
	saveRDS(bindatamats,here('data/bindatamats.rds'))
}else{
	bindatamats<-readRDS(here('data/bindatamats.rds'))
}

names(bindatamats)
#functions for norm
rownorm <-function(x) x %>%sweep(.,MARGIN=1,F='/',STAT=rowSums(.)%>%{pmax(.,min(.[.>0])/10)})
#matrix <- matrix(1,6,24)
codmerge <- function(matrix,size=3) {
	matrix%>%t%>%rowsum(.,ceiling((1:nrow(.))/size),na.rm=T)%>%t
}

#name and check
{
matgenes = longcdstrs
bindatamats%<>%lapply(function(mat)set_rownames(mat,longcdstrs))
stopifnot(bindatamats[[1]]%>%dim%>%identical(c(length(longcdstrs),as.integer(TOTBINS))))
#why are the start and end mats so similiar?
startinds = 1:(STARTWINDSIZE)
endinds = (1+STARTWINDSIZE):TOTBINS
}


stop()



# fpcovlist <- fpcovlist[names(allbamtbls)]
# if(!file.exists(here('data/fpprofilelist.rds'))){
# 	fpprofilelist <-imap(fpcovlist[mainsamps],function(sampfpcov,sampname){
# 		trsums = sampfpcov%>%map(sum)%>%purrr::reduce(.,`+`)#sum over counts for that transcript
# 		sampfpcov%>%lapply(function(rlfpcov){
# 			rlfpcov = rlfpcov/(trsums)
# 			allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
# 			cods = names(allcodlistnz)%>%str_split('\\.')%>%map_chr(1)
# 			message('.')
# 			out = rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
# 			out
# 		})
# 	})
# 	saveRDS(fpprofilelist,here('data/fpprofilelist.rds'))
# }else{
# 	fpprofilelist<-readRDS(here('data/fpprofilelist.rds'))
# }


################################################################################
########Test PPG motifs
################################################################################
cdsseqs = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%extractTranscriptSeqs(x=fafileob)
cdsseqsaa <- cdsseqs%>%translate

peptides=c('PPE','PPG')

ppGgr = peptides%>%setNames(.,.)%>%map_df(.id='peptide',function(peptide){
		matchdf = vmatchPattern(peptide,cdsseqsaa)%>%
			as.data.frame%>%
		transmute(start,end,width,seqnames=names(cdsseqsaa)[group])
})%>%GRanges

gpeptidegr = ppGgr%>%{
	x=.
	x = resize(x,1,'start',ignore.strand=TRUE)
	end(x)  = ((start(x)-1)*3)+1
	start(x)  = end(x)
	# end(x)  = start(x)
	out = mapFromTranscripts(x,cdsgrl)
	out$transcript_id = names(cdsgrl)[out$transcriptsHits]
	out$peptide = x$peptide[out$xHits]
	out
	}
tr_peptidegr = gpeptidegr%>%mapToTranscripts(exonsgrl[highcountcovtrs])
tr_peptidegr$peptide = gpeptidegr$peptide[tr_peptidegr$xHits]

peps = unique(tr_peptidegr$peptide)%>%setNames(.,.)

matcolvalues = map_df(.id='peptide',peps,function(ipep){
	psitecovlist[T]%>%
	map_df(.id='sample',~.[tr_peptidegr%>%subset(peptide==ipep)%>%resize(9,'start')]%>%
		as.matrix%>%
		colMeans%>%
		enframe('position','mean_rdens')
	)
})
if('1' %in% matcolvalues$peptide) matcolvalues$peptide <- peps[as.numeric(matcolvalues$peptide)]

#now plot
for(ipep in peps){
	plotfile<- here(paste0('plots/',ipep,'_motif_Psite_Alignment','.pdf'))
	pdf(plotfile)
	pepchrs = str_split(ipep,'')%>%unlist
	p=matcolvalues%>%
		filter(peptide==ipep)%>%
		separate(sample,c('time','assay','rep'))%>%
		ggplot(.,aes(x=position,y=mean_rdens,color=time,group=paste0(time,rep)))+
		geom_line()+
		scale_color_manual(values=stagecols)+
		scale_x_continuous(paste0('Position'),breaks=1:9,labels=c('',pepchrs[1],'','',pepchrs[2],'','',pepchrs[3],''))+
		scale_y_continuous(paste0('Mean CDSnormed Ribosome Density'))+
		geom_vline(xintercept=4,linetype='dashed')+
		ggtitle(paste0(ipep,' Motif P-site alignment'))+
		theme_bw()
	print(p)
	dev.off()
	message(normalizePath(plotfile))
}

# }
################################################################################
########
################################################################################
	




# #
# #

# # startsigs = 
# fpcovlist = fpcovlist[fpcovlist%>%names%>%str_subset('ribo')]

# reduce=purrr::reduce
# startproflist = 
# 	imap(fpcov2use,function(sampfpcov,sampname){
# 		trsums = sampfpcov%>%map(sum)%>%reduce(`+`)#sum over counts for that transcript
# 		sampfpcov%>%imap(function(rlfpcov,rl){
# 			rl=as.numeric(rl)
# 			rlfpcov = rlfpcov/(trsums)
# 			# rlfpcov = sampfpcov[[rl]]
# 			rloffset = offsets%>%filter(length==rl)%>%.$offset
# 			startwinds = cdsstarts[highcountcovtrs]%>%
# 				enframe('seqnames','start')%>%mutate(end=start)%>%
# 				GRanges%>%
# 				{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
# 				{suppressWarnings({resize(.,3,'start')%>%resize(6,'end')%>%resize(45+3,'start')%>%shift(-rloffset)})}%>%
# 				.[!is_out_of_bounds(.)]
# 			#
# 			message('.')
# 			startwinds = startwinds%>%subset(seqnames%in%names(trsums)[trsums!=0])
# 			startwindsums = rlfpcov[startwinds]
# 		})
# 	})

# sampname = 'E13_ribo_1'
# sampfpcov = fpcovlist[[1]]

# endproflist = 
# 	imap(fpcovlist['E13_ribo_1'],function(sampfpcov,sampname){
# 		trsums = sampfpcov%>%map(sum)%>%reduce(`+`)#sum over counts for that transcript
# 		sampfpcov%>%imap(function(rlfpcov,rl){
# 			rl=as.numeric(rl)
# 			rlfpcov = rlfpcov/(trsums)
# 			# rlfpcov = sampfpcov[[rl]]
# 			rloffset = offsets%>%filter(length==rl)%>%.$offset
# 			startwinds = cdsends%>%
# 				enframe('seqnames','start')%>%mutate(end=start)%>%
# 				GRanges%>%
# 				{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
# 				{suppressWarnings({resize(.,3,'end')%>%resize(6,'start')%>%resize(45+3,'end')%>%shift(-rloffset)})}%>%
# 				.[!is_out_of_bounds(.)]
# 			#
# 			message('.')
# 			startwinds = startwinds%>%subset(seqnames%in%names(trsums)[trsums!=0])
# 			startwindsums = rlfpcov[startwinds]%>%as.matrix%>%colMeans
# 		})
# 	})



# startproflist[[1]][['28']][1:9]%>%txtplot
# endproflist[[1]][['28']][(48-2-6):48]%>%txtplot

# startproflist[[1]][['29']][1:9]%>%txtplot
# endproflist[[1]][['29']][(48-2-6):48]%>%txtplot



# startproflist[[1]]%>%map(.%>%matrix(byrow=TRUE,ncol=3)%>%{cbind(.,.-.)}%>%t%>%{txtplot(.,ylim=c(1,max(.)),height=20,width=100)})

# #
# startwinds = cdsstarts%>%
# 	enframe('seqnames','start')%>%mutate(end=start)%>%
# 	GRanges%>%
# 	{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
# 	{suppressWarnings({resize(.,3,'start')%>%resize(30,'start')%>%shift(-rloffset)})}%>%
# 	.[!is_out_of_bounds(.)]
# #

# startproflist[[1]][['28']]%>%matrix(byrow=TRUE,ncol=3)%>%{cbind(.,.-.)}%>%t%>%{txtplot(.,ylim=c(min(.)*2,max(.)),height=20,width=100)}


# startproflist[[1]][['27']][1:9]%>%txtplot
# startproflist[[1]][['28']][1:9]%>%txtplot
# startproflist[[1]][['29']][1:9]%>%txtplot
# startproflist[[1]][['30']][1:9]%>%txtplot


# startproflist[[1]]%>%reduce(`+`)%>%txtplot

# endproflist[[1]]%>%reduce(`+`)%>%txtplot


################################################################################
########Now get the fp dist around the codons
################################################################################
	
# allcodlist=codmatchwindowlist%>%GRangesList%>%unlist
# cods = names(allcodlist)%>%str_split('\\.')%>%map_chr(1)
# sampfpcov=fpcovlist[[1]]
# sampname = names(fpcovlist)[1]
# rlfpcov=sampfpcov$`29`
# if(!file.exists(here('data/fpprofilelist.rds'))){


# saveRDS(fpprofilelist,here('data/fpprofilelist.rds'))

# fpprofilelist<-readRDS(here('data/fpprofilelist.rds'))

# ################################################################################
# ########profile plots
# ################################################################################
# allpsitedfs_bind = startproflist%>%map_df(.id='sample',.%>%map_df(.id='readlen',enframe,'start','signal'))%>%mutate(section='AUG')

# allpsitedfs_bind %<>% group_by(sample,start,section)%>%summarise(signal=sum(signal))
# allpsitedfs_bind %<>% ungroup%>%mutate(stage=stageconv[str_extract(sample,'[^_]+')])

# allpsitedfs_bind%<>%group_by(sample,section)%>%mutate(signal=signal / median(signal))

# allpsitedfs_bind%<>%group_by(sample,section,start)%>%mutate(signal=sum(signal))
# # allpsitedfs_bind%<>%group_by(sample,section,start)%>%mutate(signal=sum(signal))
# allpsitedfs_bind%<>%group_by(sample,section)%>%mutate(signal=signal / median(signal))
# allpsitedfs_bind %<>% ungroup%>%mutate(stage=str_extract(sample,'[^_]+'))
# #
# allpsitedfs_bind_stgrp <- allpsitedfs_bind %>% group_by(stage,section,start)%>%summarise(signal=mean(signal)) 









# {
# library(rlang)
# plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=6,w=12)
# rwplot <- allpsitedfs_bind_stgrp%>%
# 	split(.,.$section)%>%map( .%>%
# 		# slice_by(sample,8)%>%
# 		# filter(position%%3 == 0)%>%
# 		{
# 			isfirst = .$section[1]=='AUG'
# 			qplot(data=.,color=stage,x=start,y=signal,geom='blank')+
# 			geom_line()+
# 			# scale_x_continuous(name='position',
# 			# 	# limits=if(.$section[1]=='AUG') c(-30,100) else c(-100,30) ,
# 				# limits=if(isfirst) c(-30,100) else c(-100,30) ,
# 			# 	minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
# 			# # scale_x_continuous(name='position',minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
# 			# # scale_x_continuous(name='position')+
# 			# # scale_color_discrete(guide=TRUE)+
# 			# # facet_grid( sample+readlen ~ section)+
# 			# facet_grid( ~ section,scale='free_x')+
# 			# scale_color_manual(values=stagecols)+
# 			# scale_y_continuous(name='Mean Psite Count / CDS Total',limits=c(0,20))+
# 			theme_bw()
# 		}
# 	)%>%
# 	ggarrange(plotlist=.,ncol=2,common.legend=T)
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }

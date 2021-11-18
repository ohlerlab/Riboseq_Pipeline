{
base::source(here::here('src/Rprofile.R'))
base::source(here::here('src/functions.R'))

if(!exists('fpcovlist')) base::source('src/load_pos_data.R')
STARTBUFF=60
ENDBUFF=60
LOWREADLIM=25
HIGHREADLIM=35
{
exonseq = exonsgrl[highcountcovtrs]%>%extractTranscriptSeqs(x=fafileob)
allcodons=getGeneticCode()
shift<-GenomicRanges::shift
i=1
innercds = trspacecds%>%
	subset(width>(3+STARTBUFF+ENDBUFF))%>%
	resize(width(.)-STARTBUFF,'end')%>%
	resize(width(.)-ENDBUFF,'start')
FLANKCODS=15
}




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

#
innercdspos <- innercds%>%{strand(.)<-'+';.}
cds_codons <- allcodlist%>%mapToTranscripts(innercdspos)
#
if(!file.exists(here('data/cdsfpcovlist.rds'))){
	cdsfpcovlist <- ribobams%>%mclapply(mc.cores=8,function(ribobam){
		ribogr <- GenomicAlignments::readGAlignments(ribobam)
		mcols(ribogr)$readlen <-  GenomicAlignments::qwidth(ribogr)
		ribogr%<>%as("GenomicRanges")
		cov = ribogr%>%resize(1)%>%mapToTranscripts(exonsgrl[names(innercdspos)])
		cov$readlen = mcols(ribogr)$readlen[cov$xHits]
		cov%<>%subset(between(readlen,LOWREADLIM,HIGHREADLIM))
		cov%<>%resize(1,'start')
		cdscov <- cov%>%mapToTranscripts(innercdspos)
		cdscov$readlen <- cov$readlen[cdscov$xHits]
		cdscov$xHits <- NULL
		cdscov$transcriptsHits<-NULL
		split(cdscov,cdscov$readlen)%>%
			lapply(coverage)
	})
	cdsfpcovlist%<>%setNames(samples)
	saveRDS(cdsfpcovlist,here('data/cdsfpcovlist.rds'))
}else{
	cdsfpcovlist<-readRDS(here('data/cdsfpcovlist.rds'))
}

# cds_codons_nz = cds_codons%>%subset(
mcols(cds_codons)=NULL
cds_codons$codon = names(cds_codons)%>%str_split('\\.')%>%map_chr(1)
sampinds = seq_along(samples)


# if(!file.exists(here('data/cn_norm.rds'))){
# 	cn_norm <- 	lapply(dtselgenelist['all'],function(seltrs){
# 			# mcmapply(mc.cores=4,cdsfpcovlist[mainsamps[1:10]],mainsamps[1:10],F=function(cdssampfpcov,sampname){
# 			mcmapply(mc.cores=4,SIMPLIFY=F,cdsfpcovlist[mainsamps[sampinds]],mainsamps[sampinds],F=function(cdssampfpcov,sampname){
# 				cdssampfpcov%>%lapply(function(cdsrlfpcov){
# 					rustcdsrlfpcov <- cdsrlfpcov
# 					# rustcdsrlfpcov <- rustcdsrlfpcov>mean(rustcdsrlfpcov)
# 					nz_trs <- any(rustcdsrlfpcov)%>%names(.)[.]
# 					#
# 					cds_codons_nz = cds_codons%>%subset(seqnames%in%nz_trs)
# 					codtrs = seqnames(cds_codons_nz)%>%as.character
# 					cat('.')
# 					#calculate ro vals
# 					ro_cl = rustcdsrlfpcov[cds_codons_nz]%>%
# 						split(cds_codons_nz$codon)%>%
# 						lapply(as.matrix)%>%
# 						map(colMeans)
# 					#also get evals
# 					tr_rust_evals <- rustcdsrlfpcov%>%mean
# 					re_c <- tr_rust_evals[codtrs]%>%
# 						split(cds_codons_nz$codon)%>%
# 						map_dbl(mean)
# 					ro_cl <- ro_cl%>%map_df(.id='codon', enframe, 'position', 'ro_cl')
# 					re_c <- enframe(re_c, 'codon', 're_c')
# 					ro_cl%>%left_join(re_c, by='codon')
# 				})
# 			})%>%setNames(mainsamps[sampinds])
# 		})
# 	saveRDS(cn_norm,here('data/cn_norm.rds'))
# }else{
# 	cn_norm<-readRDS(here('data/cn_norm.rds'))
# }

dtselgenelist=list(all=highcountcovtrs)
if(!file.exists(here('data/cdsfpcovlist.rds'))){
	rust_roel <- lapply(dtselgenelist['all'],function(seltrs){
		mcmapply(mc.cores=4,SIMPLIFY=F,cdsfpcovlist[samples],samples,F=function(cdssampfpcov,sampname){
			cdssampfpcov%>%lapply(function(cdsrlfpcov){
					rustcdsrlfpcov <- cdsrlfpcov
					rustcdsrlfpcov <- rustcdsrlfpcov>mean(rustcdsrlfpcov)
					nz_trs <- any(rustcdsrlfpcov)%>%names(.)[.]
					#
					cds_codons_nz = cds_codons%>%subset(seqnames%in%nz_trs)
					codtrs = seqnames(cds_codons_nz)%>%as.character
					cat('.')
					#calculate ro vals
					ro_cl = rustcdsrlfpcov[cds_codons_nz]%>%
						split(cds_codons_nz$codon)%>%
						lapply(as.matrix)%>%
						map(colMeans)
					#also get evals
					tr_rust_evals <- rustcdsrlfpcov%>%mean
					re_c <- tr_rust_evals[codtrs]%>%
						split(cds_codons_nz$codon)%>%
						map_dbl(mean)
					ro_cl <- ro_cl%>%map_df(.id='codon', enframe, 'position', 'ro_cl')
					re_c <- enframe(re_c, 'codon', 're_c')
					ro_cl%>%left_join(re_c, by='codon')
			})
		})%>%setNames(samples)
	})
	saveRDS(rust_roel,here('data/rust_roel.rds'))
}else{
	rust_roel<-readRDS(here('data/rust_roel.rds'))
}

#https://www.nature.com/articles/ncomms12915#Sec10 see equation 3
# of RUST paper
frustprofilelist <- rust_roel%>%
# frustprofilelist <- cn_norm%>%
	map_df(.id='set',.%>%map_df(.id='sample',.%>%bind_rows(.id='readlen')))
#
frustprofilelist %<>% mutate(position = position - 1 - (FLANKCODS*3))
	frustprofilelist %<>% filter(!codon %in% c('TAG','TAA','TGA'))
frustprofilelist%<>%mutate(count = ro_cl/re_c)
#
frustprofilelist$sample%>%table
frustprofilelist$readlen%>%n_distinct
# frustprofilelist$sample %<>% {mainsamps[as.numeric(.)]}
# frustprofilelist$readlen %<>% as.numeric %>% names(fpcovlist[[1]])[.]
frustprofilelist$nreadlen <-frustprofilelist$readlen%>%as.numeric
frustprofilelist$readlen%<>%str_replace('^(\\d)','rl\\1')
#
kl_df<-frustprofilelist%>%
	filter(set=='all')%>%
	# group_by(readlen)%>%group_slice(4)%>%
	# group_by(position)%>%
	# group_slice(1)%>%
	# filter(position%%3 ==0)%>%
	# .$position
	# filter(position< -6)%>%
	# filter(position> -(nreadlen-6))%>%
	group_by(set,sample,nreadlen,position)%>%
	# summarise(KL=sum(ro_cl * log2(ro_cl/re_c)))
	mutate(ro_cl = ro_cl/sum(ro_cl), re_c = re_c/sum(re_c))%>%
	summarise(KL=sum(ro_cl * log2(ro_cl/re_c)))
# #
# kl_df%>%
# 	mutate(phase=position%%3)%>%
# 	group_by(nreadlen,phase)%>%group_slice(1)%>%{txtplot(.$position,.$KL)}
# #get offsets by picking top two KLs per phase, choosing rightmost as psite
# kl_offsets <- kl_df%>%
# 	mutate(phase=position%%3)%>%
# 	group_by(set,sample,nreadlen,phase)%>%
# 	filter(position< -6)%>%
# 	filter(position> -(nreadlen-6))%>%
# 	mutate(rank=rank(-KL))%>%
# 	filter(rank%in%c(1:2))%>%
# 	arrange(sample,nreadlen,phase,position)
# kl_offsets <- kl_df%>%
# 	mutate(phase=position%%3)%>%
# 	group_by(set,sample,nreadlen,phase)%>%
# 	arrange(position)%>%
# 	# mutate(codon_KL = KL+lag(KL)+lag(lag(KL)))%>%
# 	mutate(codon_KL = KL)%>%
# 	filter(position< -6)%>%
# 	filter(position> -(nreadlen-6))%>%
# 	group_by(sample,nreadlen)%>%
# 	slice(which.max(codon_KL))%>%
# 	arrange(sample,nreadlen,phase,position)
# #
# kl_offsets <- kl_offsets%>%
# 	group_by(sample,nreadlen,phase)%>%
# 	# filter(sample=='E13_ribo_1',nreadlen==25)%>%
# 	slice(which.max(position))%>%
# 	summarise(offset=-(position+3))%>%
# 	mutate(readlen=paste0('rl',nreadlen))
#
kl_offsets <- kl_df%>%group_by(nreadlen,position)%>%
	summarise(sumKL = sum(KL))%>%
	filter(position< -6)%>%
	filter(position> -(nreadlen-6))%>%
	slice(which.max(sumKL))%>%
	mutate(p_offset = position+3)
#
most_freq <- function(x) x%>% table%>%sort%>%names%>%as.numeric%>%tail(1)
#no equalize all timepoints
# kl_offsets <- kl_offsets%>%
	# separate(sample,c('time','assay','rep'),remove=F)%>%
	# filter(assay=='ribo')%>%
	# group_by(nreadlen)%>%
	# group_slice(5)%>%
	# mutate(offset = offset%>%most_freq)
#
clean_fr_sampnames<-function(x) x%>%str_replace('.*_(Poly|80S)(.*)_()','\\2_\\1ribo_\\3')
kl_offsets2plot <- kl_offsets%>%
	# separate(sample,c('time','assay','rep'),remove=F)%>%
	# filter(assay=='ribo')
	mutate(readlen = paste0('rl',nreadlen))
#

selreadlens=27:33
pdf<-grDevices::pdf
dir.create('plots')
plotfile='plots/rust_fppos_vs_codon_variance.pdf'
# pdf(plotfile,w=12,h=3*n_distinct(kl_df$nreadlen))
pdf(plotfile,w=12,h=1+3*length(selreadlens))
#plotting variance amongst codons at each point.
# sh_codprof%>%
kl_df%>%
	filter(position< -3)%>%
	filter(position> -(nreadlen-6))%>%
	separate(sample,c('fraction','genotype','rep'),remove=F)%>%
	filter(nreadlen%in%selreadlens)%>%
	{
		qplot(data=.,x=position,y=KL)+
		theme_bw()+
		facet_grid(nreadlen~genotype+fraction)+
		scale_y_continuous('RUST KL divergence')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=kl_offsets2plot%>%filter(nreadlen%in%selreadlens),aes(xintercept= p_offset-3),color=I('green'),linetype=2)+
		geom_vline(data=kl_offsets2plot%>%filter(nreadlen%in%selreadlens),aes(xintercept= p_offset),color=I('blue'),linetype=2)+
		ggtitle("RUST KL divergence vs position")
	}%>%print
dev.off()
message(normalizePath(plotfile))
#
kl_offsets2plot%>%select(p_offset,length=nreadlen)%>%write_tsv('ext_data/offsets_rustvar.tsv')





#posseldf
posseldf = bind_rows(
	kl_offsets2plot %>% filter(nreadlen%in%28:31)%>%mutate(position=position+3,site='p_site')%>%select(nreadlen,position,site),
	kl_offsets2plot %>% filter(nreadlen%in%28:31)%>%mutate(position=position,site='a_site')%>%select(nreadlen,position,site),
	kl_offsets2plot %>% filter(nreadlen%in%28:31)%>%mutate(position=position-3,site='a_p3_site')%>%select(nreadlen,position,site)
)

allcodondt = frustprofilelist%>%
	filter(set=='all')%>%
	group_by(set,sample,nreadlen,position)%>%
	filter(nreadlen %in% 28:31)%>%
	mutate(ro_cl = ro_cl/sum(ro_cl), re_c = re_c/sum(re_c))%>%
	inner_join(posseldf)%>%
	ungroup%>%
	select(-set,-readlen,-re_c,-count,-position)%>%
	dplyr::rename('RUST_score':=ro_cl)%>%
	pivot_wider(names_from='site',values_from='RUST_score')

allcodondt %>%write_tsv('tables/allcodondt.tsv')
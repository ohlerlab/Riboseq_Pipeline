source('src/Rprofile.R')
#okay so, lets get the overlapping portion of each cluster 
ovclustdf = 'ext_data/processed_data_parclip_overlapping_clusters.txt'%>%
	read_tsv
#get overlapping region
ovclustdf=ovclustdf%>%mutate_at(vars(matches('start|end')),list(as.numeric))%>%
	mutate(ovstart = pmax(start1,start2),ovend=pmin(stop1,stop2))
ovclustgr = GRanges(ovclustdf$chr1,IRanges(ovclustdf$ovstart,ovclustdf$ovend),ovclustdf$strand)
ovclustgr = keepStandardChromosomes(ovclustgr,pruning='coarse')
ovclustgrtr = ovclustgr%>%mapToTranscripts(exonsgrl)
library(liftOver)
#system('wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz')
#system('guzip hg19ToHg38.over.chain.gz')
chainpath = here('hg19ToHg38.over.chain')
ch = import.chain(chainpath)
# seqlevelsStyle(cur) = "UCSC"  # necessary
lovclustgr = liftOver(ovclustgr, ch)
lovclustgr%<>%unlist
#indeed these look right
lovclustgr$name <- paste0('ovcluster_',seq_along(lovclustgr))
lovclustgr%>%rtracklayer::export('test38.bed')

covtrs = names(fpcovlist[[1]][[1]])

tx2genemap<- gtf_gr%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id)

membgids = read_tsv('ext_data/dermott_membrane_localized.txt')%>%filter(localization_cat=='membrane')%>%.$gene_id
membtrs <- gtf_gr%>%subset(trimids(gene_id)%in%trimids(membgids))%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id)%>%
	filter(transcript_id%in%covtrs)%>%.$transcript_id
exonsgrl = keepStandardChromosomes(exonsgrl,pruning='coarse')
cdsgrl=gtf_gr%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st
metasitewindows = lovclustgr%>%resize(1,'center')%>%
	mapToTranscripts(cdsgrl[membtrs])%>%
	{strand(.)<-'+';.}%>%
	resize(101,'center')%>%
	{.=.[!is_out_of_bounds(.)]}
#only count once per gene
metasitewindows$gene_id = tx2genemap$gene_id[as.vector(match(seqnames(metasitewindows),tx2genemap$transcript_id))]
metasitewindows <- metasitewindows[match(unique(metasitewindows$xHits),metasitewindows$xHits)]
width1grs <- function(gr){
	stopifnot(Negate(is.unsorted)(gr))
	isw1 <- width(gr)==1
	broad <- gr[!isw1]
	#vector of integers - 1,2,3 for range 1-3
	narrowstarts <- unlist(as(broad@ranges,'IntegerList'))
	narrow <- {GRanges(
			rep(seqnames(broad),width(broad)),
			IRanges(narrowstarts,w=1)
		)}
	mcols(narrow) <- mcols(broad)[rep(seq_along(broad),width(broad)),,drop=F]
	sort(c(gr[isw1],narrow))
}
samples = names(fpcovlist)
i_sample=samples[1]
meta_site_profiles <- mclapply(mc.cores=8,samples,function(i_sample){
	fpcomb = fpcovlist[[i_sample]][c('28','29','30','31')]%>%Reduce(f='+')
	fpcomb = fpcomb%>%as("GRanges")%>%width1grs%>%subset(score!=0)%>%
		shift(12)%>%{.=.[!is_out_of_bounds(.)]}%>%coverage(weight='score')
	fpsums = fpcomb%>%sum
	fpdens = fpsums/(sum(width(cdsgrl))[names(fpsums)])
	library(txtplot)
	fppossums = fpcomb[metasitewindows]%>%as.matrix%>%colSums
	fpposexpts = sum(fpdens[as.vector(seqnames(metasitewindows))])
	# txtplot(fppossums/fpposexpts)
	cat('.')
	fppossums/fpposexpts
})
meta_site_profiles%<>%setNames(samples)
# mutate(pos = name )
plotfile <- here(paste0('plots/','meta_parclip','.pdf'))
pdf(plotfile)
meta_site_profiles%>%map_df(.id='sample',enframe)%>%
		mutate(name=name-51)%>%
		separate(sample,c('fraction','genotype','rep'),remove=FALSE)%>%
		ggplot(aes(x=name,y=value,color=genotype,group=sample))+
		geom_line()+
		scale_x_continuous('P site position relative to binding ov midpoint')+
		scale_y_continuous('P site Enrichment')+
		facet_grid(fraction~.)+
		theme_bw()
dev.off()
message(normalizePath(plotfile))


# mutate(pos = name )
plotfile <- here(paste0('plots/','meta_parclip_zoom','.pdf'))
pdf(plotfile)
meta_site_profiles%>%map_df(.id='sample',enframe)%>%
		mutate(name=name-51)%>%
		filter(name%>%between(-20,20))%>%
		filter(sample!='Membrane_KO_R4')%>%
		separate(sample,c('fraction','genotype','rep'),remove=FALSE)%>%
		ggplot(aes(x=name,y=value,color=genotype,group=sample))+
		geom_line()+
		scale_x_continuous('P site position relative to binding ov midpoint')+
		scale_y_continuous('P site Enrichment')+
		facet_grid(fraction~.)+
		# geom_text(data=tibble(labl=),
			# aes(label=labl,x=Inf,y=Inf),hjust=1,vjust=1)+
		theme_bw()
dev.off()
message(normalizePath(plotfile))

#so Membrane_KO_R1/4 have this huge peak.


#now plot
plotfile<- here(paste0('plots/','pca','.pdf'))
pdf(plotfile,w=14,h=14)
pca = DESeq2::plotPCA(rld,intgroup="group")
pca = pca +
  expand_limits(x = pca$data$PC1 %>%
                  range %>%
                  multiply_by(1.5)) +
  expand_limits(y = pca$data$PC2 %>%
                  range %>%
                  multiply_by(1.3))
print((pca+geom_text(aes(label=name),size=I(2),alpha=I(0.5)))%>%{.$layers=.$layers[-1];.+geom_point(size=I(0))}+ggtitle('PCA : labelled by group')+scale_color_discrete(guide=F))
dev.off()
message(normalizePath(plotfile))
print((pca+geom_text(aes(label=name),size=I(2),alpha=I(0.5)))%>%{.$layers=.$layers[-1];.+geom_point(size=I(0))}+ggtitle('PCA : labelled by group'))





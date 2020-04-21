library(GenomicFeatures)
source('src/Rprofile.R')
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}

anno <- projmemoise(function(...){rtracklayer::import(...)})(here('pipeline/gencode.vM12.annotation.gtf'))
dds <- readRDS('data/dds.rds')
#We want a function that 
cds <- anno%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sort_grl_st
exons <- anno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st 
# sorted_is_dj <- function(grl){
#   starts <- start(grl)%>%{.[ IntegerList(as.list(rep(-1,length((.)))) ) ]}
#   ends = start(grl)%>%{.[ IntegerList(as.list(-elementNROWS(.)) ) ]}
#   all(all(starts < ends))
# }
# sorted_is_dj(cds)

cdscovs <- mclapply(bams,function(bam)bamsignals::bamCount(cds%>%unlist,bampath=bam)%>%split(cds@partitioning)%>%sum)
cdscovdf <- cdscovs%>%setNames(bams%>%dirname%>%basename)%>%map_df(.id='sample',enframe,'protein_id','count')

p2giddf<-mcols(cds@unlistData)[,c('protein_id','gene_id')]%>%as.data.frame%>%distinct
bestcdscovs<-cdscovdf%>%left_join(p2giddf)%>%group_by(gene_id,protein_id)%>%mutate(medcount=median(count))%>%group_by(gene_id)%>%arrange(desc(medcount))%>%
	filter(protein_id==protein_id[1])%>%select(-medcount)
bestcdscovs %<>% filter(max(count)>32)




exprgenes <- counts(dds)%>%rowSums%>%add(1)%>%log10 %>%`>`(2)%>%names

cds <- cds[split(cds@unlistData@elementMetadata$gene_id %in% exprgenes,cds@partitioning)]

maxlen_pids <- sum(width(cds))%>%enframe('protein_id','len')%>%
    left_join(data.frame(mcols(cds@unlistData)[,c('protein_id','gene_id')])%>%distinct,by='protein_id')%>%
    group_by(gene_id)%>%slice(which.max(len))%>%
    .$protein_id

mcds = cds[maxlen_pids]

expexons <- exons[fmcols(mcds,transcript_id)]%>%trim_grl(-100,end='fp')%>%trim_grl(-100,end='tp')

#now map our start codons to the expanded transcript space
startwinds = pmapToTranscripts(mcds%>%resize_grl(1)%>%unlist,expexons)%>%
    resize(100,'start',ignore.strand=TRUE)%>%#expand them to desired window size
    resize(200,'end',ignore.strand=TRUE)%>%
    pmapFromTranscripts(expexons)#and map back to the genome

startwinds <- startwinds[fmcols_List(startwinds,hit)]

bams <- Sys.glob(here('pipeline/star/data/*/*ribo*.bam'))



sampwindcovs <- mclapply(mc.cores=20,bams,startwinds,F=projmemoise(function(bam,startwinds){
	st_psites <- get_genomic_psites(bam,startwinds%>%unlist,cutoffs)
	windcov <- st_psites%>%mapToTranscripts(startwinds)%>%coverage
	windcov
}))

mcdscores <- mcds%>%head(1e3)%>%.[sum(width(.))>24]%>%trim_grl((5*3),end='fp')%>%trim_grl((3*3),'tp')%>%{countOverlaps(.,get_genomic_psites(bams[[1]],unlist(.),cutoffs))}


longcds <- cds[sum(width(cds))>24]
psites <- longcds%>%unlist%>%GenomicRanges::reduce(.)%>%{get_genomic_psites(bams[[1]],.,cutoffs)}
startcounts <- cds%>%head(100)%>%resize_grl(15)%>%countOverlaps(psites)
endcounts <- mcds%>%head(100)%>%resize_grl(9,'end')%>%countOverlaps(psites)


fread('pipeline/feature_counts/all_feature_counts')%>%.[,-1]%>%{DESeq2::estimateSizeFactorsForMatrix(.)}

fread(str_interp('samtools view ${bams[1]}  | head -n 10000000'))

bamsignals::bamCoverage(cds,bams[1])


bestcds <- bestcdscovs$protein_id%>%unique

psites <- cds[bestcds%>%head(2)]%>%{get_genomic_psites(bams[[1]],unlist(.),cutoffs)}




coverage(psites)[cds[[1]]]


# #now plot
# plotfile<- here(paste0('plots/',,'.pdf'))
# pdf(plotfile)
# %>%
# 	ggplot(.,aes())+
# 	scale_color_discrete(name='colorname',colorvals)+
# 	scale_x_continuous(paste0('xname'))+
# 	scale_y_continuous(paste0('yname'))+
# 	ggtitle(paste0('title'))+
# 	theme_bw()
# dev.off()
# normalizePath(plotfile)

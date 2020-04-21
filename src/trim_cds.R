#load the gtf
library(GenomicFeatures)
source('src/Rprofile.R')
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
DWPT         =       function(signal){
  library(wmtsa);
  TR_lth          =       length(signal);
  if(TR_lth       <       64){
          signal  =       c(signal,rep(0,64-length(signal)));
  }
  W1              =       wavMODWPT(signal, wavelet="s4",n.levels=6);
  W2              =       wavShift(W1);
  bands           =       32:51;# use the 0.2~0.5 Hz components only.
  mx              =       matrix(0,nrow=length(bands),ncol=length(signal));
  for(i in 1:length(bands)){
          tmp     =       paste("w6.",bands[i],sep="");
          mx[i,]  =       W2$data[[tmp]];
  }
  ID_signal       =       which(signal>0);        #the positions with signal;
  mx[,-ID_signal] =       0;                      #remove noise;
  minus3nt        =       mx[-11,];
  only3nt         =       mx[11,];
  BKgrnd          =       apply(minus3nt,2,max);
  ID1             =       which(only3nt > BKgrnd);#the positions with 3nt energy higher than other frequency;
  higher3nt       =       signal;
  higher3nt[-ID1] =       0;
  out             =       higher3nt[1:TR_lth];
}
#load copy of the gtf in the pipeline folder
gtf_file = here::here('pipeline/',basename(yaml::yaml.load_file('src/config.yaml')$GTF_orig)) %T>%{stopifnot(file.exists(.))}
anno <- projmemoise(function(...){rtracklayer::import(...)})(gtf_file)
#We want a function that 
allcds <- anno%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sort_grl_st
exons <- anno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st 

STARTTRIMCODS <- yaml::yaml.load_file('src/config.yaml')$STARTCODTRIM%T>%{stopifnot(!is.null(.))}
ENDTRIMCODS <- yaml::yaml.load_file('src/config.yaml')$STOPCODTRIM %T>%{stopifnot(!is.null(.))}
TRIMLEN = 3*(STARTTRIMCODS+ENDTRIMCODS)
longlen = (TRIMLEN+30)

shortcds <- allcds[sum(width(allcds))<=longlen]

cds <- allcds[sum(width(allcds)) > longlen]

cds_trimmed <- cds%>%trim_grl((STARTTRIMCODS*3),end='fp')%>%trim_grl((ENDTRIMCODS*3),'tp')
cds_starts <- cds%>%resize_grl((STARTTRIMCODS*3),'start')
cds_ends <- cds%>%resize_grl((ENDTRIMCODS*3),'end')

cds_trimmed%>%unlist%>%rtracklayer::export(gtf_file%>%str_replace('\\.gtf$','.trimmed.cds.gtf'))
shortcds%>%unlist%>%rtracklayer::export(gtf_file%>%str_replace('\\.gtf$','.short.cds.gtf'))

unique(mcols(anno)[,c('gene_id','transcript_id')])%>%write.table(row.names=F,here('pipeline','gid2trid.txt'))
unique(mcols(anno)[,c('transcript_id','protein_id')])%>%write.table(row.names=F,here('pipeline','trid2prid.txt'))
unique(mcols(anno)[,c('gene_id','protein_id')])%>%write.table(row.names=F,here('pipeline','gid2prid.txt'))
unique(mcols(anno)[,c('gene_name','gene_id')])%>%write.table(row.names=F,here('pipeline','gnm2gid.txt'))
here('pipeline','gnm2gid.txt')%>%fread

gid2piddf<-fread(here('pipeline','gid2prid.txt'))

cutoffs <- Sys.glob('pipeline/riboseqc/data/*/_P_sites_calcs')%>%map_df(fread)%>%
  group_by(read_length,comp)%>%
  filter(read_length-cutoff > 3)%>%
  summarise(cutoffs = list(as.data.frame(table(cutoff))))%>%
  unnest(cutoffs)%>%
  mutate(cutoff = as.character(cutoff)%>%as.numeric)%>%
  arrange(desc(Freq))%>%
  slice(1)%>%
  rename(offset:=cutoff,length:=read_length)
bams <- Sys.glob(here('pipeline/star/data/*/*ribo*.bam'))%>%setNames(.,basename(dirname(.)))
bamseqnames <- seqinfo(Rsamtools::BamFile(bams[1]))@seqnames
sample = names(bams)[1]
chr = 'chrM'

mclapply(mc.cores=12,names(bams),function(sample){

  cdstrimcounts = rep(NA,length(cds_trimmed))
  cdsstartcounts = rep(NA,length(cds_trimmed))
  cdsendcounts = rep(NA,length(cds_trimmed))
  cdsshortcounts = rep(NA,length(shortcds))

  for(chr in unique(intersect(unique(seqnames(unlist(cds))),bamseqnames))){
    ischr <- as.vector(seqnames(cds_trimmed@unlistData)[start(cds_trimmed@partitioning)]==chr)
    redcds <- cds[ischr]%>%unlist%>%GenomicRanges::reduce(.)

    #get psitesj
    chrpsites <- get_genomic_psites(bams[sample],redcds,cutoffs)

    #quantify
    cdstrimcounts[ischr] = countOverlaps(cds_trimmed[ischr],chrpsites)
    cdsstartcounts[ischr] = countOverlaps(cds_starts[ischr],chrpsites)
    cdsendcounts[ischr] = countOverlaps(cds_ends[ischr],chrpsites)

    #for short cds
    shortischr <- as.vector(seqnames(shortcds@unlistData)[start(shortcds@partitioning)]==chr)
    cdsshortcounts[shortischr] <-  countOverlaps(shortcds[shortischr],chrpsites)

    message(chr)
  }

  outfolder=here('pipeline','cdscounts',sample)%T>%dir.create(.,showWarn=F,rec=TRUE)
  outfile = file.path(outfolder,'trimcounts.tsv')
  message(str_interp('writing to ${outfile}'))

  cdscountdf <- tibble(
    protein_id = names(cds_trimmed),
    trimcount = cdstrimcounts ,
    startcount = cdsstartcounts,
    endcount = cdsendcounts,
    totalcount = cdstrimcounts+cdsstartcounts+cdsendcounts
  ) %>% left_join(gid2piddf,by='protein_id') %T>% write_tsv(outfile)

  outfile = file.path(outfolder,'shortcounts.tsv')
  message(str_interp('writing to ${outfile}'))

  shortcdscountdf <- tibble(
    protein_id = names(shortcds),
    shortcount = cdsshortcounts 
  ) %>% left_join(gid2piddf,by='protein_id') %T>% write_tsv(outfile)

})

#extend exons so we can always make our windows
bakcds <- cds
cds <- bakcds

spl_mapFromTranscripts<-function(trspacegr,exons_grl){
  exons_tr<-exons_grl%>%unlist%>%mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  ov <- findOverlaps(trspacegr,exons_tr)
  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
    trspacegr_spl,
    exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}
bestcds <- cdscountdf%>%filter(totalcount>32)%>%group_by(gene_id)%>%slice(which.max(trimcount))%>%.$protein_id
bestcds <- bestcds%>%setdiff(cds%>%unlist%>%.[seqnames(.)%in%'chrM']%>%.$protein_id)
bestcds <- cds[bestcds]
#cds<-cdsbak
#cdsbak <- cds
windback <- 50
windforw <- 150
exonsexp <- exons[fmcols(bestcds,transcript_id)]%>%resize_grl(sum(width(.))+windback,'end',check=FALSE)%>%resize_grl(sum(width(.))+windforw,'start',check=FALSE)%>%.[!any(is_out_of_bounds(.))]


#now 
metaplotwinds <- pmapToTranscripts (bestcds%>%resize_grl(1),exonsexp[fmcols(bestcds,transcript_id)])%>%
  unlist%>%
  resize(windback+1,fix='end',ignore.strand=TRUE)%>%
  resize(windback+windforw,'start',ignore.strand=TRUE)%>%
  {spl_mapFromTranscripts(.,exons_grl = exonsexp[seqnames(.)])}
metaplotwinds%<>%split(.,names(.))

metaplotmats <- mclapply(bams,function(bam){
  psites <- get_genomic_psites(bam,unlist(metaplotwinds),cutoffs)
  psites%>%mapToTranscripts(metaplotwinds%>%subsetByOverlaps(chrpsites))%>%coverage%>%as.matrix
})

meatplotdf<-metaplotmats%>%setNames(names(bams))%>%map_df(.id='sample',as.data.frame)%>%gather(pos,count,-sample)
meatplotdf$pos%<>%as_factor%>%as.numeric

condition = names(bams)%>%setNames(str_extract(.,'mock|si'),.)
#now plot
plotfile<- here(paste0('plots/','metaplot','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
meatplotdf%>%
  group_by(sample,pos)%>%summarise(count=sum(count))%>%
  mutate(condition = condition[sample])%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 51+c(-30,0,30,60,90),labels=c('-30bp','AUG','30bp','60bp','90bp'))+
  scale_y_continuous(paste0('Psite Count'))+
  ggtitle(paste0('Psite distribution Metaplot'))+
  facet_grid(condition ~ . )+
  theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
plotfile<- here(paste0('plots/','metaplot_codbin','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
meatplotdf%>%
  mutate(pos = pos - (pos - 51)%%3 )%>%
  group_by(sample,pos)%>%summarise(count=sum(count))%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  mutate(condition = condition[sample])%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 51+c(-30,0,30,60,90),labels=c('-30bp','AUG','30bp','60bp','90bp'))+
  scale_y_continuous(paste0('Psite Count'))+
  ggtitle(paste0('Psite distribution Metaplot'))+
  facet_grid(condition ~ . )+
  theme_bw()
dev.off()
normalizePath(plotfile)

dspdf <- fread('ext_data/no_DSP.csv')
dspdf%<>%filter(CDS_reads_IP>32 &(CDS_reads_total>32))

bestcds@elementMetadata$isindsp <- fmcols(bestcds,gene_name) %in% dspdf$Gene_name
bestcds@elementMetadata$is_target <- fmcols(bestcds,gene_name) %in% (dspdf%>%filter(RPKM_ratio > 1.5)%>%.$Gene_name)



#now plot
plotfile<- here(paste0('plots/','metaplot_cod_targetsep','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
meatplotdf%>%group_by(sample,pos)%>%mutate(cds=1:n())%>%
  mutate(isindsp = bestcds@elementMetadata$isindsp[cds])%>%
  mutate(is_target = bestcds@elementMetadata$is_target[cds])%>%
  filter(isindsp)%>%
  ungroup%>%
  mutate(pos = pos - (pos - 51)%%3 )%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  group_by(sample,is_target,pos)%>%summarise(count=sum(count))%>%
  group_by(sample,is_target)%>%mutate(count = count/sum(count))%>%
  mutate(condition = condition[sample])%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 51+c(-30,0,30,60,90),labels=c('-30bp','AUG','30bp','60bp','90bp'))+
  scale_y_continuous(paste0('Psite Count'))+
  ggtitle(paste0('Psite distribution Metaplot'))+
  facet_grid(is_target ~ . )+
  theme_bw()
dev.off()
normalizePath(plotfile)





metaplotmats

fafile=file.path(here('pipeline',yaml::yaml.load_file('src/config.yaml')%>%.$REF_orig%>%basename%>%str_replace('.gz$','')))%T>%{stopifnot(file.exists(.))}

testatgs <- cds[bestcds]%>%head%>%resize_grl(3,'start')%>%.[elementNROWS(.)==1]%>%head%>%unlist

getSeq(FaFile(fafile),testatgs)

testatgs%>% resize(1)%>%{pmapToTranscripts (.,exonsexp[testatgs$transcript_id])}%>%
  # resize(3,fix='end',ignore.strand=TRUE)%>%
  resize(3,'start',ignore.strand=TRUE)%>%
  {spl_mapFromTranscripts(.,exons_grl = exonsexp[seqnames(.)])}%>%{getSeq(FaFile(fafile),.)}

metaplotwinds%>%.[elementNROWS(.)==1]%>%resize(150,'end')%>%resize(3,'start')%>%head%>%unlist%>%getSeq(FaFile(fafile),.)
metaplotwinds




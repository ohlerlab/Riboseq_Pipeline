
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


#cds<-cdsbak
#cdsbak <- cds
windback <- 50
windforw <- 150
stopexonsexp <- exons[fmcols(bestcds,transcript_id)]%>%
  resize_grl(sum(width(.))+windback,'start',check=FALSE)%>%
  resize_grl(sum(width(.))+windforw,'end',check=FALSE)%>%
  .[!any(is_out_of_bounds(.))]

#now 
stopmetaplotwinds <- pmapToTranscripts (bestcds%>%head%>%resize_grl(1,'end'),stopexonsexp[fmcols(bestcds,transcript_id)])%>%
  unlist%>%
  shift(1)%>%#move tot he actual stop codon
  resize(windforw +1,fix='end',ignore.strand=TRUE)%>%
  resize(windforw+ windback,'start',ignore.strand=TRUE)%>%
  {spl_mapFromTranscripts(.,exons_grl = stopexonsexp[seqnames(.)])}
stopmetaplotwinds%<>%split(.,names(.))
stopmetaplotmats <- mclapply(bams,function(bam){
  psites <- get_genomic_psites(bam,unlist(stopmetaplotwinds),cutoffs)
  psites%>%mapToTranscripts(stopmetaplotwinds%>%subsetByOverlaps(chrpsites))%>%coverage%>%as.matrix
})

stopmeatplotdf<-stopmetaplotmats%>%setNames(names(bams))%>%map_df(.id='sample',as.data.frame)%>%gather(pos,count,-sample)
stopmeatplotdf$pos%<>%as_factor%>%as.numeric




condition = names(bams)%>%setNames(str_extract(.,'mock|si'),.)
#now plot
plotfile<- here(paste0('plots/stop_','metaplot','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
stopmeatplotdf%>%
  group_by(sample,pos)%>%summarise(count=sum(count))%>%
  mutate(condition = condition[sample])%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 1+windforw+c(-30,0,30,60,90,120),labels=c('-120bp','-90bp','-60bp','-30bp','Stop','30bp'))+
  scale_y_continuous(paste0('Psite Count'))+
  ggtitle(paste0('Psite distribution Metaplot'))+
  facet_grid(condition ~ . )+
  theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
plotfile<- here(paste0('plots/stop_','metaplot_codbin','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
stopmeatplotdf%>%
  mutate(pos = pos - (pos - 1+windforw)%%3 )%>%
  group_by(sample,pos)%>%summarise(count=sum(count))%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  mutate(condition = condition[sample])%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 1+windforw+c(-30,0,30,60,90,120),labels=c('-120bp','-90bp','-60bp','-30bp','Stop','30bp'))+
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
plotfile<- here(paste0('plots/stop_','metaplot_cod_targetsep','.pdf'))
dirname(plotfile)%>%dir.create
pdf(plotfile)
stopmeatplotdf%>%group_by(sample,pos)%>%mutate(cds=1:n())%>%
  mutate(isindsp = bestcds@elementMetadata$isindsp[cds])%>%
  mutate(is_target = bestcds@elementMetadata$is_target[cds])%>%
  filter(isindsp)%>%
  ungroup%>%
  mutate(pos = pos - (pos - 1+windforw)%%3 )%>%
  group_by(sample)%>%mutate(count = count/sum(count))%>%
  group_by(sample,is_target,pos)%>%summarise(count=sum(count))%>%
  group_by(sample,is_target)%>%mutate(count = count/sum(count))%>%
  mutate(condition = condition[sample])%>%
  ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
  geom_line()+
  scale_x_continuous(condition,breaks = 1+windforw+c(-30,0,30,60,90,120),labels=c('-120bp','-90bp','-60bp','-30bp','Stop','30bp'))+
  scale_y_continuous(paste0('Psite Count'))+
  ggtitle(paste0('Psite distribution Metaplot'))+
  facet_grid(is_target ~ . )+
  theme_bw()
dev.off()
normalizePath(plotfile)



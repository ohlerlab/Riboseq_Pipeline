library(GenomicFeatures)
source(here::here('src/Rprofile.R'))
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200,comps=c('chrM'='chrM')) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  mcols(reads)$length <- qwidth(reads)
  #
  if(is.null(offsets)){
    mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    useqnms <- as.character(unique(seqnames(reads)))%>%setdiff(names(comps))

    compmap = safe_hashmap(c(useqnms,names(comps)),c(rep('nucl',length(useqnms)),comps))
    stopifnot(all(compmap$values()%in%offsets$comp))
    
    mcols(reads)$offset <-
      data.frame(length=mcols(reads)$length,
        compartment=compmap[[as.character(seqnames(reads))]])%>%
      safe_left_join(offsets,allow_missing=TRUE)%>%.$offset
  }
  #
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  #
  # mcols(reads)$length <- width(reads)
  reads%<>%subset(!is.na(offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  mcols(psites)$length <- mcols(reads)$length
  psites
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


#NOTE - this was the end of trim_cds before, might have dependencies.

#extend exons so we can always make our windows

cutoffs <- Sys.glob('pipeline/riboseqc/data/*/_P_sites_calcs')%>%map_df(fread)%>%
  group_by(read_length,comp)%>%
  filter(read_length-cutoff > 3)%>%
  summarise(cutoffs = list(as.data.frame(table(cutoff))))%>%
  unnest(cutoffs)%>%
  mutate(cutoff = as.character(cutoff)%>%as.numeric)%>%
  arrange(desc(Freq))%>%
  slice(1)%>%
  rename(offset:=cutoff,length:=read_length)

bakcds <- cds
cds <- bakcds
bams <- Sys.glob(here('pipeline/star/data/*/*ribo*.bam'))%>%setNames(.,basename(dirname(.)))
bamseqnames <- seqinfo(Rsamtools::BamFile(bams[1]))@seqnames
sample = names(bams)[1]
chr = 'chrM'

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
countdfs <- Sys.glob('pipeline/cdscounts/*/*trimcounts*')%>%
  setNames(.,basename(dirname(.)))%>%
  map_df(.id='sample',fread)
highcountgenes <- countdfs%>%group_by(gene_id)%>%filter(any(totalcount>32))%>%.$gene_id
bestcds <- countdfs%>%filter(gene_id%in%highcountgenes)%>%group_by(gene_id)%>%slice(which.max(trimcount))%>%.$protein_id
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


metaplotmats_lensep <- mclapply(mc.cores=20,bams,function(bam){
  psites <- get_genomic_psites(bam,unlist(metaplotwinds)%>%reduce,cutoffs)
  message('psites read, creating matrices for rls')
  lapply(unique(psites$length)%>%setNames(.,.),function(ilength){
    psites%>%.[.$length==ilength]%>%
      mapToTranscripts(metaplotwinds%>%subsetByOverlaps(psites))%>%
        coverage%>%
        {set_rownames(as.matrix(.),names(.))}
  })
})

metaplotdf_ls <- metaplotmats_lensep%>%map_df(.id='sample',.%>%map_df(.id='length',.%>%as.data.frame%>%rownames_to_column('protein_id')%>%gather(pos,count,-protein_id)))
metaplotdf_ls$pos%<>%as_factor%>%as.numeric

meatplotdf <- metaplotdf_ls%>%group_by(sample,pos,protein_id)%>%summarise(count=sum(count))

# metaplotmats <- mclapply(bams,function(bam){
#   psites <- get_genomic_psites(bam,unlist(metaplotwinds),cutoffs)
#   psites%>%mapToTranscripts(metaplotwinds%>%subsetByOverlaps(psites))%>%coverage%>%{set_rownames(as.matrix(.),names(.))}
# })

# meatplotdf<-metaplotmats%>%setNames(names(bams))%>%map_df(.id='sample',as.data.frame)%>%gather(pos,count,-sample)
# meatplotdf$pos%<>%as_factor%>%as.numeric



## Get the categories for the metaplots
library(GO.db)
refseqdf <- fread(here('pipeline/ensembl2reseq.tsv'))
catfiles <- Sys.glob(here('ext_data/manu_gene_cats/*'))
godb <- AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))
go_termtext<-tibble(file=catfiles)%>%
  mutate(GOID=file%>%basename%>%str_extract('GO_\\d+')%>%str_replace("_",':'))%>%
  left_join(godb)%>%
  {setNames(.$TERM,.$GOID)}
go_termtext[]%<>%str_replace('via.*','via..')
names(go_termtext) = basename(catfiles)
go_termtext<-ifelse(is.na(go_termtext),names(go_termtext)%>%str_replace('.txt',''),go_termtext)
#
cats <- catfiles%>%setNames(.,go_termtext)%>%map_df(.id='cat',fread,header=F)%>%set_colnames(c('cat','refseq_mrna'))%>%
  mutate(refseq_mrna = str_replace(refseq_mrna,'\\.\\d+',''))
#
cats%<>%left_join(refseqdf)%>%rename('gene_id':=ensembl_gene_id)
#
catlist = cats%>%{split(.[[3]],.[[1]])}
catlist[['all']]=fmcols(bestcds,gene_id)%>%str_replace('\\.\\d+$','')

str(catlist)

condition = names(bams)%>%setNames(str_extract(.,'mock|si'),.)
profilebreaks <- c( -rev(seq(0,windback-1,by=30)),seq(0,windforw-1,by=30)[-1])
profilebreaks%<>%setNames(profilebreaks%>%paste0(.,'bp')%>%str_replace('^0bp','AUG'))
profilebreaks <- profilebreaks+windback+1
pid2gidmap <- safe_hashmap(fmcols(cds,protein_id),fmcols(cds,gene_id)%>%str_replace('\\.\\d+$',''))


sizemap <- sum(width(bestcds))%>%{safe_hashmap(names(.),.)}
normdf <- countdfs%>%filter(protein_id%in%names(bestcds))%>%transmute(sample,protein_id,density=totalcount/sizemap[[protein_id]])

#now plot
ovsetting <- c('overlaid','nonoverlaid')

for(catnm in names(catlist)){
for(ovnm in ovsetting){
  gs2use = catlist[[catnm]]

#facet or not
 facetfun =  if(ovnm=='nonoverlaid') facet_grid( condition~ . ) else NULL
#now plot
plotfile<- here(paste0('plots/','metaplot_',ovnm,'_',catnm,'.pdf'))
pdf(plotfile)
posdf <- meatplotdf%>%
    filter(pid2gidmap[[protein_id]]%in%catlist[[catnm]])%>%
    safe_left_join(normdf)
    # filter(totalcount>32)%>%
plotdf<-posdf   %>%mutate(count = count / density)%>% group_by(sample,pos)%>% summarise(count=sum(count))%>%
    mutate(condition = condition[sample])%>%
    group_by(sample)%>%mutate(count = count/sum(count))%>%
    group_by(condition,pos)%>%summarise(count=mean(count))
p <- plotdf%>%    ggplot(.,aes(x=pos,y=count,color=condition))+
    geom_line()+
    scale_x_continuous(condition,breaks = profilebreaks)+
    scale_y_continuous(paste0('Psite Count'))+
    ggtitle(str_interp(paste0('Psite distribution Metaplot\nCategory ${catnm}\n(${n_distinct(posdf$protein_id)}) genes) ')))+
    facetfun+
    theme_bw()
print(p)
dev.off()
normalizePath(plotfile)%>%message

}
}

# normalizePath(plotfile)
################################################################################
########Messing with nbPCA, models, tf and other such nonsense
################################################################################
reticulate::use_python("~/work/miniconda3/envs/rtensorflow/bin/python")
reticulate::use_condaenv("rtensorflow")
# tensorflow::install_tensorflow(version=1.4)

source('Applications/scNBMF/scNBMF.R')
source('src/scNBMF.R')
#BiocManager::install(c('igraph'))

profs <- matrix(c(c(1,1,1,1,1,1,1,1,1),
c(1,1,1,1,5,1,1,1,1)),ncol=2) 

bmixes <- t(matrix(byrow=T,ncol=2,c(3,0,4,0,5,0,5,2,3,2)))
mixes=list(bmixes)%>%rep(100)%>%do.call(what=cbind)

simnbpcadat <- t(mixes) %*% t(profs)

write.table(col.names=F,row.names=F,t(simnbpcadat),sep=',','tmp.tsv',quote=F)
scNBMF(input_file='tmp.tsv',batch_size=50,num_iter=500,ndim=3)

x <- sim.negbin(c(4,5,10),3,10,12)

library(glmpca)
res<-glmpca(simnbpcadat,9)
txtplot(res$factors[[1]])
txtplot(res$factors[[2]])

simnbpcadat%>%dim
res[[1]]$loadings%>%dim

vsd((simnbpcadat))


highexprprots<-countdfs%>%filter(sample==sample[1])%>%arrange(desc(totalcount))%>%slice(1:500)%>%.$protein_id

posdf <- meatplotdf%>%
    filter(protein_id%in%highexprprots)%>%
    safe_left_join(normdf)%>%
    filter(sample==sample[1])

splglm <- glm.nb(data=posdf,count ~ as_factor(pos) + density ,link='log')

library(txtplot)

#now plot
plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)
predict(newdata=posdf%>%ungroup%>%distinct(pos)%>%mutate(density=1),splglm,se.fit=TRUE)%>%
  {data.frame(l=.$fit - 1.96*.$se.fit ,est=.$fit, u = .$fit + 1.96*.$se.fit ,pos=posdf$pos%>%unique)}%>%
  ggplot(.,aes(ymin=l,ymax=u,y=est,x=pos))+
  geom_line()+
  geom_ribbon(alpha=I(0.2))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('pos'))+
  scale_y_continuous(paste0('negative binomial prediction Rprofile.R'))+
  ggtitle(paste0('title'))+
  theme_bw()
dev.off()
normalizePath(plotfile)

plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)
  posdf%>%ungroup%>%mutate(resid=splglm$residuals)%>%
  ggplot(.,aes(y=resid,x=pos,group=protein_id))+
  geom_line(alpha=I(0.1))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('pos'))+
  scale_y_continuous(paste0('negative binomial prediction Rprofile.R'))+
  ggtitle(paste0('title'))+
  theme_bw()
dev.off()
normalizePath(plotfile)


## smoothing codons

# #now plot
# plotfile<- here(paste0('plots/','metaplot_codbin','.pdf'))
# dirname(plotfile)%>%dir.create
# pdf(plotfile)
# meatplotdf%>%
#   mutate(pos = pos - (pos - 51)%%3 )%>%
#   group_by(sample,pos)%>%summarise(count=sum(count))%>%
#   group_by(sample)%>%mutate(count = count/sum(count))%>%
#   mutate(condition = condition[sample])%>%
#   ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
#   geom_line()+
#   scale_x_continuous(condition,breaks = 51+c(-30,0,30,60,90),labels=c('-30bp','AUG','30bp','60bp','90bp'))+
#   scale_y_continuous(paste0('Psite Count'))+
#   ggtitle(paste0('Psite distribution Metaplot'))+
#   facet_grid(condition ~ . )+
#   theme_bw()
# dev.off()
# normalizePath(plotfile)

# dspdf <- fread('ext_data/no_DSP.csv')
# dspdf%<>%filter(CDS_reads_IP>32 &(CDS_reads_total>32))

# bestcds@elementMetadata$isindsp <- fmcols(bestcds,gene_name) %in% dspdf$Gene_name
# bestcds@elementMetadata$is_target <- fmcols(bestcds,gene_name) %in% (dspdf%>%filter(RPKM_ratio 1.5)%>%.$Gene_name)

## Seperating out TE targets

# #now plot
# plotfile<- here(paste0('plots/','metaplot_cod_targetsep','.pdf'))
# dirname(plotfile)%>%dir.create
# pdf(plotfile)
# meatplotdf%>%group_by(sample,pos)%>%mutate(cds=1:n())%>%
#   mutate(isindsp = bestcds@elementMetadata$isindsp[cds])%>%
#   mutate(is_target = bestcds@elementMetadata$is_target[cds])%>%
#   filter(isindsp)%>%
#   ungroup%>%
#   mutate(pos = pos - (pos - 51)%%3 )%>%
#   group_by(sample)%>%mutate(count = count/sum(count))%>%
#   group_by(sample,is_target,pos)%>%summarise(count=sum(count))%>%
#   group_by(sample,is_target)%>%mutate(count = count/sum(count))%>%
#   mutate(condition = condition[sample])%>%
#   ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
#   geom_line()+
#   scale_x_continuous(condition,breaks = 51+c(-30,0,30,60,90),labels=c('-30bp','AUG','30bp','60bp','90bp'))+
#   scale_y_continuous(paste0('Psite Count'))+
#   ggtitle(paste0('Psite distribution Metaplot'))+
#   facet_grid(is_target ~ . )+
#   theme_bw()
# dev.off()
# normalizePath(plotfile)


#this checks that the stop codon is in fact included in the cds
#if it isn't, you'll want to 
stop_in_cds <- bestcds%>%sample(20)%>%resize_grl(3,'end')%>%strandshift(0)%>%
  unlist%>%subset(width==3)%>%
  {getSeq(FaFile(fafile),.)}%>%translate%>%
  as.character%>%is_in('*')
stopifnot(all(stop_in_cds)|all(!stop_in_cds))
stopcodonshift = if(all(stop_in_cds)) -2 else 1


windback <- 50
windforw <- 150
stopexonsexp <- exons[fmcols(bestcds,transcript_id)]%>%
  resize_grl(sum(width(.))+windback,'start',check=FALSE)%>%
  resize_grl(sum(width(.))+windforw,'end',check=FALSE)%>%
  .[!any(is_out_of_bounds(.))]

stopmetaplotwinds <- pmapToTranscripts (bestcds%>%resize_grl(1,'end'),stopexonsexp[fmcols(bestcds,transcript_id)])%>%
    unlist%>%
    shift(stopcodonshift)%>%
    resize(windforw +1,fix='end',ignore.strand=TRUE)%>%
    resize(windforw+ windback,'start',ignore.strand=TRUE)%>%
    {spl_mapFromTranscripts(.,exons_grl = stopexonsexp[seqnames(.)])}
stopmetaplotwinds%<>%split(.,names(.))

stopmetaplotmats <- mclapply(bams,function(bam){
  psites <- get_genomic_psites(bam,unlist(stopmetaplotwinds),cutoffs)
  psites%>%mapToTranscripts(stopmetaplotwinds%>%subsetByOverlaps(psites))%>%coverage%>%{set_rownames(as.matrix(.),names(.))}
})

stopmeatplotdf<-stopmetaplotmats%>%setNames(names(bams))%>%map_df(.id='sample',as.data.frame)%>%gather(pos,count,-sample)
stopmeatplotdf$pos%<>%as_factor%>%as.numeric

condition = names(bams)%>%setNames(str_extract(.,'mock|si'),.)




stopprofilebreaks <- c( -rev(seq(0,windforw-1,by=30)),seq(0,windback-1,by=30)[-1])
stopprofilebreaks%<>%setNames(stopprofilebreaks%>%paste0(.,'bp')%>%str_replace('^0bp','Stop'))
stopprofilebreaks <- stopprofilebreaks+windforw+1

profilebreaks <- c( -rev(seq(0,windback-1,by=30)),seq(0,windforw-1,by=30)[-1])
profilebreaks%<>%setNames(profilebreaks%>%paste0(.,'bp')%>%str_replace('^0bp','AUG'))
profilebreaks <- profilebreaks+windback+1

profbrklist <- list(start=profilebreaks,stop=stopprofilebreaks)




#now plot
ovsetting <- c('overlaid','nonoverlaid')
for(ovnm in ovsetting){
#facet or not
 facetfun =  if(ovnm=='nonoverlaid') facet_grid( condition~ . ) else NULL
#now plot
plotfile<- here(paste0('plots/','stopmetaplot_',ovnm,'.pdf'))
pdf(plotfile)
stopmeatplotdf%>%
    group_by(sample,pos)%>%summarise(count=sum(count))%>%
    mutate(condition = condition[sample])%>%
    group_by(sample)%>%mutate(count = count/sum(count))%>%
    ggplot(.,aes(x=pos,y=count,color=condition,group=sample))+
    geom_line()+
    scale_x_continuous(condition,breaks = profilebreaks)+
    scale_y_continuous(paste0('Psite Count'))+
    ggtitle(paste0('Psite distribution Metaplot'))+
    facetfun+
    theme_bw()
dev.off()
normalizePath(plotfile)%>%message
}








fafile=file.path(here('pipeline',yaml::yaml.load_file('src/config.yaml')%>%.$REF_orig%>%basename%>%str_replace('.gz$','')))%T>%{stopifnot(file.exists(.))}
library(Rsamtools)

getSeq(FaFile(fafile),stopmetaplotwinds%>%unlist%>%subset(width==3))

testatgs%>% resize(1)%>%{pmapToTranscripts (.,exonsexp[testatgs$transcript_id])}%>%
  # resize(3,fix='end',ignore.strand=TRUE)%>%
  resize(3,'start',ignore.strand=TRUE)%>%
  {spl_mapFromTranscripts(.,exons_grl = exonsexp[seqnames(.)])}%>%{getSeq(FaFile(fafile),.)}

metaplotwinds%>%.[elementNROWS(.)==1]%>%resize(150,'end')%>%resize(3,'start')%>%head%>%unlist%>%getSeq(FaFile(fafile),.)
metaplotwinds


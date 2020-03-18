library(tidyverse)
library(GenomicAlignments)
library(Rsamtools)
library(data.table)
library(magrittr)
library(Matrix)
library(qlcMatrix)
library(txtplot)
library(assertthat)

acgtn <- c('A','C','T','G','N')

##read our design file and load necessary parameters
argv=here::here('src/config.yaml')%T>%{stopifnot(file.exists(.))}
design_file <- commandArgs(trail=TRUE)[1]
if((length(design_file)==0) | is.na(design_file) ) design_file<-normalizePath(file.path('src/config.yaml'),mustWork=TRUE) 
message(paste0("Reading the design file at ",design_file,'\n\n'))
#now load the and meta data
design_list <- yaml::yaml.load_file(design_file)
#and put the parameters from our design file in the global environment
parameters <- design_list$parameters
assert_that(file.exists(parameters$root))
for(param in names(parameters)) {
  message(paste0('Assigning ',param,' = ',parameters[[param]]))
  assign(param,parameters[[param]])
}
message('\n\n')

seq <- getSeq(Rsamtools::FaFile('pipeline/masked_fasta/masked.fa'))

#load annotation
anno<-rtracklayer::import(design_list$GTF_orig) 
exons<-anno%>%subset(type=='exon')
exongrl<-exons%>%split(.,.$gene_id)

letterFrequency(seq[gaps(exons[1:10])],acgtn,as.prob=T)
letterFrequency(seq[identity(exons[1:10])],acgtn,as.prob=T)

exont<-rtracklayer::import('pipeline/masked_fasta/exons.gtf')
exont<-rtracklayer::import('pipeline/masked_fasta/exons.gtf')

oligonucleotideFrequency(seq[identity(exons[1:10])],1)
?oligonucleotideFrequency

Biostrings::oligonucleotideFrequency(seq[[1]],1,as.prob=T,alphabet=acgtn)

bam = 'pipeline/star/data/RNA_LEP/RNA_LEP.bam'
ribobam = 'pipeline/star/data/RPF_LEP/RPF_LEP.bam'

mapqthresh=200


reads = readGAlignments(BamFile(bam),param=ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),tagFilter=list('NH'=1),what="rname",mapqFilter=mapqthresh),use.names=T)
unreads = readGAlignments(BamFile(bam),param=ScanBamParam(tagFilter=list('NH'=1),what="rname",mapqFilter=mapqthresh),use.names=T)
reads = readGAlignments(BamFile(bam),param=ScanBamParam(tagFilter=list('NH'=2:10),what="rname"),use.names=T,mapqFilter=mapqthresh)



fcounts<-fread('pipeline/feature_counts/all_feature_counts')
bam

fcounts[,'RNA_LEP']%>%sum/1e6
fcounts[,'RPF_LEP']%>%sum/1e6



exongrl%>%countOverlaps(reads,ignore.strand=F)%>%sum

unreads%>%overlapsAny(GenomicRanges::reduce(exongrl))%>%sum%>%divide_by(1e6)

mcols(reads)$exonov <- overlapsAny(reads,exons)


readlentables<-qwidth(reads)%>%split(seqnames(reads))%>%lapply(table)

readlentables%>%lapply(enframe)%>%bind_rows(.id='chr')%>%mutate(comp=str_replace(chr,'chr\\d+','nucl'))%>%group_by(comp,name)%>%summarise(n=sum(value))%>%
	filter(name%in% c(18:40))%>%
	spread(name,n)

#the number of gnes reads align to - is it, e.g., mostly 2? Then genome duplication might be overcome by merging the duplicate genes.
rovnumtally<-data.frame(exonov = mcols(reads)$exonov,rname = names(reads))%>%group_by(rname)%>%summarise(exonhits=sum(exonov))%>%group_by(exonhits)%>%tally
rovnumtally%>%mutate(frac=n/sum(n))










ov <- findOverlaps(reads,exongrl)
ov %<>% as.data.frame
rnameints <- names(reads)%>%as.factor%>%as.numeric
ov$queryHits <- rnameints[ov$queryHits]
ov <- ov%>%arrange(queryHits,subjectHits)





head(ov)

ovmat <- Matrix(0,nrow=max(ov$queryHits),ncol=max(ov$subjectHits),sparse=TRUE)
ovmat[as.matrix(ov)]<-1
ovmat[1,3225:3226]

ovmat <- ovmat[,colSums(ovmat)>5]

length(exongrl)

sharedreadmat <- t(ovmat) %*% ovmat

diag(sharedreadmat)%>%log10%>%txtdensity

#so now it's the fraction of the row reads shared
sharedreadmatnorm <- sharedreadmat/diag(sharedreadmat)
diag(sharedreadmatnorm)<-0



rowmaxs <- sharedreadmatnorm %>% rowMax
rowmaxs  %>% as.vector %>% txtdensity

rowSums(sharedreadmatnorm > 0.7)%>%table %>%{ ./sum(.)}%>%round(3)

(sharedreadmatnorm > 0.7)%>%which%>%head(1)

sharedreadmatnorm
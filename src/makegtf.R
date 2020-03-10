library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(dplyr)
library(magrittr)
library(assertthat)
library(Rsamtools)

conflist <- yaml::yaml.load_file('src/config.yaml')
REF = FaFile(conflist$REF_orig)
gtf<-conflist$GTF_orig
annogr <- import(gtf)

assert_that('gene_id' %in% colnames(mcols(annogr)))
assert_that('protein_id' %in% colnames(mcols(annogr)))
assert_that('transcript_id' %in% colnames(mcols(annogr)))
if('gene_biotype' %in% colnames(mcols(annogr))) annogr$gene_type=annogr$gene_biotype
assert_that('gene_type' %in% colnames(mcols(annogr)))

system(str_interp('samtools faidx ${REF$path}'))

unique(seqnames(annogr))%in%seqnames(seqinfo(REF))
setdiff(unique(seqnames(annogr)),seqnames(seqinfo(REF)))
setdiff(seqnames(seqinfo(REF)),unique(seqnames(annogr)))

test <- annogr%>%subset(strand=='+')%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sample(10)

test <- test%>%lapply(head,1)%>%GRangesList%>%unlist

getSeq(REF,resize(test,3))


db <- makeTxDbFromGFF(gff)

db%>%export(gff%>%str_replace('gff(.gz)?','gtf'))



impgtf <- import(gff%>%str_replace('gff(.gz)?','gtf'))

library(Rsamtools)
REF=FaFile('/fast/work/groups/ag_ohler/dharnet_m/Mitosis_Riboseq/ext_data/ZmB73_AGPv1_genome.fasta.bgz')


getSeq(,)

allseq=getSeq(REF)

names(allseq) <- paste0('chr',allseq%>%names)
allseq%>%writeXStringSet('ext_data/Zea_mays/Ensembl/AGPv4/Sequence/WholeGenomeFasta/genome.chrnm.fa')


seqinfo(REF)
annogr%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%head%>%{extractTranscriptSeqs(x=REF,transcripts=.)}

getSeq(REF,annogr%>%subset(type=='CDS')%>%subset(protein_id==unique(protein_id)[3])%>%.[1])

allseq[annogr%>%subset(type=='CDS')%>%subset(protein_id==unique(protein_id)[1])%>%.[1]]
allseq[annogr%>%subset(type=='CDS')%>%subset(protein_id==unique(protein_id)[2])%>%tail(1)%>%.[1]]%>%reverseComplement
allseq[annogr%>%subset(type=='CDS')%>%subset(protein_id==unique(protein_id)[1])%>%.[1]]


exonsBy(db)
cdsBy(db)
transcriptsBy(db)

'TxID:53764' %in% annogr$Name


testcds <- impgtf%>%subset(type=='CDS') %>%head(20)

testcds$Parent[1]

cdsBy(db)[1]

overlapsAny(cdsBy(db)[1],)






geneids <- annogr%>%subset(type=='gene')%>%mcols%>%as.data.frame%>%dplyr::select(gene_id = ID)%>%.$gene_id
annogr$gene_id<-NA
annogr$gene_id[annogr$type=='gene'] <- geneids

annogr$type%<>%recode('mRNA'='transcript')
stopifnot ( all(all(annogr%>%subset(type=='transcript')%>%.$Parent %in% ( annogr$gene_id))))

annogr$Parent%>%elementNROWS%>%table
annogr$gene_id[annogr$type=='transcript'] <- unlist(annogr$Parent[annogr$type=='transcript'])
annogr$transcript_id<-NA
annogr$transcript_id[annogr$type=='transcript'] <- unlist(annogr$ID[annogr$type=='transcript'])

stopifnot(all(all(annogr%>%subset(type=='exon')%>%.$Parent%>%unlist %in% (annogr$transcript_id))))

isexon <- annogr$type=='exon'
annogr$exon_id<-NA
annogr$transcript_id[isexon] <- unlist(annogr$Parent[isexon])
annogr$exon_id[isexon] <- unlist(annogr$ID[isexon])


stopifnot(all(all(annogr%>%subset(type=='CDS')%>%.$Parent%>%unlist %in% (annogr$transcript_id))))

iscds <- annogr$type=='CDS'
annogr$protein_id<-NA
annogr$transcript_id[iscds] <- unlist(annogr$Parent[iscds])
annogr$protein_id[iscds] <- unlist(annogr$ID[iscds])

tr2giddf <- annogr%>%subset(type=='transcript')%>%mcols%>%as.data.frame%>%dplyr::distinct(transcript_id,gene_id)

annogr$gene_id <- data.frame(transcript_id=annogr$transcript_id)%>%left_join(tr2giddf)%>%.$gene_id


annogr%>%export(gff%>%str_replace('gff(.gz)?','gtf'))


library(GenomicFeatures)
source('src/Rprofile.R')

anno <- projmemoise(rtracklayer::import)('data/gencode.vM12.annotation.gtf')

cds <- anno%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sort_grl_st
exons <- anno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st
cdsdj <- cds%>%isDisjoint
# 
# sorted_is_dj <- function(grl){
#   starts <- start(grl)%>%{.[ IntegerList(as.list(rep(-1,length((.)))) ) ]}
#   ends = start(grl)%>%{.[ IntegerList(as.list(-elementNROWS(.)) ) ]}
#   all(all(starts < ends))
# }
# sorted_is_dj(cds)

exprgenes <- counts(dds)%>%rowSums%>%add(1)%>%log10 %>%`>`(2)%>%names

cds <- cds[split(cds@unlistData@elementMetadata$gene_id %in% exprgenes,cds@partitioning)]


maxlen_pids <- sum(width(cds))%>%enframe('protein_id','len')%>%
    left_join(data.frame(mcols(cds@unlistData)[,c('protein_id','gene_id')])%>%distinct,by='protein_id')%>%
    group_by(gene_id)%>%slice(which.max(len))%>%
    .$protein_id


pmapToTranscripts(mcds%>%resize_grl(1),exons[fmcol(mcds,transcript_id)]%>%trim_grl(-100,end='fp'))%>%
    resize(100,'start')%>%
    resize(200,'end')%>%
  

fmcol <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcol(mcds,transcript_id)



debugger()

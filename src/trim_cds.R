#load the gtf
library(rtracklayer)
gtf = here('../cortexomics/my_gencode.vM12.annotation.gtf')
STARTTRIMCODS <- 15
ENDTRIMCODS <- 5
if(!exists('gtf_gr')) gtf_gr<-import(con=gtf,format='gtf')
# gtf_grbak <- gtf_gr
# gtf_gr<-gtf_grbak

# gtf_gr <- gtf_grbak[gtf_grbak$transcript_id%in%unique(gtf_grbak$transcript_id)[666:676],]
# gtf_gr<-sort(gtf_gr)

exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')



trimlen = (STARTTRIMCODS*3) + (STARTTRIMCODS*3)
#define longe enough ids
longenough_ids = cds%>%split(.,.$protein_id)%>%width%>%sum%>%.[.>trimlen]%>%names
#now subset
cds %<>% subset(protein_id %in% longenough_ids)
exons %<>% subset(transcript_id %in% unique(cds$transcript_id))
genes <- gtf_gr %>% subset(type=='gene') %>% subset(gene_id %in% exons$gene_id)

#select some cds regions and corresponding transciropts
testcds <- cds%>%split(.,.$protein_id)
testtrs <- exons%>%subset(transcript_id%in%unlist(testcds)$transcript_id)%>%split(.,.$transcript_id)
#get start codons
startcods <- testcds%>%revElements(.,any(strand(.)=='-'))%>%lapply('[',1)%>%GRangesList%>%unlist%>%resize(1,'start')
endcods <- testcds%>%revElements(.,any(strand(.)=='+'))%>%lapply('[',1)%>%GRangesList%>%unlist%>%resize(1,'end')
#map our start codons to transcript space
library(GenomicFeatures)
matching_trs <- testtrs[startcods$transcript_id]
startcods_trspace <- pmapToTranscripts(startcods,matching_trs)
endcods_trspace <- pmapToTranscripts(endcods,matching_trs)
#this function maps granges in transcript space back to the genome, poentially creating more
#than n ranges for n in put ranges, if they cross exon boundaries
spl_mapFromTranscripts<-function(trspacegr,exons_grl){
  exons_tr<-exons_grl%>%unlist%>%setNames(paste0('exon_',seq_along(.)))%>%mapToTranscripts(exons_grl)
  ov <- findOverlaps(trspacegr,exons_tr)
  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
  trspacegr_spl,
  # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
  exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}

cds_trspace <- startcods_trspace
end(cds_trspace) <- end(endcods_trspace)

start(cds_trspace) <- start(cds_trspace) + (STARTTRIMCODS*3)
end(cds_trspace) <- end(cds_trspace) - (ENDTRIMCODS*3)

trimmed_cds <- cds_trspace%>%
  spl_mapFromTranscripts(testtrs)%>%#map back to the genome
  split(.,.$xHits)%>%
  unlist(use.names=TRUE)

mcols(trimmed_cds) <- testcds[trimmed_cds$xHits]%>%.[IntegerList(as.list(rep(1,length(.))))]%>%unlist%>%
  {mcols(.)[,c('gene_id','transcript_id','protein_id','gene_name','transcript_name','type','phase','tag','gene_type')]}
names(trimmed_cds)<-NULL

trimmed_cds%>%export(gtf%>%str_replace('.gtf$','.trimmed.gtf'))

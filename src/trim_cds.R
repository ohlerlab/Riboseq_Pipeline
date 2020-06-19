#load the gtf
source('../src/Rprofile.R')
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}

#load copy of the gtf in the pipeline folder
#gtf_file = here::here('pipeline',basename(yaml::yaml.load_file('../src/config.yaml')$GTF_orig))%>%{stopifnot(file.exists(.))}
gtf_file = here::here('pipeline',basename(yaml::yaml.load_file('../src/config.yaml')$GTF_orig))

anno <- projmemoise(function(...){rtracklayer::import(...)})(gtf_file)
#We want a function that 
allcds <- anno%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sort_grl_st
exons <- anno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st 

STARTTRIMCODS <- yaml::yaml.load_file('../src/config.yaml')$STARTCODTRIM%T>%{stopifnot(!is.null(.))}
ENDTRIMCODS <- yaml::yaml.load_file('../src/config.yaml')$STOPCODTRIM%T>%{stopifnot(!is.null(.))}
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





# metaplotmats

# fafile=file.path(here('pipeline',yaml::yaml.load_file('src/config.yaml')%>%.$REF_orig%>%basename%>%str_replace('.gz$','')))%T>%{stopifnot(file.exists(.))}

# testatgs <- cds[bestcds]%>%head%>%resize_grl(3,'start')%>%.[elementNROWS(.)==1]%>%head%>%unlist

# getSeq(FaFile(fafile),testatgs)

# testatgs%>% resize(1)%>%{pmapToTranscripts (.,exonsexp[testatgs$transcript_id])}%>%
#   # resize(3,fix='end',ignore.strand=TRUE)%>%
#   resize(3,'start',ignore.strand=TRUE)%>%
#   {spl_mapFromTranscripts(.,exons_grl = exonsexp[seqnames(.)])}%>%{getSeq(FaFile(fafile),.)}

# metaplotwinds%>%.[elementNROWS(.)==1]%>%resize(150,'end')%>%resize(3,'start')%>%head%>%unlist%>%getSeq(FaFile(fafile),.)
# metaplotwinds




suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,DESeq2))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,xtail))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,here))
select<-dplyr::select
args <- c(
  countfile=here('pipeline/feature_counts/all_feature_counts'),
  # uORFcountfile='SaTAnn/uORFs.feature_counts',
  outdir= here('pipeline/xtail')
)
CONDITIONVAR = 'fraction'

args[] = commandArgs(trailingOnly=TRUE)
for(i in names(args)) assign(i,args[i])

featurecountsagg <- data.table::fread(countfile)
# uorfcounts <- fread(uORFcountfile)

samplepars = fread(here('src/sample_parameter.csv'))

samples = fread(here('src/sample_parameter.csv'))%>%
  # filter(is.na(fraction))%>%
  filter(!str_detect(sample_id,'test'))%>%
  .$sample_id
featurecountsagg = featurecountsagg%>%select(feature_id,one_of(samples))
# uorfcounts = uorfcounts%>%select(feature_id,one_of(samples))

# featurecountsagg <- featurecountsagg%>%rbind(uorfcounts)

    #samples to use
xtailcounts =   featurecountsagg
xtailcounts = xtailcounts[!apply(xtailcounts,1,.%>%is.na%>%any),]


rnasamples<-str_subset(samplepars$sample_id,'rnaseq')
ribosamples<-str_subset(samplepars$sample_id,'ribo')

conditionvect <-  samplepars%>%filter(sample_id %in% rnasamples)%>%.[[CONDITIONVAR]]
conditionvect_ribo <-  samplepars%>%filter(sample_id %in% ribosamples)%>%.[[CONDITIONVAR]]
stopifnot(conditionvect == conditionvect_ribo)
conditions=unique(conditionvect)
stopifnot(length(conditions)==2)

mrnatab <- xtailcounts%>%as.data.frame%>%.[,rnasamples]%>%set_rownames(xtailcounts$feature_id)
ribotab <- xtailcounts%>%as.data.frame%>%.[,ribosamples]%>%set_rownames(xtailcounts$feature_id)

xtailres <- xtail::xtail(mrnatab,ribotab,conditionvect,threads=20)

xtailtable <- xtailres$resultsTable%>%rownames_to_column%>%set_colnames(
  c("feature_id","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", paste0(conditions[1],"_log2TE"), 
  paste0(conditions[2],"_log2TE"), "log2FC_TE_v2", "pvalue_v2", "log2fc", 
  "p_value", "adj_p_value")
)

xtailfile = file.path(outdir,paste0('xtail_',CONDITIONVAR,'.tsv'))

sizefacts <- xtailcounts%>%select(-feature_id)%>%{DESeq2::estimateSizeFactorsForMatrix(.)}

basemeandf<-sweep(xtailcounts%>%select(-feature_id),MARGIN=2,STAT=sizefacts,FUN='/')%>%rowMeans%>%setNames(xtailcounts$feature_id)%>%
  enframe('feature_id','base_mean')

xtailtable%<>%left_join(basemeandf)%>%select(feature_id,base_mean,everything())


xtailtable$padj<-xtailtable$adj_p_value
xtailtable$log2FoldChange<-xtailtable$log2fc
xtailtable$log2fc_se<-NA


write_tsv(xtailtable,xtailfile)


# save.image(file.path('prepimage.R'))


#let's make some toy images. Genes



# colnames(xtailres)%<>%str_replace_all('[\\(\\s\\)]','_')
# 

# annotation_gr <- rtracklayer::import("~/bih_cluster/projects/cubit/current/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3")
# 
# ribodiffres<-left_join(
#   ribodiffres,
#   annotation_gr%>%mcols%>%as.data.frame%>%select(gene_id,gene_name)%>%distinct(gene_id,gene_name),
#   by=c('geneIDs'='gene_id')
# )%>%  select(-disper,-pval,-X8)%>%
#   select(gene_name,everything())
# 
# ribodiffressig = ribodiffres %>% filter(padj<0.05,!is.nan(padj))
# 
# ribodiffressig%>%arrange(log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# ribodiffressig%>%arrange(-log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# 
# ribodiffres %<>% mutate(l2fc_te = `log2FC_TE(P0 vs E13)`)
# ribodiffres %<>% mutate(log10_mean_te = log10(TEE13+TEP0) / 2)
# 
# ribodiffres$aveLogCPM <-
#   edgeR::aveLogCPM(ribodiffcounts%>%select(-Entry)%>%as.matrix%>%as.integer)%>%
#   setNames(ribodiffcounts$Entry)%>%
#   .[ribodiffres$geneIDs]
#   
# ribodiffres%>%{qplot(data=.,x =aveLogCPM,y=l2fc_te,color=padj<0.05)}

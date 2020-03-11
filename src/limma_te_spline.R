
library(limma)
library(DESeq2)


isribo <- colData(dds)$assay=='ribo'
ribocounts <- counts(dds)[,isribo]
rnacounts <- counts(dds)[,!isribo]

riboislow = apply(ribocounts,1,function(x) all(x==0))
rnaislow = apply(rnacounts,1,function(x) all(x==0))

eitherlow = riboislow | rnaislow


#ribolib <- 
# voom%>%args
ribovoom <- voom(ribocounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=colData(dds)[isribo,]))

rnavoom <- voom(rnacounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=colData(dds)[isribo,]))

allvoom <- voom(counts(dds)[!eitherlow,],design = model.matrix(~ assay*ns(nstage,1),data=colData(dds)[,]))

allvoom$weights[,isribo] <- ribovoom$weight
allvoom$weights[,!isribo] <- rnavoom$weight


allvoomfit <- lmFit(allvoom)
allvoom_ebayes<-eBayes(allvoomfit)






tecontrastvect = allvoom_ebayes$coef%>%colnames%>%str_detect('assayribo:')
allvoom_ebayes$coef%>%colnames
colnames(allvoom_ebayes)[tecontrastvect]

# topTable(fit.contrast(allvoom_ebayes,tecontrastvect),coef = tecontrastvect,n=Inf)%>%head

topTable(eBayes(contrasts.fit(allvoom_ebayes,tecontrastvect)),n=Inf)%>%.$adj.P.Val%>%`<`(0.05)%>%table

topTable(eBayes(contrasts.fit(allvoom_ebayes,tecontrastvect)),n=Inf)%>%.$adj.P.Val%>%`<`(0.05)%>%table

topTable(eBayes(contrasts.fit(allvoom_ebayes,tecontrastvect)),n=Inf)%>%.$adj.P.Val%>%`<`(0.05)%>%table




topTable(allvoom_ebayes,coef = tecontrastvect,n=Inf)%>%head


resultsNames(dds)

 dds <- '/fast/work/groups/ag_ohler/dharnet_m/Mitosis_Riboseq/pipeline/run_rseq/subanalyses/All/TE_change_linear/intermediate_results/deseq2_object.rds'%>%readRDS
counts(dds)
 
results(dds[apply(counts(dds),1,median)>0,],as.numeric(tecontrastvect))%>%subset(padj<0.05)

res <- results(dds,contrast = as.numeric(tecontrastvect))

res%>%as.data.frame%>%filter(log2FoldChange > 5)

res%>%subset(padj<0.05)


topTable(eBayes(contrasts.fit(lmFit(allvoom),tecontrastvect)),n=Inf)%>%.$adj.P.Val%>%`<`(0.05)%>%table



library(xtail)



xtail%>%args



library(xtail)


xtail_m2<-xtail(ribocounts[!eitherlow,],rnacounts[!eitherlow,],condition=ifelse(colnames(ribocounts)%>%str_detect('MII'),'MII','Other'))




xtail_m2[[1]]%>%subset(pvalue.adjust<0.05)


library(txtplot)


ribocounts%>%rowMedians%>%add(1)%>%log10%>%txtdensity



################################################################################
########Spline to linear space
################################################################################
ribovoom <- voom(ribocounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=colData(dds)[isribo,]))
rnavoom <- voom(rnacounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=colData(dds)[isribo,]))
allvoom <- voom(counts(dds)[!eitherlow,],design = model.matrix(~ assay*ns(nstage,3),data=colData(dds)[,]))

allvoomfit <- lmFit(allvoom)

# itime_modelmat<-model.matrix(~ 1 + assay+assay:(ribo+MS),data=mscountvoomdesign%>%mutate(time=factor(tps[add(time,3)])))
# itime_modelmat%<>%set_colnames(itime_modelmat%>%colnames%>%str_replace('riboTRUE','TE')%>%str_replace('MSTRUE','MS_dev'))

SPLINE_N = allvoomfit$coef%>%colnames%>%str_extract('ns\\(.*?\\)')%>%str_extract('\\d+(?=\\))')%>%tail(1)%>%as.numeric

#for transforming back to time space from spline space
splinemat <- t(ns(1:5,SPLINE_N))
splinezeros <- matrix(0,ncol=5,nrow=SPLINE_N)

voomeffects<-allvoomfit$coef%>%colnames
ntps=5

effect_is_time <- allvoom_ebayes$coef%>%colnames%>%str_detect(negate=TRUE,'ns\\(')
stopifnot(length(rle(effect_is_time))==2)#all main effects, then time
nosplinen <- voomeffects%>%str_detect(negate=TRUE,'ns\\(')%>%which%>%tail(1)
maineffzeros <- matrix(0,ncol=ntps,nrow=nosplinen)
tps = colData(dds)$type%>%unique

#contrast matrices transform back to time space
alltimeeff <- rbind(
	maineffzeros,
	t(ns(1:5,SPLINE_N)),
	splinezeros
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('all_',tps))
#
timeTEeffect <- rbind(
	maineffzeros,
	splinezeros,
	t(ns(-2:2,SPLINE_N))
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('TE_',tps))

ntps<-length(tps)
sntps <- seq_along(tps)
i_n <- 1

itimecontrasts=alltimeeff

stepwise_contrasts <- lapply(list(all=alltimeeff,TE=timeTEeffect),function(itimecontrasts){
	out = lapply(1:(ntps),function(i_n){
		nexttp = (i_n%%ntps)+1
		cname <- colnames(itimecontrasts)[nexttp]%>%str_replace('_','from_')
		(itimecontrasts[,nexttp,drop=FALSE] - (itimecontrasts[,i_n,drop=FALSE])) %>%
		set_colnames(cname)
	})%>%do.call(cbind,.)
})%>%do.call(cbind,.)

((1:5) %% 5)+1

allvoom

itime <- eBayes(contrasts.fit(allvoomfit,timeTEeffect[,-1][,]))
itime%>%topTable(number=Inf)%>%rownames_to_column%>%mutate(sig = adj.P.Val<0.05)%>%group_by(sig)%>%tally

stepwisespline <- eBayes(contrasts.fit(allvoomfit,stepwise_contrasts[,-c(1:5)]))

stepwisespline%>%topTable(number=Inf)%>%rownames_to_column%>%mutate(sig = adj.P.Val<0.05)%>%group_by(sig)%>%tally

stepwisespline$coef%>%head
stepwisespline$coef%>%princomp%>%summary
stepwisespline$coef%>%princomp%>%.$loadings

splinedds<-DESeq(DESeqDataSet(dds,design = ~ assay*ns(nstage, 3)))

results(splinedds,stepwise_contrasts[,5+1])%>%.$padj%>%`<`(0.05)%>%table
results(splinedds,stepwise_contrasts[,5+2])%>%.$padj%>%`<`(0.05)%>%table
results(splinedds,stepwise_contrasts[,5+3])%>%.$padj%>%`<`(0.05)%>%table
results(splinedds,stepwise_contrasts[,5+4])%>%.$padj%>%`<`(0.05)%>%table
results(splinedds,stepwise_contrasts[,5+5])%>%.$padj%>%`<`(0.05)%>%table

# results(splinedds,timeTEeffect[,-1][,1])%>%.$padj%>%`<`(0.05)%>%table
# results(splinedds,timeTEeffect[,-1][,2])%>%.$padj%>%`<`(0.05)%>%table
# results(splinedds,timeTEeffect[,-1][,3])%>%.$padj%>%`<`(0.05)%>%table
# results(splinedds,timeTEeffect[,-1][,4])%>%.$padj%>%`<`(0.05)%>%table








################################################################################
########Let's do limma with stepwise factors
################################################################################
library(DESeq2)
library(limma)
colData(dds)
stepcoldata <- cbind(colData(dds),map(unique(dds[['nstage']]), ~ dds[['nstage']]>= .)%>%simplify2array%>%set_colnames(unique(dds[['type']])))

ribovoom <- voom(ribocounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=stepcoldata[isribo,]))
rnavoom <- voom(rnacounts[!eitherlow,],design = model.matrix(~ ns(nstage,3),data=stepcoldata[isribo,]))
allvoom <- voom(counts(dds)[!eitherlow,],design = model.matrix(~ assay*(ZYG+PACH+MII),data=stepcoldata[,]))












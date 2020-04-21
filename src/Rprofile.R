#options(repos = BiocManager::repositories());packrat::init(options = list(ignored.packages = c('sleuth','xtail','SaTAnn','RiboseQC','ORFquant','rseq','rseqdata','proDD','riboWaltz','colorout')))

try(silent=T,{library(colorout)})
library(Biostrings)
library(checkmate)
library(memoise)
library(assertthat)
library(stringr)
library(tidyverse)
library(magrittr)
library(checkmate)
# library(conflicted)
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(biomaRt))
suppressMessages(library(testthat))
library(zeallot)
library(splines)
library(GenomicRanges)
library(limma)
library(broom)
library(hashmap)




safe_hashmap<-setRefClass("Safe_Rcpp_Hashmap",
      contains="Rcpp_Hashmap",
      inheritPackage=TRUE
)
setMethod('[[','Safe_Rcpp_Hashmap',function (x, i, j, default,...){
  hashmapname = ''
    .local <- function (x, i, j, ..., exact = TRUE)
    {
        x$`[[`(i)
    }
    out <- .local(x, i, j, ...)
    if(missing(default)){
      if(any(is.na(out))){
        keymissingtxt = as.character(i[is.na(out)])%>%head%>%{ifelse(nchar(.)>15,paste0(substr(.,0,13),'...'),.)}%>%paste0(collapse='...')
        stop(paste0('Keys missing from safe hashmap: ',hashmapname,':',keymissingtxt))        
      }
    } else if(length(default)==1){
      out[is.na(out)] <- default
    }else{
      out[is.na(out)] <- default[is.na(out)]
    }
    return(out)
})

# #
# conflict_prefer('setdiff','BiocGenerics')
# conflict_prefer('rowMedians','Biobase')
# conflict_prefer('setequal','S4Vectors')
# conflict_prefer("between", "dplyr")
# conflict_prefer("intersect", "BiocGenerics")
# conflict_prefer("lag", "dplyr")
matches <- dplyr::matches
filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice
qs<-checkmate::qassert


# ###memoise
project_cache=here::here('R_cache')%T>%dir.create(showWarnings=F)
message('current cache size:')
message(system(str_interp('du -h ${project_cache}'),intern=T))
mycache=memoise::cache_filesystem(project_cache)
message(str_interp('Number of keys: ${length(mycache$keys())}'))

myclearcache=function() system(str_interp('rm -rf ${project_cache}'))


mymemoise <- function(f){
  if(!is.memoised(f)){
    memoise(f,cache=mycache)
    } else{ 
      f
  }
}
projmemoise<-mymemoise
addfileinf <- function(file){
  attr(file,'fileinfo')<-file.info(file)
  file
}

if(!interactive()) mymemoise=identity
  gigsused <- function(x)system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"),intern=TRUE)%>%str_extract('\\d+')%>%as.numeric%>%divide_by(1e6)
  message('memory in use ',gigsused())


# rm(foomat)
# gc(reset=TRUE,full=TRUE)
# message('memory in use ',gigsused())


# #
# foo<-function(x){message('foooooobar called')}
# myf <- mymemoise(foo)
# myf(2)

# Q
#  2 %>% mymemoise(function(x,.foo=foo){.foo();message('foooo'); x +1})(.)
# # e <- as.call(c(as.name("{"),quote(message('foo_ins')),body(foo)[-1]))
# body(foo) <- e

# foo

# mymemoise(function(x,.foo=foo)

safe_filter <- function(...){
  filtered = filter(...)
  assert_that(nrow(filtered)>0)
  filtered
}

#' safe_left_join
#' @description left join that fails if a row in x is either duplicated or
#'   unmatched.
#' @param x table to join
#' @param y table to join
#' @param by a character vector of column names to join by.
#' @param verbose Default is TRUE.
#' @export


safe_left_join = function (x, y, by = NULL, verbose = TRUE,allow_missing=FALSE,allow_dups=FALSE) {
  rows_start = nrow(x)

  if (is.null(by)) {
    by = intersect(names(x), names(y))
  } else {
    by = as.character(by)
  }

  y[["..1.."]] = 1
  x = left_join(x, y, by)

  if(!allow_dups){
    if (nrow(x) > rows_start) {
      stop("Rows have been duplicated in 'safe' left join")
    }
  }
  if(!allow_missing){
    if (any(ind <- is.na(x[["..1.."]]))) {
      sample = sample(which(ind), min(10, sum(ind)))
      examples = distinct(x[sample, by, drop = FALSE])
      if (verbose) print(examples)
      stop(sprintf("Failed to match %d rows in x.", sum(ind)))
    }
  }

  x[["..1.."]] = NULL

  x

}

is_offchr<-function(gr,si){
  if(is(gr,'GenomicRangesList')){
   (end(gr) > split(seqlengths(gr)[as.character(unlist(seqnames(gr)))],gr@partitioning) ) %in% TRUE
  }else{
    seqinfo(gr)<-si
    end(gr) > seqlengths(gr)[as.character(seqnames(gr))]

  }
}
is_out_of_bounds <- function(gr,si = seqinfo(gr)){
  start(gr)<1 | is_offchr(gr,si) 
}

get_all_obsizes <- function(){.GlobalEnv%>%names%>%discard(~is.function(get(.)))%>%setNames(.,.)%>%map(~object.size(get(.)))}

allobjsizes<-get_all_obsizes()
allobjsizes%<>%enframe
allobjsizes$value%<>%unlist
allobjsizes$value%<>%divide_by(1e6)
allobjsizes%>%arrange(desc(value))
# obsizes <- get_all_obsizes()

# obsizes%>%enframe('object','size')%>%mutate(size=as.numeric(size)/1e6)%>%arrange(desc(size))


purely <- function(fun,throw_error=TRUE,allow_functions=FALSE){
  
  funname <- rlang::quo_text(enquo(fun))

  globalobs <- .GlobalEnv%>%names

  globalobs <-   discard(globalobs,~is.function(get(.,envir=.GlobalEnv)) && (!identical(environment(get(.,envir=.GlobalEnv)),.GlobalEnv )))

  if(!allow_functions){
    globalobs <-   discard(globalobs,~is.function(get(.,envir=.GlobalEnv)) )
  }
 

  keep(globalobs,~is.function(get(.,envir=.GlobalEnv)) && (!identical(environment(get(.,envir=.GlobalEnv)),.GlobalEnv )))

  assert_that(! '%>%' %in% globalobs)

  funargs <- names(formals(fun))

  environment(fun) <- new.env()

  if(throw_error) messagefun = stop else messagefun = warning

  for(obnm in globalobs) {

    warnenv<-new.env()
    warnenv[['obnm']]<-obnm
    warnenv[['val']]<- force(get(obnm,envir=.GlobalEnv))
    delayedAssign(obnm,
      {
        messagefun(force(paste0('Object: ',obnm,' is being evaluated, but is not a formal argument of ',funname)));
        val
      },
      eval.env=warnenv,
      assign.env=environment(fun))
  }
  
  return(fun)
}

# myfun<-function(a){ a+1+b}

# a<-1
# b<-10
# myfun(a)
# pf<-purely(myfun)

# purely(myfun)(1)

# environment(pf)$get_frame_entropy



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}
load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


extract_oneof <- function(strings,ids){
  
  matchlist <- map(strings,~str_extract(pattern = ids,string = .))
  
  matchnum <- matchlist%>%map(~sum(!is.na(.)))
  
  stopifnot(all(matchnum < 2 ))
  stopifnot(all(matchnum > 0))

  matches <- matchlist%>%map_chr(keep,Negate(is.na))

  matches
}

extract_id <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = sampleids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}
load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


extract_oneof <- function(strings,ids){
  
  matchlist <- map(strings,~str_extract(pattern = ids,string = .))
  
  matchnum <- matchlist%>%map(~sum(!is.na(.)))
  
  stopifnot(all(matchnum < 2 ))
  stopifnot(all(matchnum > 0))

  matches <- matchlist%>%map_chr(keep,Negate(is.na))

  matches
}

DT2GR = function(dt,seqinf=si,checksi=TRUE){

  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(c('seqnames','start')%in%colnames(dt))
  stopifnot(c('width')%in%colnames(dt)|c('end')%in%colnames(dt))
  if(checksi){stopifnot(dt[['seqnames']] %in% seqlevels(seqinf))
  }else{seqinf=NULL}
  
  hasend=FALSE
  haswidth=FALSE

  if('end' %in% colnames(dt) ){
    stopifnot (dt[['end']] %>% `>`(0) %>%all)
    hasend=TRUE
  }
  if('width' %in% colnames(dt) ){
    stopifnot (dt[['width']] %>% `>`(0) %>%all)
    haswidth=TRUE
  }
  
  stopifnot(dt[['start']] %>% is.numeric)
  stopifnot(hasend|haswidth )
  
  if(haswidth & ! hasend ){
    dt[['end']]  = dt[['start']]+dt[['width']]-1 
  } 
  if(hasend ){

  } 

  #strand
  if(! 'strand' %in% colnames(dt)){
    dt[['strand']] ='*'
  }

  stopifnot(dt[['strand']] %in% c('+','-','*'))
  



  mdatcols = colnames(dt) %>% setdiff(c('seqnames','start','width','strand','end')) 
  #create basic granges object
  if(checksi){
    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']],seqinfo=seqinf)
  }else{    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])}

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

    stopifnot(all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')))

  gr
}


mergeGR2DT = function(mergegr){
  grcols = mergegr%>%vapply(is,TRUE,'GRanges')%>%which
  mergedt=cbind(
    mergegr[[grcols[[1]]]]%>%GR2DT,
    mergegr[[grcols[[2]]]]%>%GR2DT
  )
  cnames = colnames(mergedt)
  cnames[duplicated(cnames)]%<>%paste0('.1')
  setnames(mergedt,cnames)
  mergedt
}

GR2DT = function(gr){
  if(is.data.table(gr)) {
    warning('already a data table')
    return(gr)
  }
  #all columns must be vectors
  for(i in seq_len(ncol(mcols(gr)))){
    mcols(gr)[[i]] = as.vector(mcols(gr)[[i]])
  }

  dt = as.data.frame(gr,row.names=NULL,optional=FALSE)%>%as.data.table
  dt$strand= as.character(dt$strand)
  # setkey(dt,seqnames,strand,start,end)

  stopifnot( all(colnames( mcols(gr)) %in% colnames(dt) )  )
  dt
}



load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}


# reads_tr<-oldenv$cdsread_trmap
fp<-function(gr) ifelse(strand(gr)=='-',end(gr),start(gr))
tp<-function(gr) ifelse(strand(gr)=='-',start(gr),end(gr))
fpend<-function(x)resize(x,1,'start')
tpend<-function(x)resize(x,1,'end')

# gr1<-GRanges(c('chr1:5-6:+'))
# gr2<-GRanges(c('chr1:50-51:+','chr1:40-51:+'))

downstream_dist_till<-function(gr1,gr2){
  (fp(gr2)[precede(gr1,gr2)] - tp(gr1)) * ( ifelse(strand(gr1)=='+',1,-1))
}

upstream_dist_till<-function(gr1,gr2){
  (tp(gr2)[follow(gr1,gr2)] - fp(gr1)) * ( ifelse(strand(gr1)=='+',-1,1))
}

istpmost<-function(cds,groupvar='transcript_id'){ 
  ids <- seq_along(cds)
  tpmostids<-data_frame(id=ids,end=end(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(end*ifelse(strand=='-',-1,1)))%>%.$id
  ids %in% tpmostids
}
isfpmost<-function(cds,groupvar='transcript_id'){
  ids <- seq_along(cds)
  fpmostids<-data_frame(id=ids,start=start(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(start*ifelse(strand=='-',1,-1)))%>%.$id
  ids %in% fpmostids
}



clip_start <- function(x,n) resize(x,width(x)-n,fix='end')
clip_end <- function(x,n) resize(x,width(x)-n,fix='start')
setstrand<-function(x) {strand(x)<-Rle('+') 
  x}


# BiocManager::install('GenomicRanges')
testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:5:-','chr1:40-51:+'))),5)
testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:14:+','chr1:40-51:+'))),4)
testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:12:-','chr1:40-51:+'))),2)
testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:11:+','chr1:6:+','chr1:40-51:+'))),4)
testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:14:-','chr1:40-51:+'))),4)




take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
     if(length(x)<25){
          remain<-50-length(x)
          x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
     }

     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
     
     resSpec2<-dropFreqs(resSpec1,0.29,0.39)
     
     closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
     
     freq_max_3nt<-resSpec1$freq[closestfreqind]
     Fmax_3nt<-resSpec1$mtm$Ftest[closestfreqind]
     spect_3nt<-resSpec1$spec[closestfreqind]
     return(c(Fmax_3nt,spect_3nt))
     
}


ftestvect<-function(psit,k=24,bw=12){
  sl<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
  vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = sl)
  pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
  return(c(vals[2],pval))
}

#' Plot Heatmap
#'
#' Input matrix will be normalized (scaled and centered) by row. Then, values
#' smaller than z_min are set to z_min, likewise values larger than z_max are
#' set to z_max.
#'
#' @param mat a matrix of numbers
#' @param z_min all values in scaled matrix smaller than z_min are set to z_min
#' @param z_max all values in scaled matrix larger than z_max are set to z_max
#' @param dist_fun a distance function used for hierarchical clustering
#' @param title string, the title of the heatmap
#'
#' @return a list (see ?gplots::heatmap.2 for details)
#' @export
plot_heatmap <- function(mat, title = "", z_min = -Inf, z_max = Inf, dist_fun = NULL){
  assert_that(is.matrix(mat))
  assert_that(is.numeric(mat))
  scaled_mat <- t(scale(t(mat)))
  
  dual_scaled_mat <- pmin(pmax(scaled_mat, z_min), z_max)
  # pmin returns the minima of two vectors by position; the shorter vector gets recycled
  # pmax returns the maxima of two vectors by position; the shorter vector gets recycled
  if (is.null(dist_fun)) {
    cor_dist <- function(mat) {
      my_dist <- as.dist(1 - cor(t(mat), use = 'pairwise.complete.obs'))
      return(my_dist)
    }
    dist_fun <- cor_dist
  }else{
    assert_that(is.function(dist_fun))
  }
  
  gplots::heatmap.2(dual_scaled_mat[,10:1],
                    trace = 'none',
                    scale = 'none',
                    distfun  = cor_dist,
                    margins = c(8,4),
                    srtCol = 45,
                    labRow = F,
                    na.color = 'grey',
                    main = title,
                    col = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(128) ))
}

#' Plot Heatmap of Top Fluctuating Genes
#'
#' Takes a data frame with regularized log transformed counts, sorts genes by
#' decreasing standard deviation and includes the top number of fluctuating
#' features in a heatmap.
#'
#' @param num integer, number of top fluctuating features to be included in
#'   heatmap
#' @param rld_df a data frame with variables feature_id, sample_name,
#'   reg_log_count
#' @param ... parameters are passed to plot_heatmap
#'
#' @return a list (see ?gplots::heatmap.2 for details)
#' @export
plot_heatmap_fluc_features <- function(num, rld_df, ...){
  assert_that(is.count(num))
  assert_that(is.data.frame(rld_df))
  assert_that(all(has_name(rld_df, c("feature_id", "sample_name", "reg_log_count"))))
  
  rld_disp_wide <-
    rld_df %>%
    dplyr::select(feature_id, sample_name, reg_log_count) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::mutate(sd_reg_log_count = sd(reg_log_count)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(sample_name, reg_log_count) %>%
    dplyr::arrange(desc(sd_reg_log_count))
  
  rld_mat <-
    rld_disp_wide %>%
    dplyr::slice(1:num) %>%
    dplyr::select(-sd_reg_log_count) %>%
    as.data.frame(stringsAsFactors = F) %>%
    tibble::column_to_rownames("feature_id")  %>%
    as.matrix()
  
  plot_heatmap(rld_mat, paste0('top ', num, ' fluctuating features'), ...)

}

library(tidyverse)
library(GenomicFeatures)

grl = GRangesList(list(c(
  GRanges('a:3-6:+',foo=1),
  GRanges('a:8-10:+'),
  GRanges('a:13-15:+')
),  c(GRanges('a:3-6:-'),
      GRanges('a:8-10:-'),
      GRanges('a:13-15:-',bar=2)
)))%>%setNames(letters[1:2])


resize_grl_startfix<-function(grl,width){
  #what follows is some slightly black magic using S4 vectors
  #Integerlist which showings how much we'd need to trim that exon to get to to the desired transcript length
  trim =  cumsum(width(grl)) - width 
  #Where trim is greater than the exon width, we drop it
  drop = trim >=  width(grl)
  grl = grl[!drop]
  #vector showing location of the new 3' end of each transcript
  newends = cumsum(elementNROWS(grl))
  #vector with the amount we need to trim each new 3' end by
  endtrims=trim[IntegerList(as.list(elementNROWS(grl)))]@unlistData
  #finally, use these to trim
  grl@unlistData[newends] <- resize(grl@unlistData[newends], width(grl@unlistData[newends]) - endtrims  )
  grl
  
}

str_order_grl<-function(grl){order( start(grl)*(((strand(grl)!='-')+1)*2 -3) )}
sort_grl_st <- function(grl)grl[str_order_grl(grl),]
resize_grl_endfix <- function(grl,width){
  grl = invertStrand(grl)%>%sort_grl_st
  
  grl = resize_grl_startfix(grl,width)
  invertStrand(grl)%>%sort_grl_st
}
resize_grl <- function(grl,width,fix='start',check=TRUE){
  stopifnot(all(width>0))
  assert_that(all(all(diff(str_order_grl(grl))==1) ),msg = "grl needs to be 5'->3' sorted")
  if(fix=='start'){
    grl = resize_grl_startfix(grl,width)
  }else if(fix=='end'){
    grl = resize_grl_endfix(grl,width)
  }else if(fix=='center'){
    grlwidths = sum(width(grl)) 
    diffs = (width - grlwidths)
    
    grl = resize_grl_startfix(grl,grlwidths - floor(diffs/2))
    grl = resize_grl_endfix(grl,grlwidths - ceiling(diffs/2))
    
  }
  if(check){
    startstoolow <- any(start(grl)<=0)
    if(any(startstoolow)){
      stop(str_interp("${sum(startstoolow)} ranges extended below 1 .. e.g. ${head(which(startstoolow,1))}"))
    }
    grlseqs <- as.vector(unlist(use.names=F,seqnames(grl)[IntegerList(as.list(rep(1,length(grl))))]))
    endstoohigh <- any((end(grl)>seqlengths(grl)[grlseqs])%in%TRUE)
    if(any(endstoohigh)){
      stop(str_interp("${sum(endstoohigh)} ranges extended below above seqlength .. e.g. ${head(which(endstoohigh,1))}"))
    }
  }
  grl
}
trim_grl <- function(grl,bp,end='tp'){
  if(end=='tp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='start')
  }else if(end=='fp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='end')
  }else {
    stop("end should be fp or tp")
  }
}

setGeneric('apply_psite_offset',function(offsetreads,offset) strandshift(offsetreads,offset))
setMethod('apply_psite_offset','GAlignments',function(offsetreads,offset){
  if(is.character(offset)){
    offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
  } 
  if(length(offset)==1) offset = rep(offset,length(offsetreads))
  isneg <-  as.logical(strand(offsetreads)=='-')
  offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
  ends <- qwidth(offsetreads[isneg])-offset[isneg]
  offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 
  offsetreads
})

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift) shift(gr , ifelse( strand(gr)=='-',- shift,shift))

get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  #
  if(is.null(offsets)){
  mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    mcols(reads)$offset <- 
      data.frame(length=qwidth(reads),compartment='nucl')%>%
      safe_left_join(offsets,allow_missing=TRUE)%>%.$offset
  }
  #
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  # 
  mcols(reads)$length <- width(reads)
  reads%<>%subset(!is.na(offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  mcols(psites)$length <- mcols(reads)$length   
  psites
}
#exonsexp->exons_grl
# cds[bestcds]->trspacegr
spl_mapFromTranscripts <- function(trspacegr,exons_grl){

  exons_tr<-exons_grl%>%unlist%>%mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  # grlcs = exons_grl%>%sort_grl_st%>%width%>%cumsum%>%IntegerList
  
  # grlcs%>%head%>%{pc(.,rep(1,length(.)))}
  
  # grlcs%>%{.[IntegerList(as.list(-elementNROWS(.)))]}

  # lastelems <- elementNROWS(grlcs)%>%as.list%>%IntegerList%>%setNames(NULL)
  # pc(rep(1,length(grlcs)),grlcs+1)[-1 * lastelems]


  # exons_grl%>%unlist%>%.[exons_tr[subjectHits(ov)]$xHits]
  ov <- findOverlaps(trspacegr,exons_tr)
  # trspacegr[testpid]

  # exons_grl[seqnames(trspacegr)]

  # exons_grl%>%unlist%>%setNames(paste0('exon_',seq_along(.)))%>%


  # mergeByOverlaps(trspacegr,exons_tr)

  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
  trspacegr_spl,
  # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
  exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}

#metaplots




# 
# resize_grl(grl,sum(width(grl))+1,fix='end')
# 
# cdsgrl<-cds%>%split(.,.$protein_id)
# cdsgrl%>%width%>%sum%>%min
# cdsgrl <- cdsgrl%>%sort
# 
# out<-resize_grl(cdsgrl,3)
# 
# 
# cdsgrl[1]%>%trim_grl(-2)

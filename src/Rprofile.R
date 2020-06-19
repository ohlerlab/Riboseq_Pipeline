#options(repos = BiocManager::repositories());packrat::init(options = list(ignored.packages = c('sleuth','xtail','SaTAnn','RiboseQC','ORFquant','rseq','rseqdata','proDD','riboWaltz','colorout')))

try(silent=T,{library(colorout)})
message('loading libraries')
# library(Biostrings)
# library(checkmate)
# library(stringr)
# library(conflicted)
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(biomaRt))
suppressMessages(library(testthat))
suppressMessages(library(memoise))
suppressMessages(library(checkmate))
suppressMessages(library(GenomicRanges))
suppressMessages(library(hashmap)) ##has to be installed like this devtools::install_github("nathan-russell/hashmap")
suppressMessages(library(zeallot))
suppressMessages(library(ggpubr))
suppressMessages(library(gplots))
suppressMessages(library(xtail)) ## devtools::install_github('xryanglab/xtail')
# 
# library(zeallot)
# # library(splines)
# library(GenomicRanges)
# library(limma)
# library(broom)
# library(hashmap)

matches <- dplyr::matches
filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice
qs<-checkmate::qassert

get_all_obsizes <- function(){.GlobalEnv%>%names%>%discard(~is.function(get(.)))%>%setNames(.,.)%>%map(~object.size(get(.)))}


#make an object like a fast, vectorized python dictionary, that FAILS if you look up something
#that isn't in it.
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




# ###memoise
project_cache=here::here('R_cache')
if(!exists('project_cache'))project_cache=tempdir()
message(system(str_interp('du -h ${project_cache}'),intern=T))
mycache=memoise::cache_filesystem(project_cache)

myclearcache=function() system(str_interp('rm -rf ${project_cache}'))

#memoise to the folder
projmemoise<- function(f){
  if(!is.memoised(f)){
    memoise(f,cache=mycache)
    } else{ 
      f
  }
}
#for file operations, memoises not just on the name but on the mod date etc.
addfileinf <- function(file){
  attr(file,'fileinfo')<-file.info(file)
  file
}

if(!interactive()) mymemoise=identity
  gigsused <- function(x)system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"),intern=TRUE)%>%str_extract('\\d+')%>%as.numeric%>%divide_by(1e6)
  message('memory in use ',gigsused())

  #' safe_left_join
  #' @description filter that fails if it returns no rows
  #'   unmatched.
  #' @param .data the table
  #' @param ... one or more expressions to be passed on to filter as conditions
  #' @export
  #' 
  safe_filter <- function(.data,...){
    filt_df = filter(.data, ...)
    assertthat::assert_that(nrow(filt_df)>0,msg = str_interp("filter of ${nrow(.data)} row data frame results in zero rows"))
    filt_df
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

#for checking granges objects are on their proper chromosome
is_offchr<-function(gr,si){
  seqinfo(gr)<-si
  end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
}
is_out_of_bounds <- function(gr,si = seqinfo(gr)){
  start(gr)<1 | is_offchr(gr,si) 
}

# allobjsizes<-get_all_obsizes()
# allobjsizes%<>%enframe
# allobjsizes$value%<>%unlist
# allobjsizes$value%<>%divide_by(1e6)
# allobjsizes%>%arrange(desc(value))
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

load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


#useful funtion for converting back and forth between GRanges objects
#and data tables
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

#Granges utilities functions
fp<-function(gr) ifelse(strand(gr)=='-',end(gr),start(gr))
tp<-function(gr) ifelse(strand(gr)=='-',start(gr),end(gr))
fpend<-function(x)resize(x,1,'start')
tpend<-function(x)resize(x,1,'end')
setstrand<-function(x) {strand(x)<-Rle('+'); x}

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
  grl = invertStrand(grl)
  grl = invertStrand(grl)%>%sort_grl_st
  
  grl = resize_grl_startfix(grl,width)
  invertStrand(grl)%>%sort_grl_st
}
resize_grl <- function(grl,width,fix='start',check=TRUE){
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
  
  gplots::heatmap.2(dual_scaled_mat[,],
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

### try to flirt with the DESeq2 functions
plotPCA_allPC = function(object, intgroup="condition", ntop=500, returnData=FALSE,XPC='PC1',YPC = 'PC2')
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(pca$x, group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x=XPC, y=YPC, color="group")) + geom_point(size=3) + 
    xlab(paste0(XPC,': ',round(percentVar[as.integer(str_sub(XPC,-1))] * 100),"% variance")) +
    ylab(paste0(YPC,': ',round(percentVar[as.integer(str_sub(YPC,-1))] * 100),"% variance")) +
    coord_fixed()
}
# 
# grl = GRangesList(list(c(
#   GRanges('a:3-6:+',foo=1),
#   GRanges('a:8-10:+'),
#   GRanges('a:13-15:+')
# ),  c(GRanges('a:3-6:-'),
#       GRanges('a:8-10:-'),
#       GRanges('a:13-15:-',bar=2)
# )))%>%setNames(letters[1:2])




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

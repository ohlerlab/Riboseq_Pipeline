library(here)
library(data.table)
library(magrittr)
library(tidyverse)
library(GenomicFeatures)
library(dplyr)
library(assertthat)
select<-dplyr::select


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


inclusiontable<-function(a,b){
  library(rlang)
  aname = rlang::quo_text(enquo(a))
  bname = rlang::quo_text(enquo(b))
  all = BiocGenerics::union(a,b)
  outtab = table(all %in% a,all %in% b)
  dimnames(outtab)%<>%setNames(c(aname,bname))
  outtab
}

make_quantcompplot <- function(compdf, col1, col2, fname){
  require(LSD)
  base::source(here('Applications/LSD/R/LSD.heatscatter.R'))
  require(broom)
  col1<-enquo(col1)
  col2<-enquo(col2)
  corlabel = compdf%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
    summarise(tidy(cor.test(!!col1, !!col2)))
  corlabel = corlabel%>%
    mutate(
      pformat=format(p.value,format='e',digits=4),
      pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
      labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
  #
  nlabel=tibble(labl=paste0('N=',nrow(compdf)))
  pdf(fname)
  gplot = heatscatter(ggplot=TRUE,
      compdf[[quo_name(col1)]],compdf[[quo_name(col2)]])+
    scale_x_continuous(quo_name(col1))+
    scale_y_continuous(quo_name(col2))+
    ggtitle(basename(fname))+
    geom_text(show.legend=F,data=corlabel,
      hjust=1,vjust=1,x= Inf,y=Inf,aes(label=labl))+
    geom_text(show.legend=F,data=nlabel,
      hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))
  dev.off()
  pdf(fname)
  gplot
  print(gplot)
  dev.off()
  message(normalizePath(fname))
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


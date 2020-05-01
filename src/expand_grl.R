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
  grl@unlistData[newends] %<>% resize(grl@unlistData[newends], width(grl@unlistData[newends]) - endtrims  )
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


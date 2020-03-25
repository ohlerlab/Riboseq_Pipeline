library(tidyverse)
library(data.table)
library(magrittr)

#read the index file
indexfasta <- 'pipeline/tRNA_rRNA_index/tRNA_rRNA_index.fa'
contam_headers <- readLines(indexfasta) %>% str_extract('>.*')
contam_headers_df <- data.frame(header = contam_headers)
contam_headers_df %<>% filter(!is.na(header))
#extract unique ids from the contamination fasts
contam_headers_df %<>%mutate(contam_id = str_extract(header,'(?<=>)\\d+') %>% as.numeric)	
contam_headers_df %>% head
#all ids extracted successfully
stopifnot((contam_headers_df[['contam_id']]==(1:nrow(contam_headers_df))))

#read in our idx files
idxfiles <- Sys.glob('pipeline/filter_reads/*/*.idxtmp')

idxfiles%<>%setNames(.,basename(dirname(.)))
idxfile <- idxfiles[1]

samples_idx <- idxfiles%>%map_df(.id='sample',function(idxfile){
	sample_idx <- fread(idxfile) %>%
		arrange(- V3)%>%		set_colnames(c('seq','length','reads','unmapped'))
	#extract our id from the headers
	sample_idx%<>%filter(seq!='*')
	sample_idx%<>%mutate(contam_id = str_extract(seq,'^\\d+')%>%as.numeric)
	#all ids extracted successfully
	stopifnot(sample_idx$contam_id%>%is.na%>%not%>%all)
	#now join headers
	sample_idx%<>%left_join(contam_headers_df,by='contam_id')
	#joined to headers successfully
	stopifnot(sample_idx$header%>%is.na%>%not%>%all)
	#sam as total of read colum in idx, actually
	# sample_idx%<>%mutate(frac = reads / alignlogvect['total'])
	#now categorize
	sample_idx <- sample_idx%>%
		mutate(mt_rRNA =header%>%str_detect(regex(ignore_case=T,'chrM.*MT-RNR')))%>%
		mutate(rRNA =header%>%str_detect(regex(ignore_case=T,'rRNA|(ribosomal RNA)')))%>%
		mutate(tRNA =header%>%str_detect(regex(ignore_case=T,'(transfer RNA)|tRNA')))%>%
		mutate(miRNA =header%>%str_detect(regex(ignore_case=T,'_MIR\\d+\\\\w|_MIR\\d+\\-\\d$|_MIR\\d+$|microRNA|miR\\-|MIRLET\\w\\dMIR\\d+|miRt\\d\\w')))%>%
		mutate(snoRNA =header%>%str_detect(regex(ignore_case=T,'(small nucleolar)|snoRNA|SNHG\\d+|_snR\\d{0,2}|SNOR|snoU\\d+|snoR\\d{1,4}|(Zm|_|HvU)U\\d{1,2}|MTU\\d{1,2}'		)))%>%
		mutate(snRNA =header%>%str_detect(regex(ignore_case=T,'(small nuclear)|LINC00910|RNU\\d+$|RNU\\d(\\-\\d)')))%>%
		mutate(SCARNA =header%>%str_detect(regex(ignore_case=T,'SCARNA\\d+$')))%>%
		mutate(RPI_barcode =header%>%str_detect(regex(ignore_case=T,'RPI\\d+$')))%>%
		mutate(spikein =header%>%str_detect(regex(ignore_case=T,'ERCC')))%>%
		mutate(marker =header%>%str_detect(regex(ignore_case=T,'marker')))
		
	contamtypecols = sample_idx%>%select(-(seq:header))%>%colnames
	
	noncat=filter_at(sample_idx,vars(all_of(contamtypecols)),all_vars(!.))

    noncat%>%select(-one_of(contamtypecols))%>%arrange(desc(reads))%>%head

	#test most things have a category above
	#
	uncatfrac <- sum(noncat$reads)/sum(sample_idx$reads)
	if(uncatfrac > 0.05){
		message(capture.output(noncat%>%arrange(-reads)%>%select(-one_of(contamtypecols))%>%head)%>%paste0(col='\n'))
		stop('more than 5% of contaminants fall in the unkown category')
	}
	if(sample_idx[,c('tRNA','rRNA','miRNA','snoRNA','spikein')]%>%rowSums%>%`<=`(1)%>%all%>%not){
		stop('some of these are categorized in two bins, you need to adjust')
	}
	sample_idx$category <- sample_idx[,contamtypecols]%>%mutate(uncategorized=TRUE)%>%apply(1,which)%>%map_dbl(head,1)%>%{c(contamtypecols,'uncategorized')[.]}

	# sample_idx <- sample_idx%>%
	# 	mutate(category = case_when(
	# 		rRNA ~ 'rRNA',
	# 		tRNA ~ 'tRNA',
	# 		miRNA ~ 'miRNA',
	# 		snoRNA ~ 'snoRNA',
	# 		spikein ~ 'spikein',
	# 		TRUE ~ 'uncategorized'
	# 	))
	sample_idx
})

samples_alignlogstats <- idxfiles%>%map_df(.id='sample',function(idxfile){
	#fraction of reads

	repfile <- dirname(idxfile)%>%list.files(full=T,patt='Log.final.out')

	replines<-repfile%>%readLines%>%.[-length(.)]

	alignlogvect <- c(
		total=replines%>%str_subset('Number of input reads')%>%str_extract('\\d+')%>%as.numeric,
		noncontam = replines%>%str_subset('Number of reads unmapped')%>%str_extract('\\d+')%>%as.numeric%>%sum,
		one_contam = replines%>%str_subset('Uniquely mapped reads number')%>%str_extract('\\d+')%>%as.numeric%>%sum,
		mult_contam = replines%>%str_subset('Number of reads mapped to (multiple|too many)')%>%str_extract('\\d+')%>%as.numeric%>%sum
	)
	# alignlogvect/1e6
	# #make sure most or none align
	# stopifnot(
	# 	between(alignlogvect['unpaired']/alignlogvect['total'],1-0.05,1) |
	# 	between(alignlogvect['unpaired']/alignlogvect['total'],0,.05 ) 
	# )
	alignlogvect%>%enframe
})

#get mapped, noncontam reads for each sample
mapped_reads_df <- Sys.glob('pipeline/star/reports/*/*bamstats*')%>%
	setNames(.,basename(dirname(.)))%>%
	map(.%>%readLines%>%head(20)%>%str_extract('(?<=reads mapped:\t)\\d+')%>%.[!is.na(.)]%>%as.numeric)%>%enframe('sample','reads')%>%
	unnest(reads)%>%
	select(sample,reads)%>%
	mutate(category='noncontam_mapped')

samples_idx$category%>%table

samples_idx_map <- samples_idx%>%bind_rows(mapped_reads_df)

totalsdf <- samples_alignlogstats%>%filter(name=='total')%>%select(sample,total=value)
noncontamsdf <- samples_alignlogstats%>%filter(name=='noncontam')%>%select(sample,noncontam=value)

sampls_catsum_df<-samples_idx_map%>%
	group_by(sample,category)%>%
	summarise(reads=sum(reads))%>%
	arrange(-reads)%>%
	left_join(totalsdf)

sampls_catsum_df%<>%mutate(totalfrac = reads / total)

sampls_catsum_df$category%>%table
sampls_catsum_df_wide <- sampls_catsum_df%>%
	select(category,totalfrac,sample)%>%
	spread(category,totalfrac)%>%as.data.frame%>%
	left_join(totalsdf)%>%
	left_join(noncontamsdf)%>%
	mutate(noncontam = noncontam/total)%>%
	select(total,noncontam,noncontam_mapped,rRNA,tRNA,everything())

sampls_catsum_df_wide


sampls_catsum_df_wide$noncontam_mapped/sampls_catsum_df_wide$noncontam

abs_sampls_catsum_df_wide <- sampls_catsum_df%>%select(category,reads,sample)%>%spread(category,reads)%>%as.data.frame
abs_sampls_catsum_df_wide%<>%left_join(sampls_catsum_df%>%distinct(sample,total))

dir.create(here::here('tables'))
sampls_catsum_df_wide%>%mutate_if(.predicate=is.numeric,list(~ round(.,3)))%>%write_tsv('tables/sampls_catsum_df_wide.tsv')
abs_sampls_catsum_df_wide%>%write_tsv('tables/abs_sampls_catsum_df_wide.tsv')


normalizePath('tables/sampls_catsum_df_wide.tsv')
normalizePath('tables/abs_sampls_catsum_df_wide.tsv')

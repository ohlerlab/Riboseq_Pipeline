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
		mutate(snoRNA =header%>%str_detect(regex(ignore_case=T,'(small nucleolar)|snoRNA|_snR\\d{0,2}|SNOR|snoU\\d+|snoR\\d{1,4}|(Zm|_|HvU)U\\d{1,2}|MTU\\d{1,2}'		)))%>%
		mutate(snRNA =header%>%str_detect(regex(ignore_case=T,'(small nuclear)|LINC00910')))%>%
		mutate(SCARNA =header%>%str_detect(regex(ignore_case=T,'SCARNA\\d+$')))%>%
		mutate(RPI_barcode =header%>%str_detect(regex(ignore_case=T,'RPI\\d+$')))%>%
		mutate(spikein =header%>%str_detect(regex(ignore_case=T,'ERCC')))
	contamtypecols = sample_idx%>%select(-(seq:header))%>%colnames
	
	noncat=filter_at(sample_idx,vars(all_of(contamtypecols)),all_vars(!.))

    noncat%>%select(-one_of(contamtypecols))%>%arrange(desc(reads))%>%head

	#test most things have a category above
	#
	uncatfrac <- sum(noncat$reads)/sum(sample_idx$reads)
	if(uncatfrac > 0.05){
		stop('more than 5% of contaminants fall in the unkown category')
	}
	if(sample_idx[,c('tRNA','rRNA','miRNA','snoRNA','spikein')]%>%rowSums%>%`<=`(1)%>%all%>%not){
		stop('some of these are categorized in two bins, you need to adjust')
	}
	#
	sample_idx <- sample_idx%>%
		mutate(category = case_when(
			rRNA ~ 'rRNA',
			tRNA ~ 'tRNA',
			miRNA ~ 'miRNA',
			snoRNA ~ 'snoRNA',
			spikein ~ 'spikein',
			TRUE ~ 'uncategorized'
		))
	sample_idx
})

samples_alignlogstats <- idxfiles%>%map_df(.id='sample',function(idxfile){
	#fraction of reads
	repfile <- idxfile%>%str_replace('.idxtmp','.alignreport.log')
	replines<-repfile%>%readLines%>%.[-length(.)]
	alignlogvect<-replines%>%str_extract('\\d+')%>%setNames(c('total','unpaired','noncontam','one_contam','mult_contam'))%>%sapply(as.numeric)#%>%enframe('stat','reads')
	alignlogvect/1e6
	#make sure most or none align
	stopifnot(
		between(alignlogvect['unpaired']/alignlogvect['total'],1-0.05,1) |
		between(alignlogvect['unpaired']/alignlogvect['total'],0,.05 ) 
	)
	alignlogvect%>%enframe
})

#get mapped reads for each sample
mapped_reads_df <- Sys.glob('pipeline/star/reports/*/*bamstats*')%>%setNames(.,basename(dirname(.)))%>%map(.%>%readLines%>%head(20)%>%str_extract('(?<=reads mapped:\t)\\d+')%>%.[!is.na(.)]%>%as.numeric)%>%enframe('sample','reads')%>%unnest%>%
	mutate(category='mapped')
samples_idx_map <- samples_idx%>%bind_rows(mapped_reads_df)

sampls_catsum_df<-samples_idx_map%>%group_by(sample,category)%>%summarise(reads=sum(reads))%>%arrange(-reads)%>%
	left_join(samples_alignlogstats%>%filter(name=='total')%>%select(sample,total=value))
sampls_catsum_df%<>%mutate(totalfrac = reads / total)

sampls_catsum_df_wide <- sampls_catsum_df%>%select(category,totalfrac,sample)%>%spread(category,totalfrac)%>%as.data.frame
sampls_catsum_df_wide$noncontam_unmapped <- sampls_catsum_df_wide[,-1]%>%rowSums%>%{1-.}

abs_sampls_catsum_df_wide <- sampls_catsum_df%>%select(category,reads,sample)%>%spread(category,reads)%>%as.data.frame
abs_sampls_catsum_df_wide%<>%left_join(sampls_catsum_df%>%distinct(sample,total))

dir.create(here::here('tables'))
sampls_catsum_df_wide%>%write_tsv('tables/sampls_catsum_df_wide.tsv')
abs_sampls_catsum_df_wide%>%write_tsv('tables/abs_sampls_catsum_df_wide.tsv')


normalizePath('tables/sampls_catsum_df_wide.tsv')

# if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
source(here::here("src/Rprofile.R"))
source(here::here("src/functions.R"))
library(tximport)

ribostanfiles <- Sys.glob("ribostan/data/rna*/quant.sf")
sampleinfo <- read_csv("../src/sample_config.csv")
ribosamples <- sampleinfo %>%
    filter(isriboseq) %>%
    .$sample_id
rnasamples <- sampleinfo %>%
    filter(!isriboseq) %>%
    .$sample_id
salmonfiles <- paste0("salmon/data/", rnasamples, "/quant.sf")
ribostanfiles <- paste0("ribostan/", ribosamples, "/", ribosamples, ".ribostan.tsv")
stopifnot(all(file.exists(salmonfiles)))
stopifnot(all(file.exists(ribostanfiles)))

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE, defaults = list(
    gtf = "pipeline/gencode.v37.primary_assembly.annotation.gtf",
    fafile = "pipeline/GRCh38.primary_assembly.genome.fa"
))

# Turn arguments into R variables
keys <- attachLocally(args)
cat("Command-line arguments attached to global environment:\n")
print(keys)
str(mget(keys, envir = globalenv()))
# }
#
fafile <- Rsamtools::FaFile(fafile)

gtf_gr <- rtracklayer::import(gtf)

#
countdatafiles <- c(salmonfiles, ribostanfiles)

################################################################################
######## Now load annotation data
################################################################################
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
tx2genetbl <- gtf_gr %>%
    mcols() %>%
    as.data.frame() %>%
    filter(!is.na(gene_id), !is.na(transcript_id)) %>%
    distinct(transcript_id, gene_id) %>%
    select(transcript_id, gene_id)
trid2gid <- tx2genetbl %>%
    {
        setNames(.$gene_id, .$transcript_id)
    }
salmontrs <- salmonfiles[1] %>%
    read_tsv() %>%
    .$Name %>%
    str_extract("[^|]+")
salmoncds <- salmontrs %>% intersect(gtf_gr %>% subset(type == "CDS") %>% .$transcript_id)
# cdswidth = gtf_gr%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%width%>%sum
# cdsrangedf <- tibble(Name=salmoncds,annolength = cdswidth[salmoncds])
# cdsrangedf$Name%<>%trimids
# trs = cdsrangedf$Name
trs <- salmontrs

# we have to trim the ids for orfquant
tx2genetbl$transcript_id %<>% trimids
tx2genetbl$gene_id %<>% trimids


################################################################################
######## Collect transript level info for top transcripts
################################################################################
tximport(ribostanfiles)
filereadfunc <- function(file) {
    if (file %>% str_detect("ribotrans")) {
        message(str_interp("reading riboemfile ${file}"))
        ribotransexprtbl <- read_tsv(file)
        outtable <- tibble(
            Name = ribotransexprtbl$tr_id,
            Length = ribotransexprtbl$cds_len,
            EffectiveLength = ribotransexprtbl$cds_len,
            NumReads = ribotransexprtbl$read_count
        ) %>% mutate(
            TPM = (NumReads / EffectiveLength) %>%
                {
                    1e6 * . / sum(.)
                }
        )
    } else {
        outtable <- file %>%
            read_tsv()
        outtable <- outtable %>% mutate(Name = str_extract(Name, "[^|]+"))
    }
    outtable$Name %<>% trimids
    full_outtable <- cdsrangedf %>% safe_left_join(allow_missing = TRUE, outtable)
    # stopifnot(full_outtable$NumReads%>%is.na%>%mean%>%`>`(0.8))
    full_outtable$NumReads %<>% replace_na(0)
    full_outtable$TPM %<>% replace_na(0)
    full_outtable$Length %<>% {
        .[is.na(.)] <- full_outtable$annolength[is.na(.)]
        .
    }
    full_outtable$EffectiveLength %<>% {
        .[is.na(.)] <- full_outtable$Length[is.na(.)]
        .
    }
    full_outtable$NumReads %>%
        sum(na.rm = T) %>%
        `>`(0)
    full_outtable %>% select(Name, Length, EffectiveLength, TPM, NumReads)
}

tx_countdata <- tximport(
    files = countdatafiles,
    ignoreTxVersion = TRUE,
    tx2gene = tx2genetbl,
    type = "salmon",
    countsFromAbundance = "scaledTPM",
    ignoreAfterBar = TRUE,
    importer = filereadfunc
)

tx_countdata$counts %>%
    apply(2, sum, na.rm = T) %>%
    divide_by(1e6)

iso_tx_countdata <- tximport(
    files = countdatafiles,
    txOut = TRUE,
    ignoreTxVersion = TRUE,
    tx2gene = tx2genetbl,
    type = "salmon",
    countsFromAbundance = "scaledTPM",
    importer = filereadfunc
)

stopifnot(iso_tx_countdata$abundance %>% rownames() %>% setequal(trs))

tx_countdata %>% saveRDS("r_data/tx_countdata.rds")
iso_tx_countdata %>% saveRDS("r_data/iso_tx_countdata.rds")
iso_tx_countdata %>% saveRDS("r_data/gtf_gr.rds")
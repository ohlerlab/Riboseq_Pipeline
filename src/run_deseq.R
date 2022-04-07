library(tximport)
library(DESeq2)
library(assertthat)

source(here::here("src/Rprofile.R"))
source(here::here("src/functions.R"))

if (!exists("gtf_gr")) gtf_gr <- readRDS(here("pipeline/r_data/gtf_gr.rds"))
if (!exists("tx_countdata")) tx_countdata <- readRDS(here("pipeline/r_data/tx_countdata.rds"))

# table with gene 2 gname transformation
gid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(gene_id, gene_name) %>%
  select(gene_id, gene_name)

# read in the sample_parameter file from the piepline
sampdf <- here("config/sample_config.csv") %>%
  read_csv() %>%
  select(-one_of(dropcols))
stopifnot(colnames(tx_countdata$counts) %in% sampdf$sample_id)

# use the factor ordering in the sample_parameters file - by default it gets
sampdf <- sampdf %>%
  as.data.frame() %>%
  set_rownames(.$sample_id) %>%
  mutate(group <- sample_id %>%
    str_replace("(_)?\\d+$", "") %>%
    str_replace("\\+", "pls") %>%
    str_replace("\\-", "neg")) %>%
  mutate_if(is.character, as_factor)
if (is.null(sampdf$sample_name)) sampdf$sample_name <- sampdf$sample_id

for (i in seq_len(ncol(sampdf))) {
  colnm <- colnames(sampdf)[[i]]
  if (colnm %in% c("sample_id", "sample_name")) next
  if (is.factor(sampdf[[i]])) {
    message(paste0("factor levels for ", colnm, ":"))
    message(paste(collapse = ",", levels(sampdf[[i]])))
  }
}


# read the counts data - this takes a while sometimes
# file.remove(here('pipeline/r_data/dds.rds'))
if (!file.exists(here("pipeline/r_data/dds.rds"))) {
  # fractional counts sometimes result from e.g. ribomap - fix this
  randomround <- function(x) floor(x) + rbinom(length(x), 1, x %% 1)
  fractionalcounts <- (tx_countdata$counts %% 1 == 0) %>% any()
  if (fractionalcounts) {
    for (i in seq_len(ncol(tx_countdata$counts))) {
      tx_countdata$counts[, i] %<>% randomround
    }
  }
  samples <- colnames(tx_countdata$counts[highcountgenes, ])
  # create a granges object with the coordinates in the first fcounts file
  dds <- DESeqDataSetFromTximport(tx_countdata,
    colData = sampdf[match(samples, sampdf$sample_id), ],
    design = ~group
  )
  rownames(dds) <- tx_countdata$counts %>% rownames()
  rowRanges(dds) <- gtf_gr %>%
    subset(type == "exon") %>%
    split(., trimids(.$gene_id)) %>%
    .[rownames(dds)]
  rowData(dds) <- gtf_gr %>%
    subset(type == "gene") %>%
    split(., trimids(.$gene_id)) %>%
    .[rownames(dds)] %>%
    unlist() %>%
    mcols()
  # save deseq object
  rownames(dds) <- rowData(dds)$gene_id
  rownames(colData(dds)) <- colData(dds)$sample_id
  saveRDS(dds, here("pipeline/r_data/dds.rds"))
} else {
  dds <- readRDS(here("pipeline/r_data/dds.rds"))
}
cat("successfully loaded  countdata and annotation")
stopifnot(exists("sampdf"))


# you might want to use vst if there's a lot of data, for speed
dir.create(here("data"))
# file.remove('pipeline/r_data/normcounts.rds')
if (!file.exists(here("pipeline/r_data/normcounts.rds"))) {
  normfunc <- if (ncol(dds) > 20) DESeq2::vst else DESeq2::rlog

  normcounts <- normfunc(dds)
  rownames(normcounts) <- rownames(dds)
  saveRDS(normcounts, here("pipeline/r_data/normcounts.rds"))
} else {
  normcounts <- readRDS(here("pipeline/r_data/normcounts.rds"))
}
cat("successfully normalized  countdata")

stopifnot(exists("dds"))
stopifnot(exists("normcounts"))
stopifnot(exists("sampdf"))
{
  # do a trial run of DESeq to get the contrast names
  design <- as.formula("~ assay * (subunit+induced:subunit)")
  design(dds) <- design
  testdds <- DESeq(head(dds, 2 * ncol(dds)))
  resnames <- resultsNames(testdds)
  resnames
}


################################################################################
######## Construct Model
################################################################################

resnames <- resultsNames(testdds)

contrasts <- list(
  "assay_ribo_vs_total",
  "assayribo.subunit4E.inducedyes",
  "assayribo.subunit4G1.inducedyes",
  "assayribo.subunit4G2.inducedyes",
  "assayribo.subunit4G3.inducedyes"
) %>% setNames(., .)

# check the contrasts
for (contrastname in names(contrasts)) {
  contrast <- contrasts[[contrastname]]
  if (is.character(contrast) & length(contrast) == 1) {
    assert_that(contrast %in% resnames)
    contrast <- as.numeric(resnames == contrast)
    contrasts[[contrastname]] <- contrast
  }
  message(contrastname)
  try({
    results(testdds, contrast)
  })
}

cat("tested contrasts")



# you might want to use vst if there's a lot of data, for speed
dir.create(here("data"))
# file.remove('pipeline/r_data/resultslist.rds')
if (!file.exists(here("pipeline/r_data/resultslist.rds"))) {
  design(dds) <- design

  dds <- DESeq(dds, betaPrior = F)
  # and get the contrasts we need
  resultslist <- lapply(contrasts, results, object = dds)
  names(resultslist) <- names(contrasts)
  for (i in seq_along(resultslist)) resultslist[[i]]$gene_id <- rowData(dds)$gene_id
  saveRDS(resultslist, here("pipeline/r_data/resultslist.rds"))

  ddstodf <- function(ddsdf) {
    results(ddsLRT) %>%
      as.data.frame() %>%
      rownames_to_column("gene_id")
  }
  ddsLRT <- DESeq(dds, test = "LRT", full = design(dds), reduced = ~assay)
  changegenes <- results(ddsLRT) %>%
    ddstodf() %>%
    subset() %>%
    filter(padj < 0.05) %>%
    .$gene_id

  saveRDS(changegenes, here("pipeline/r_data/changegenes.rds"))

  if (do_xtail) {
    run_xtail <- function() Sys.glob(here("pipeline", "xtail/*")) %>% setNames(., basename(.) %>% str_extract(regex("(?<=xtail_).*(?=\\.tsv)")))
    # xtailcontrasts = 'condition'
    # xtailfiles = run_xtail(dds)
    xtailfiles <- run_xtail()
    xtailres <- map_df(.id = "contrast", xtailfiles, fread)
    xtailres$isxtail <- TRUE
    xtailres$contrast <- paste0("xtailTE_", xtailres$contrast)

    ddsreslist <- resultslist[!str_detect(names(resultslist), "xtail")]
    # xtailres%<>%rename(log2FoldChange := log2fc)
    xtailres %<>% rename(gene_id := feature_id)
    xtailres %<>% rename(baseMean := base_mean)
    resultslist <- c(ddsreslist, split(xtailres, xtailres$contrast))
  }


  # resultslist =
} else {
  resultslist <- readRDS(here("pipeline/r_data/resultslist.rds"))
  changegenes <- readRDS(here("pipeline/r_data/changegenes.rds"))
}

cat("successfully calculated contrasts:")
message(names(resultslist) %>% paste(collapse = " "))
################################################################################
# This script demonstrates how to use the RNA-Seq and ATAC-Seq data we 
# collected to map eQTL and caQTL.
# To run, download the figshare data from the following URL:
# https://doi.org/10.6084/m9.figshare.12233570 and unzip in the 
# "figshare_data" directory
################################################################################
library(tidyverse)
library(assertthat)
library(sva)            # for ComBat
library(qtl2)           # https://kbroman.org/qtl2/
library(qtl2convert)    # https://github.com/rqtl/qtl2convert
n.cores <- 1            # could use multiple CPUs for speed

# Load datasets that will be used in mapping
gmap <- readRDS("../data/gmap.Rds")
covar_atac <- read.table("../data/covar_atac.tsv")
covar_rna <- read.table("../data/covar_rna.tsv")
probs <- readRDS("../figshare_data/genotype_probs.Rds")
rna <- read.table("../figshare_data/counts_rna_norm_DO.tsv.gz")    # Gene expression is upper quartile normalized
atac <- read.table("../figshare_data/counts_atac_norm_DO.tsv.gz")  # Peak accessibility counts are TMM-normalized

# Only examine samples present in both RNA and ATAC datasets
matches <- read_tsv("../data/sample_matching.tsv")
matches$id <- make.unique(matches$RNA)
assert_that(all(matches$ATAC %in% colnames(atac)))
assert_that(all(matches$RNA %in% colnames(rna)))
assert_that(are_equal(rownames(probs), matches$top_match))
rownames(probs) <- matches$id
rna <- rna[, matches$RNA]
colnames(rna) <- matches$id
atac <- atac[, matches$ATAC]
colnames(atac) <- matches$id
covar_rna <- covar_rna[matches$RNA, , drop=FALSE]
rownames(covar_rna) <- matches$id
covar_atac <- covar_atac[matches$ATAC, , drop=FALSE]
rownames(covar_atac) <- matches$id

# ComBat normalize for RNA library prep batch.
dat <- log(rna + 1)
mod <- model.matrix(~ sex, data=covar_rna)
rnaComBat <- ComBat(dat=dat, batch=covar_rna$libraryprep, mod=mod, 
    par.prior=TRUE, prior.plots=FALSE)
covar_rna <- covar_rna[, "sex", drop=FALSE]

# RankZ normalize.
# Samples in rows, gene expression/chromatin accessibility in columns
rankZ <- function (x) {
    x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
    qnorm(x)
}
rnaZ <- apply(rnaComBat, 1, rankZ)
atacZ <- apply(atac, 1, rankZ)

# Prepare data.frame giving physical and genetic positions of each marker
markers <- dimnames(probs)[[3]]
map_dat <- qtl2::read_csv("../data/ref_genome_grid_69k.csv")[markers, ]
map <- qtl2:::split_map(map_dat)
map_dat$marker <- rownames(map_dat)
map_dat$pos_cM <- as.numeric(map_dat$cM)
map_dat$pos_bp <- as.integer(sapply(strsplit(map_dat$marker, "_"), "[[", 2))

# Convert genotype probabilities to qtl2 format
# Calculate kinship matrix
probs2 <- probs_doqtl_to_qtl2(probs, map_dat,
    chr_column="chr", pos_column="pos_cM", marker_column="marker")
kinship_loco <- calc_kinship(probs2, "loco", cores=n.cores)

# Do mapping in batches of genes/ATAC peaks.
# Save only QTL peaks with LOD > 5
batchmap <- function(nbatch, normdata, covar, thrA=5, thrX=5, ...) {
    ss <- round(seq(0, ncol(normdata), length.out=nbatch + 1))
    peaks <- list()
    for (i in 1:nbatch) {
        start <- ss[i] + 1
        end <- ss[i + 1]
        out <- scan1(probs2, normdata[, start:end, drop=FALSE],
            kinship_loco, addcovar=covar, cores=n.cores, ...)
        peaks[[i]] <- find_peaks(out, gmap, drop=1.5,
            threshold=thrA, thresholdX=thrX)   # returns a long & tidy dataset
    }
    do.call('rbind', peaks) %>% select(-lodindex) %>%
        rename(phenotype=lodcolumn, peak_chr=chr, peak_cM=pos)
}

# Map a few genes and ATAC peaks in this test script
# ... Should finish in a few minutes.
features <- 1:10
eQTL <- batchmap(nbatch=2, normdata=rnaZ[, features], covar=covar_rna)
caQTL <- batchmap(nbatch=2, normdata=atacZ[, features], covar=covar_atac)

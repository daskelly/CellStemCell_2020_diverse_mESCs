################################################################################
# This script gives an example showing how to perform a mediation analysis. 
# In this example we will mediate the chromatin accessibility peaks in a caQTL 
# hotspot using as candidate mediators all transcripts whose genes are located 
# within the hotspot boundaries.
################################################################################
library(intermediate)   # See https://github.com/simecek/intermediate

# First run map_qtl.r in this directory
# or load the results from that script
mediation_dat <- rnaZ   # rankZ-normalized transcript abundances

# Load sample peaks
sample_peaks <- read_tsv("../data/sample_atac_peaks.tsv")
transcripts_in_hotspot <- scan("../data/sample_transcripts_in_hotspot.txt", what="")
# Include intercept in the matrix of covariates
covar_atacI <- as.data.frame(covar_atac)
covar_atacI$intercept <- 1
covar_atacI <- as.matrix(covar_atacI)
# Load gene annotations:
annot <- read_tsv('../data/gene_annotation_Ensembl_v82.tsv', col_types='cciicc') %>%
    mutate(id=ensembl_gene_id, chr=chrom, symbol=mgi_symbol,
    mid=(start+end)/2, pos=mid/1e6) %>% # for intermediate
    filter(id %in% transcripts_in_hotspot) %>%
    arrange(chrom, mid)
assert_that(all(annot$id %in% colnames(mediation_dat)))
mediators <- mediation_dat[, annot$id, drop=FALSE]

# Do the mediation
mediation_results_list <- list()
for (i in 1:nrow(sample_peaks)) {
    target <- atacZ[, sample_peaks$id[i], drop=FALSE]
    # map without including candidate mediator in order to get top marker
    before <- scan1(probs2[, sample_peaks$peak_chr[i]], target,
        addcovar=covar_atac[rownames(target), , drop=FALSE])
    closest_marker <- rownames(before)[which.max(before)]
    qtl_geno <- probs[, LETTERS[2:8], closest_marker]
    # Do a mediation scan with all candidate mediators (transcripts in hotspot)
    med <- mediation.scan(target=target[, 1], mediator=mediators,
        annotation=annot, covar=covar_atacI[rownames(target), , drop=FALSE],
        qtl.geno=qtl_geno, method="double-lod-diff", verbose=FALSE)
    mediation_results_list[[i]] <- select(med, id, symbol, LOD) %>%
        rename("mediator_id"="id", "mediator_symbol"="symbol", "after_LOD"="LOD") %>%
        mutate(before_LOD=max(as.numeric(before)), target_id=sample_peaks$id[i])
}
mediation_results <- bind_rows(mediation_results_list) %>%
    mutate(LOD_diff=before_LOD - after_LOD)
# results is a data.frame with results of testing each transcript as a 
# candidate mediator for each peak in sample_peaks.

# Find best mediator for each peak. Tally up the number of times
# each transcript is the top mediator.
group_by(mediation_results, target_id) %>%
    top_n(1, LOD_diff) %>% group_by(mediator_symbol) %>%
    tally() %>% filter(n > 5) %>% arrange(desc(n)) %>% print()

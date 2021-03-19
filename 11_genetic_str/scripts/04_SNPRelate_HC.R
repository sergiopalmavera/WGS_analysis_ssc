# Checking example before running SNPRelate on my data
## According to https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate
## To clarify axes labeling: https://support.bioconductor.org/p/119389/

# Load packages
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

genofile <- seqOpen("../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.gds")

ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile))

sample_info <- read.table("../../sample_info/internal_external_sample_info.tsv", sep = "\t", header = TRUE)
head(sample_info)

# arrange sample info ids as in dendrogram
idx <- match(ibs.hc$sample.id, sample_info$sample_id)
sample_info <- sample_info[idx,]

# Make sure that samples have the same order before assigning groups
data.frame(x=ibs.hc$sample.id,y=sample_info$sample_id)

rv <- snpgdsCutTree(ibs.hc, samp.group=as.factor(sample_info$breed))
rv

saveRDS(rv$dendrogram, "../output/dendrogram_full_snp_set.rds")

#saveRDS(as.factor(sample_info$Linie), "../data/Linie.rds")

sessionInfo()

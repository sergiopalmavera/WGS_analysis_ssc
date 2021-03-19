# Checking example before running SNPRelate on my data
## According to https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#installation-of-the-package-snprelate

# Load packages
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

# convert VCF to GDS
vcf <- "../../10_FinalVCF/output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.vcf"
seqVCF2GDS(vcf, "../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_header_fixed.gds")


#seqVCF2GDS("tst.vcf", "tst.gds")



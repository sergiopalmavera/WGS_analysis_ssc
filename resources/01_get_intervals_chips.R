# Collecting chip information 

# Accoring to latest genome reference publication (https://doi.org/10.1093/gigascience/giaa051), there are 4 chips:
#(1) Affymetrix Axiom™ Porcine Genotyping Array
#(2) Illumina PorcineSNP60 
#(3) Gene Seek Genomic Profiler Porcine—HD beadChip
#(3) Gene Seek Genomic Profiler Porcine v2—LD Chip


######################################################################
# Load Affimetrix sites ("Axiom Porcine Genotyping Array")
# It contains 658,692 markers 
# Website: https://www.thermofisher.com/order/catalog/product/550588
# It targets Pietrain, Duroc, Large White and Landrace among others
# Reference genome: susScr11
# Download requires registration (serpalma.v@gmail.com sergio1983)
#######################################################################

library(dplyr)
library(vroom)

affi_dat <- vroom("resources/Axiom_PigHD/Axiom_PigHD_v1.na35.r4.a2.annot_NoHeader.csv")

affi_dat$Genome %>% unique()

affi_dat <- affi_dat %>% 
  dplyr::select(Chromosome, `Physical Position`) %>% 
  dplyr::rename(chr = Chromosome, pos = `Physical Position`)


###################################################################################################################################
# Illumina PorcineSNP60
# Genome version: 11.1
# Duroc, Landrace, Pietran, and Large White included https://www.illumina.com/products/by-type/microarray-kits/porcine-snp60.html
# Data (manifest file) in 
# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/porcinesnp60/porcinesnp60_v2productfiles/PorcineSNP60v2_15031945_C2_csv.zip
####################################################################################################################################


illumina60K_dat <- vroom("resources/PorcineSNP60v2/PorcineSNP60v2_15031945_C2.csv", skip = 7)

illumina60K_dat$GenomeBuild %>% unique()

# https://support.illumina.com/bulletins/2016/05/infinium-genotyping-manifest-column-headings.html
# MapInfo: Chromosomal coordinates of the SNP.

illumina60K_dat <- illumina60K_dat %>% dplyr::select(Chr, MapInfo) %>% dplyr::rename(chr = Chr, pos = MapInfo)
    
########################################################################################################################################################  
# GenSeek Illlumina chips
# Targets main porcine breeds, but it did not specified which ones.
# based on genome version 10.2
# https://support.illumina.com/array/array_kits/geneseek-ggp-arrays/downloads.html
# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/geneseek-ggp/geneseek-ggp-porcine-hd-manifest-file-csv.zip
# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/geneseek-ggp/geneseek-ggp-porcine-ld-manifest-file-csv.zip
########################################################################################################################################################

ggp_hd <- vroom("resources/GeneSeek/geneseek-ggp-porcine-hd-manifest-file-csv/GGP_HD_Porcine.csv", skip = 7)
ggp_ld <- vroom("resources/GeneSeek/geneseek-ggp-porcine-ld-manifest-file-csv/GGP_LD_Porcine.csv", skip = 7)

ggp_hd$GenomeBuild %>% unique()
ggp_ld$GenomeBuild %>% unique()

# make bed file for liftover
# https://www.biostars.org/p/359535/
bind_rows(
  ggp_hd %>% dplyr::select(Chr, MapInfo),
  ggp_ld %>% dplyr::select(Chr, MapInfo)
  ) %>% 
  unique() %>% 
  mutate(tmp = MapInfo -1) %>% # coordinates are 1-based. bed uses 0-based coordinates. That's why we have to subtract 1 for the start coordinate.
  arrange(Chr, tmp, MapInfo) %>% 
  mutate(tmp2 = paste0(Chr, ":", tmp, "-", MapInfo)) %>% 
  dplyr::select(tmp2) %>% 
  write.table(file = "resources/GeneSeek/ggp_hd_ld_v10.2.bed",quote = FALSE, row.names = FALSE, col.names = FALSE)

#...liftover did not work. I will only work with axiom and 60K chips






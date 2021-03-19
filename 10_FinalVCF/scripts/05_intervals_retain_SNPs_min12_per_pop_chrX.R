#args=(commandArgs(TRUE))
#print(args)
#eval(parse(text=args[[1]]))
#paste("Processing chromosome:", chr)

# passing characters to args does not work - running X alone

chr <- "X"


library(dplyr)
library(vroom)
library(stringr)
options(scipen=999)

# Define vcf file
tab <- paste0("../output/cohort_biallelicSNPs_HardFiltered_WithMissingness_chr", chr,".tsv")
print(tab)

sample_info <- vroom("../../sample_info/internal_external_sample_info.tsv")
dim(sample_info)

# Load genotype table
gt_tab <- vroom(tab)
dim(gt_tab)

# Prepare col names
names(gt_tab) <- str_remove(names(gt_tab), ".GT") %>% str_remove("-L1")

# Match order of table and sample info
idx1 <- match(sample_info$sample_id, names(gt_tab))

tmp1 <- gt_tab[,idx1]
gt_tab <- gt_tab %>% 
   dplyr::select(CHROM,POS) %>% 
   bind_cols(tmp1) 

# Colnames match sample info?
identical(sample_info$sample_id,names(tmp1))

# Get sites with at least N non-missing samples
idx2 <- gt_tab %>% 
  dplyr::select(-CHROM,-POS) %>% 
  apply(1, function(x){

    #x=gt_tab %>% dplyr::select(-CHROM,-POS) %>% .[1,] %>% unlist()

    # count calls (non ./.) per group and check if N is >= 12 samples
    xx <- tapply(x, sample_info$breed, function(gt) sum(gt != "./.") >= 12)
    
    # check if all groups pass the min N per group without missingness 
    all(xx)
  })

sum(idx2)

# Define output file name
out_nm <- paste0("../output/keep_snps_min_12_nonmiss_per_pop_chr", chr, ".intervals")

# Get and export intervals for GATK
gt_tab %>% 
  dplyr::select(CHROM,POS) %>% 
  .[idx2,] %>% 
  mutate(tmp = paste0(CHROM,":",POS,"-",POS)) %>% 
  dplyr::select(tmp) %>% 
  write.table(out_nm, quote = F, col.names = F, row.names = F)

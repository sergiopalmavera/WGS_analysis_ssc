# Load Affimetrix sites ("Axiom Porcine Genotyping Array")
# It contains 658,692 markers 
# Website: https://www.thermofisher.com/order/catalog/product/550588
# Download requires registration (serpalma.v@gmail.com sergio1983)

library(dplyr)
library(vroom)

affi_dat <- read.csv("Axiom_PigHD_v1.na35.r4.a2.annot_NoHeader.csv")

affi_dat <- affi_dat %>% 
  dplyr::select(Chromosome, Physical.Position) %>% 
  mutate(locus = paste0(Chromosome,":",Physical.Position)) %>% 
  dplyr::select(locus)

# Export 
write.table(affi_dat, "affi_loci.tab", quote = F, col.names = F, row.names = F)  
    
  

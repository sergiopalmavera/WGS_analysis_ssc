library(dplyr)
library(ggplot2)
library(stringr)


fls <- list.files("../output", pattern = "20Kb_thinned.recode.*.plink.ld", full.names = TRUE)

dat <- lapply(fls, function(fl){
	        d <- read.table(fl, header = TRUE, stringsAsFactors = FALSE) 
	        d <- mutate(d, dist = BP_B - BP_A)
	        return(d)
		})

names(dat) <- fls %>% 
  basename() %>% 
  str_remove("cohort_biallelicSNPs_HardFiltered_WithMissingness_filtered_20Kb_thinned.recode.") %>% 
  str_remove(".plink.ld")

dat2 <- dat %>% bind_rows(.id = "breed")

dat2 <- dat2 %>% 
  mutate(breed = ifelse(breed == "duroc", "DU", breed),
         breed = ifelse(breed == "pietrain", "PI", breed),
         breed = ifelse(breed == "landrace", "LR", breed),
         breed = ifelse(breed == "large_white", "LW", breed),
         breed = ifelse(breed == "euro_wild_boar", "EWB", breed))


dat2$breed <- factor(dat2$breed, levels = c("EWB","LW","LR","PI","DU"))

png("../output/LD_trends_max5Mb.png", width = 1500, height = 1500, res = 300, units = "px")
ggplot(data = dat2, aes(x=dist/1000, y=R2, color = breed)) +
  geom_smooth(method = "auto", se = FALSE) +
    xlab("Pairwise distance (Kb)") +
    ylab(expression(LD~(r^{2}))) +
    theme_bw(base_size = 12)+
    scale_color_brewer(palette = "Dark2") +
    theme(legend.title=element_blank(), plot.title = element_text(hjust = 0.5))
dev.off()


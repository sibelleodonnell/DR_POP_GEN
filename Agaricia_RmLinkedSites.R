##### SETTING SITES FOR IBS #####
tag_site <- read.table("tag_cut1_site.txt", header = FALSE, sep = "\t")
colnames(tag_site)[colnames(tag_site) == "V1"] <- "Tag"
colnames(tag_site)[colnames(tag_site) == "V2"] <- "Site"
library(dplyr)
unique_tag_site <- tag_site %>%
  distinct(Tag, .keep_all = TRUE)
write.table(unique_tag_site, file = "unique_tag_site.txt", sep = "\t", row.names = FALSE)
#sftp unique_tag_site.txt back into terminal for setting sites in ANGSD


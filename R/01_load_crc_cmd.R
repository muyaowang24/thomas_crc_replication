library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

se_crc <- curatedMetagenomicAnalyses::makeSEforCondition(
  "CRC",
  removestudies = "HMP_2019_ibdmdb",
  dataType = "relative_abundance"
)

meta_crc <- as.data.frame(colData(se_crc))
abund_crc <- assay(se_crc)

write.csv(meta_crc, "data_processed/meta_crc.csv", row.names = FALSE)
saveRDS(abund_crc, "data_processed/abund_crc.rds")
saveRDS(se_crc, "data_processed/se_crc.rds")

rowdat <- as.data.frame(rowData(se_crc))

se_crc_species <- se_crc
meta_crc_species <- as.data.frame(colData(se_crc_species))
abund_crc_species <- assay(se_crc_species)

meta_crc_plot <- meta_crc %>% 
  mutate(study_condition = factor(study_condition, levels = c("control", "CRC")))

p_study_comp <- ggplot(meta_crc_plot, aes(x = study_name, fill = study_condition)) +
  geom_bar(position = "fill") +
  coord_flip() +
  labs(
    title = "Case-control composition across CRC studies",
    x = "Study",
    y = "Proportion",
    fill = "Condition"
  ) +
  theme_bw()

p_study_comp

ggsave(
  filename = "output/crc_study_composition.png",
  plot = p_study_comp,
  width = 9,
  height = 6,
  dpi = 300
)

saveRDS(se_crc_species, "data_processed/se_crc_species.rds")
write.csv(meta_crc_species, "data_processed/meta_crc_species.csv", row.names = FALSE)
saveRDS(abund_crc_species, "data_processed/abund_crc_species.rds")
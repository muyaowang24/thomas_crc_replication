library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

se_crc_species <- readRDS("data_processed/se_crc_species.rds")
meta_crc <- as.data.frame(SummarizedExperiment::colData(se_crc_species))
abund_crc <- SummarizedExperiment::assay(se_crc_species)
taxa <- as.data.frame(SummarizedExperiment::rowData(se_crc_species))

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

meta_crc$sample_id <- rownames(meta_crc)

abund_df <- as.data.frame(abund_crc)
abund_df$feature_id <- rownames(abund_df)

abund_long <- abund_df %>%
  pivot_longer(
    cols = -feature_id,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  left_join(meta_crc, by = "sample_id") %>%
  left_join(
    taxa %>%
      rownames_to_column("feature_id"),
    by = "feature_id"
  ) %>%
  mutate(
    abundance_log = log10(abundance + 1e-5),
    study_condition = factor(study_condition, levels = c("control", "CRC"))
  ) %>%
  filter(study_condition %in% c("control", "CRC"))

run_one_study <- function(dat) {
  dat %>%
    group_by(feature_id, species) %>%
    summarise(
      n_control = sum(study_condition == "control"),
      n_crc = sum(study_condition == "CRC"),
      mean_control = mean(abundance_log[study_condition == "control"], na.rm = TRUE),
      mean_crc = mean(abundance_log[study_condition == "CRC"], na.rm = TRUE),
      median_control = median(abundance_log[study_condition == "control"], na.rm = TRUE),
      median_crc = median(abundance_log[study_condition == "CRC"], na.rm = TRUE),
      diff_mean = mean_crc - mean_control,
      p_value = tryCatch(
        wilcox.test(
          abundance_log[study_condition == "CRC"],
          abundance_log[study_condition == "control"]
        )$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      direction = case_when(
        diff_mean > 0 ~ "enriched_in_crc",
        diff_mean < 0 ~ "depleted_in_crc",
        TRUE ~ "no_difference"
      )
    )
}

study_results <- abund_long %>%
  group_split(study_name) %>%
  setNames(unique(abund_long$study_name)) %>%
  imap_dfr(~ {
    res <- run_one_study(.x)
    res$study_name <- .y
    res
  }) %>%
  select(
    study_name, feature_id, species,
    n_control, n_crc,
    mean_control, mean_crc,
    median_control, median_crc,
    diff_mean, p_value, fdr, direction
  )

write.csv(
  study_results,
  "results/crc_species_univariate_by_study.csv",
  row.names = FALSE
)

top_hits <- study_results %>%
  filter(!is.na(fdr)) %>%
  group_by(study_name) %>%
  arrange(fdr, desc(abs(diff_mean)), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top_hits,
  "results/crc_species_univariate_top20_by_study.csv",
  row.names = FALSE
)

consistency_table <- study_results %>%
  mutate(sig = !is.na(fdr) & fdr < 0.05) %>%
  group_by(species) %>%
  summarise(
    n_studies_tested = n(),
    n_sig = sum(sig),
    n_enriched_crc = sum(sig & direction == "enriched_in_crc"),
    n_depleted_crc = sum(sig & direction == "depleted_in_crc"),
    mean_diff = mean(diff_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig), desc(abs(mean_diff)))

write.csv(
  consistency_table,
  "results/crc_species_consistency_table.csv",
  row.names = FALSE
)

top_consistent <- consistency_table %>%
  filter(n_sig > 0) %>%
  slice_head(n = 30)

p_consistency <- ggplot(
  top_consistent,
  aes(x = reorder(species, n_sig), y = n_sig)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top recurrent CRC-associated species across studies",
    x = "Species",
    y = "Number of studies with FDR < 0.05"
  ) +
  theme_bw()

p_consistency

ggsave(
  filename = "output/top_recurrent_crc_species.png",
  plot = p_consistency,
  width = 10,
  height = 8,
  dpi = 300
)
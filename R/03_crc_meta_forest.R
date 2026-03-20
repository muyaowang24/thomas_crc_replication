library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(meta)

study_results <- read.csv("results/crc_species_univariate_by_study.csv", stringsAsFactors = FALSE)

top_species <- c(
  "Parvimonas micra",
  "Gemella morbillorum",
  "Peptostreptococcus stomatis",
  "Dialister pneumosintes",
  "Fusobacterium nucleatum"
)

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
    taxa %>% tibble::rownames_to_column("feature_id"),
    by = "feature_id"
  ) %>%
  mutate(
    abundance_log = log10(abundance + 1e-5),
    study_condition = factor(study_condition, levels = c("control", "CRC"))
  ) %>%
  filter(
    study_condition %in% c("control", "CRC"),
    species %in% top_species
  )

effect_table <- abund_long %>%
  group_by(study_name, species) %>%
  summarise(
    n_control = sum(study_condition == "control"),
    n_crc = sum(study_condition == "CRC"),
    mean_control = mean(abundance_log[study_condition == "control"], na.rm = TRUE),
    mean_crc = mean(abundance_log[study_condition == "CRC"], na.rm = TRUE),
    sd_control = sd(abundance_log[study_condition == "control"], na.rm = TRUE),
    sd_crc = sd(abundance_log[study_condition == "CRC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pooled_sd = sqrt(((n_control - 1) * sd_control^2 + (n_crc - 1) * sd_crc^2) / (n_control + n_crc - 2)),
    smd = (mean_crc - mean_control) / pooled_sd,
    se_smd = sqrt((n_control + n_crc) / (n_control * n_crc) + (smd^2) / (2 * (n_control + n_crc - 2)))
  )

write.csv(effect_table, "results/crc_species_effect_sizes.csv", row.names = FALSE)

run_meta_one_species <- function(sp) {
  dat <- effect_table %>% filter(species == sp)
  
  m <- metagen(
    TE = smd,
    seTE = se_smd,
    studlab = study_name,
    data = dat,
    sm = "SMD",
    method.tau = "REML",
    method.random.ci = "HK"
  )
  
  meta_summary <- data.frame(
    species = sp,
    k = m$k,
    pooled_effect = m$TE.random,
    pooled_se = m$seTE.random,
    lower = m$lower.random,
    upper = m$upper.random,
    p_value = m$pval.random,
    i2 = m$I2,
    tau2 = m$tau^2
  )
  
  png(
    filename = file.path("output", paste0(gsub(" ", "_", sp), "_forest.png")),
    width = 1200,
    height = 900,
    res = 150
  )
  
  forest(
    m,
    sortvar = dat$smd,
    xlab = "Standardized mean difference (CRC vs control)",
    leftcols = c("studlab"),
    leftlabs = c("Study"),
    rightcols = c("effect", "ci"),
    main = sp
  )
  
  dev.off()
  
  list(meta_summary = meta_summary, model = m)
}

meta_results <- lapply(top_species, run_meta_one_species)

meta_summary_table <- bind_rows(lapply(meta_results, function(x) x$meta_summary)) %>%
  arrange(p_value)

write.csv(meta_summary_table, "results/crc_species_meta_summary.csv", row.names = FALSE)

p_meta <- ggplot(meta_summary_table, aes(x = reorder(species, pooled_effect), y = pooled_effect)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  labs(
    title = "Pooled cross-cohort effects for top recurrent CRC-associated species",
    x = "Species",
    y = "Random-effects pooled SMD"
  ) +
  theme_bw()

p_meta

ggsave(
  filename = "output/crc_species_meta_summary.png",
  plot = p_meta,
  width = 10,
  height = 7,
  dpi = 300
)
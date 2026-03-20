library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(pROC)
library(glmnet)
library(SummarizedExperiment)

se_crc_species <- readRDS("data_processed/se_crc_species.rds")
meta_crc <- as.data.frame(colData(se_crc_species))
abund_crc <- assay(se_crc_species)
taxa <- as.data.frame(rowData(se_crc_species))

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

meta_crc$sample_id <- rownames(meta_crc)

common_species <- taxa$species
common_species[is.na(common_species) | common_species == ""] <- rownames(taxa)
rownames(abund_crc) <- make.unique(common_species)

X <- t(log10(abund_crc + 1e-5))
X <- as.data.frame(X)
X$sample_id <- rownames(X)

dat <- meta_crc %>%
  select(sample_id, study_name, study_condition) %>%
  filter(study_condition %in% c("control", "CRC")) %>%
  mutate(y = ifelse(study_condition == "CRC", 1, 0)) %>%
  left_join(X, by = "sample_id")

feature_names <- setdiff(colnames(dat), c("sample_id", "study_name", "study_condition", "y"))

run_lodo <- function(test_study) {
  train_dat <- dat %>% filter(study_name != test_study)
  test_dat  <- dat %>% filter(study_name == test_study)
  
  x_train <- as.matrix(train_dat[, feature_names])
  y_train <- train_dat$y
  
  x_test <- as.matrix(test_dat[, feature_names])
  y_test <- test_dat$y
  
  cvfit <- cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = 1,
    nfolds = 10,
    type.measure = "auc"
  )
  
  pred_prob <- as.numeric(predict(cvfit, newx = x_test, s = "lambda.min", type = "response"))
  
  roc_obj <- roc(y_test, pred_prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  data.frame(
    test_study = test_study,
    n_test = nrow(test_dat),
    n_control = sum(test_dat$y == 0),
    n_crc = sum(test_dat$y == 1),
    auc = auc_val
  )
}

study_list <- sort(unique(dat$study_name))

lodo_results <- map_dfr(study_list, run_lodo)

write.csv(lodo_results, "results/crc_lodo_auc_results.csv", row.names = FALSE)

p_auc <- ggplot(lodo_results, aes(x = reorder(test_study, auc), y = auc)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Leave-one-dataset-out AUC across CRC studies",
    x = "Held-out study",
    y = "AUC"
  ) +
  theme_bw()

p_auc

ggsave(
  filename = "output/crc_lodo_auc.png",
  plot = p_auc,
  width = 9,
  height = 6,
  dpi = 300
)

mean_auc <- mean(lodo_results$auc, na.rm = TRUE)
median_auc <- median(lodo_results$auc, na.rm = TRUE)

summary_table <- data.frame(
  mean_auc = mean_auc,
  median_auc = median_auc,
  min_auc = min(lodo_results$auc, na.rm = TRUE),
  max_auc = max(lodo_results$auc, na.rm = TRUE)
)

write.csv(summary_table, "results/crc_lodo_auc_summary.csv", row.names = FALSE)
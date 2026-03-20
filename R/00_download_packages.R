if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "curatedMetagenomicData",
  "TreeSummarizedExperiment",
  "SummarizedExperiment"
))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("waldronlab/curatedMetagenomicAnalyses")
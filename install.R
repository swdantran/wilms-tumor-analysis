# Install dependencies for the Wilms Tumor RNA-seq analysis

cran_packages <- c(
  "ggplot2",
  "tidyverse",
  "magick",
  "Cairo",
  "ggrepel",
  "textshaping",
  "dplyr",
  "circlize",
  "pheatmap",
  "RColorBrewer",
  "stringr"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "DESeq2",
  "airway",
  "ComplexHeatmap",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "EnhancedVolcano",
  "clusterProfiler",
  "apeglm"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
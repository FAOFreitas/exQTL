install.packages("rcmdcheck")
library(rcmdcheck)
rcmdcheck("C:/Users/felip/OneDrive/Esalq/Doutorado/Complementares/exQTL")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Gviz")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", force = TRUE)

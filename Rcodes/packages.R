# Package names
packages <- c("ggplot2", "readxl", "dplyr", "tidyr", "ggfortify", 
              "reshape2", "knitr", "lubridate", "readr", "psy", "car", "patchwork", 
              "imputeMissings", "RcmdrMisc", "questionr", "vcd", "multcomp",
              "KappaGUI", "rcompanion", "FactoMineR", "factoextra", "corrplot",
              "ltm", "goeveg", "corrplot", "FSA", "MASS", "scales", "nlme", "psych",
              "ordinal", "lmtest", "ggpubr", "dslabs", "stringr", "assist", 
              "ggstatsplot", "forcats", "styler", "remedy",  
              "CellID", "tidyverse", "pander", "cluster", "gprofiler2", "parallel", 
              "SeuratObject", "Seurat", "htmltools")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
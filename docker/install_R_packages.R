################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################
# Specify the list of packages that need to be imported ########################
list.of.packages <- c("BiocManager",
                      "optparse", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr",
                      "rcompanion",
                      "argparse",
                      "scatterpie", 
                      "scales",
                      "rstatix",
                      "remotes",
                      "PPCI",
                      "writexl",
                      "devtools",
                      "svglite",
                      "ggpubr"
)

bioc.packages <- c("Seurat",
                   "SingleCellExperiment", 
                   "celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "scRepertoire", 
                   "sctransform", 
                   "progeny",
                   "powerTCR",
                   "scRepertoire",
                   "org.Hs.eg.db",
                   "org.Mm.eg.db",
                   "DESeq2",
                   "diptest",
                   "methylKit",
                   "dmrseq",
                   "bsseq",
                   "DSS",
                   "DMRcate"
                   )

# Check if packages are installed ##############################################

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages, update = FALSE, ask = TRUE)

# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)

# INSTALL FROM GITHUB MANUALLY
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
devtools::install_github("cole-trapnell-lab/monocle3")

# INSTALL CLUSTER PROFILER MANUALLY, FIX ERROR CONTACTING KEGG DATABASE.
if ("clusterProfiler" %in% installed.packages()[, "Package"] == FALSE) {
  devtools::install_github("YuLab-SMU/clusterProfiler", upgrade = "never")
} else {
  if (packageVersion("clusterProfiler") != "4.7.1.3" ){
    remove.packages("DOSE")
    remove.packages("clusterProfiler")
    devtools::install_github("YuLab-SMU/clusterProfiler", upgrade = "never")
  }
}
library("clusterProfiler")
# EOF ##########################################################################

################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################
# Specify the list of packages that need to be imported ########################
list.of.packages <- c("BiocManager",
                      "optparse", 
                      "tidyverse", 
                      "ggplot2", 
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr",
                      "argparse",
                      "PPCI",
                      "writexl",
                      "devtools",
                      "svglite",
                      "ggpubr"
)

bioc.packages <- c("Seurat", 'VGAM', 'DDRTree', 'HSMMSingleCell', 'combinat', 'fastICA', 'densityClust', 'limma', 'qlcMatrix', 'pheatmap', 'proxy', 'slam', 'viridis', 'biocViews')

# Check if packages are installed ##############################################

install.packages("/home/storage/offline_pkgs/BiocGenerics_0.30.0.tar.gz", type = "source", repos = NULL)
install.packages("sparsesvd")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages, update = FALSE, ask = TRUE)

install.packages("/home/storage/offline_pkgs/qlcMatrix_0.9.7.tar.gz", type = "source", repos = NULL)
install.packages("/home/storage/offline_pkgs/monocle_2.12.0.tar.gz", type = "source", repos = NULL)

# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)
# EOF ##########################################################################

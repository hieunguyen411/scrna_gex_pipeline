gc()
rm(list = ls())

set.seed(42)

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(path.to.pipeline.src, "s10_PROGENy_pathway_analysis.R"))
source("/home/hieunguyen/CRC1382/src_2023/APanyot/00_import_libraries.R")
source("/home/hieunguyen/CRC1382/src_2023/APanyot/00_helper_functions.R")
library(progeny)

analysis.round <- "1st_round"
pct_mito <- 25

path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

PROJECT <- "APanyot_rerun_cellranger_with_Mu_rerun"

input.datadir <- file.path(path.to.outdir, sprintf("%s/%s/pct_mito_%s", PROJECT, analysis.round, pct_mito))

path.to.output <- file.path(input.datadir, "data_analysis")
path.to.01.output <- file.path(path.to.output, "01_output")

path.to.05.output <- file.path(path.to.output, "05_output")
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(path.to.01.output, "APanyot_dataset.output.s8.recluster_CellType.rds"))

##### INPUT ARGS
cluster.name <- "celltype" # default value
progeny.params <- list(scale = FALSE, 
                      species = "Mouse",
                      top = 500,
                      perm = 1,
                      return_assay = TRUE)

s.obj <- s10.PROGENy_pathway_analysis(s.obj = s.obj, 
                                      path.to.output = ".", 
                                      progeny.params = progeny.params,
                                      cluster.name = cluster.name, 
                                      save.RDS.s10 = FALSE, 
                                      PROJECT = PROJECT,
                                      scran_or_seurat = "scran")

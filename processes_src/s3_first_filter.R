s3.filter <- function(s.obj,
                      PROJECT,
                      path.to.output,
                      save.RDS.s3, 
                      nFeatureRNAfloor, 
                      nFeatureRNAceiling,
                      nCountRNAfloor, 
                      nCountRNAceiling,
                      pct_mitofloor, 
                      pct_mitoceiling,
                      pct_ribofloor, 
                      pct_riboceiling,
                      ambientRNA_thres = 0.5){
  
  # s.obj <- subset(s.obj, subset = nFeature_RNA > nFeatureRNAfloor &
  #                   nFeature_RNA < nFeatureRNAceiling&
  #                   nCount_RNA > nCountRNAfloor&
  #                   nCount_RNA < nCountRNAceiling &
  #                   percent.mt > pct_mitofloor&
  #                   percent.mt < pct_mitoceiling &
  #                   percent.ribo > pct_ribofloor &
  #                   percent.ribo < pct_riboceiling)
  
  s.obj <- subset(s.obj, subset = percent.mt < pct_mitoceiling)
  s.obj <- subset(s.obj, subset = AmbientRNA < ambientRNA_thres)
  
  if (save.RDS.s3 == TRUE){
    dir.create(file.path(path.to.output, "s3_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s3_output", 
                             paste0(PROJECT, ".output.s3.rds")))
  }
  return(s.obj)
}
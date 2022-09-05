s7.cellFactorRegressOut <- function(s.obj,
                                    path.to.output, 
                                    save.RDS.s7,
                                    PROJECT){
  all.genes <- rownames(x = s.obj)
  
  s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
  s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
  
  g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
  g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
  
  
  ### Scale data add exclude
  s.obj$CC.Difference <- s.obj$S.Score - s.obj$G2M.Score
  s.obj <- ScaleData(s.obj, vars.to.regress = c("G2M.Score", 
                                                "S.Score", 
                                                "percent.mt", 
                                                "percent.ribo",
                                                "nCount_RNA"), features = rownames(s.obj))
  # s.obj <- RunPCA(s.obj, features = c(s.genes, g2m.genes), nfeatures.print = 10, reduction.name="CELLCYCLED_PCA")
  if (save.RDS.s7 == TRUE){
    dir.create(file.path(path.to.output, "s7_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s7_output", 
                             paste0(PROJECT, ".output.s7.rds")))
  }
  return(s.obj)
}

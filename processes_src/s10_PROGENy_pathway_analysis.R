s10.PROGENy_pathway_analysis <- function( s.obj, 
                                          path.to.output, 
                                          progeny.params,
                                          cluster.name, 
                                          save.RDS.s10,
                                          scran_or_seurat = "scran",
                                          PROJECT){
  #### FUNCTION FORKED FROM https://github.com/CostaLab/scrna_seurat_pipeline/blob/master/data_factory.R
  Idents(s.obj) <- cluster.name
  
  s.obj <- progeny::progeny(s.obj, scale = progeny.params$scale, 
                                    organism = progeny.params$species, 
                                    top = progeny.params$top, 
                                    perm = progeny.params$perm, 
                                    return_assay = progeny.params$return_assay)
  s.obj <- ScaleData(s.obj, features = row.names(s.obj))
  DefaultAssay(s.obj) <-'progeny'
  all.pathways <- rownames(s.obj@assays$progeny)
  
  Idents(s.obj) <- cluster.name
  
  if (is.factor(s.obj@meta.data[, cluster.name])){
    s.obj@meta.data[, cluster.name] <- droplevels(s.obj@meta.data[, cluster.name])
  }else{
    c_names <-unique(s.obj@meta.data[, cluster.name])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    s.obj@meta.data[, cluster.name] <- factor(s.obj@meta.data[, cluster.name],
                                                      levels=sort(help_sort_func(c_names)))
  }
  

  if (scran_or_seurat == "scran"){
    diff.pathways <- data.frame()
    all.clusters <- unique(s.obj[[cluster.name]][[cluster.name]])
                           
    for (i in sort(all.clusters)){
      # get a vector storing all cluster ids
      g <- as.character(s.obj@meta.data[, cluster.name])
      # replace any cluster other than i to "others"
      g[g != i] <- "others"
      # convert the vector g to factor-vector
      g <- factor(g, levels = c(i, "others"))
      
      # run the function FindMarkers from scran (why?)
      tmp <- scran::findMarkers(as.matrix(s.obj@assays$progeny@data), g)[[1]] %>%
        as.data.frame() 
      # note that at this point we keep all pathways, FDR has not been done yet!
      # calculate the EFFECT SIZE of WILCOXON test
      r <- sapply(all.pathways, 
                  function(pw) rcompanion::wilcoxonR(as.vector(s.obj@assays$progeny@data[pw,]), g))
      names(r) <- to_vec(for (item in names(r)) str_replace(item, "[.]r", "")[[1]][[1]]) 
      # add the effect size vector r to the dataframe tmp
      tmp[names(r), "r"] <- r
      tmp <- tmp[names(r), ] %>% rownames_to_column("pathway")
      tmp$cluster <- i
      diff.pathways <- rbind(diff.pathways, tmp)
    }
    diff.pathways$p_val_adj <- p.adjust(diff.pathways$p.value, method = "fdr")
    diff.pathways <- subset(diff.pathways, select = c(pathway, p.value, p_val_adj, cluster, summary.logFC))
    colnames(diff.pathways) <- c("pathway", "p_val", "p_val_adj", "cluster", "logFC")
  } else {
    diff.pathways <- FindAllMarkers(s.obj, assay = "progeny", slot = "data", test.use = "wilcox") %>%
      subset(p_val_adj <= 0.05)
    row.names(diff.pathways) <- NULL
    colnames(diff.pathways) <- c("p_val", "logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "pathway")
    diff.pathways <- subset(diff.pathways, select = c(pathway, p_val, p_val_adj, cluster, logFC))
  }
  
  # save the differential pathway activity score analysis to the main seurat object
  s.obj@misc$progeny.diff.pathways <- diff.pathways
  
  if (save.RDS.s10 == TRUE){
    dir.create(file.path(path.to.output, "s10_output"), showWarnings = FALSE)
    saveRDS(object = cluster.markers, 
            file = file.path(path.to.output, "s10_output", 
                             paste0(PROJECT, ".PROGENy.s10.rds")))
  }
  return(s.obj)
}
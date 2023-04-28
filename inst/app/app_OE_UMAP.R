


observeEvent(input$run_UMAP_CS_batch, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    # envv$Remove_comps = c(13, 40)
    SDAres <- envv$SDAres
    if(length(envv$Remove_comps) > 50) {
      # suffix <- paste(envv$Remove_comps, collapse = "")
      
      suffix <- paste0("LargeSetOfComps_", length(envv$Remove_comps))
      
    } else if(length(envv$Remove_comps) > 30 & length(envv$Remove_comps) <= 50 ) {
      
      suffix <- paste0(findIntRuns(as.numeric(unlist(strsplit(envv$Remove_comps, ",")))), collapse="")
      
    } else if (length(envv$Remove_comps) <= 30) {
      suffix <- paste(envv$Remove_comps, collapse = "_")
    }
    
    UMAP_n.iter <- as.numeric(input$UMAPiter) # UMAP_n.iter = 1000
    # UMAP_n.iter = 300
    
    UMAPspread <- as.numeric(input$UMAPspread) # UMAP_n.iter = 1000
    
    
    # UMAP_pp <- as.numeric(input$UMAPpp) # UMAP_pp = 50
    
    choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    selected <- setdiff(choice, envv$Remove_comps)
    
    
    if(envv$Origin == "folder"){
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
    }
    if(envv$Origin == "prime"){
      head.path <- "sda_results"
      base.path <- ""
      #TODO
    }
    
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_UMAP_CellScore_QCfil_",UMAP_n.iter,"_",UMAPspread, "_", suffix, ".rds"))){
      
      envv$InfoBox_sub = "Starting UMAP with cell scores - batch-removal.. wait"
     
      UMAP_CS_batch = scCustFx:::RunUMAP.Matrix(DGEmat=  SDAres$scores[,selected],
                                         n_threads = num_threads,
                                         assay = NULL,
                                         n.neighbors = 60, #40
                                         n.components = 2L,
                                         metric = "cosine",
                                         n.epochs = UMAP_n.iter,
                                         learning.rate = 1.0,
                                         min.dist = 0.2,
                                         spread = UMAPspread,
                                         set.op.mix.ratio = 1.0,
                                         local.connectivity = 1L,
                                         repulsion.strength = 1,
                                         negative.sample.rate = 5,
                                         a = NULL,
                                         b = NULL,
                                         seed.use = 66,
                                         metric.kwds = NULL,
                                         angular.rp.forest = FALSE,
                                         reduction.key = 'UMAPSDA_',
                                         verbose = TRUE)
      
      
      
      envv$InfoBox_sub = "Saving UMAP with cell scores - batch-removal.. wait"
      saveRDS(UMAP_CS_batch, file=paste0(envv$path2SDA_dyn, "/", head.path,"_UMAP_CellScore_QCfil_",UMAP_n.iter,"_",UMAPspread, "_", suffix, ".rds"))
      
      
      
    } else {
      UMAP_CS_batch <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_UMAP_CellScore_QCfil_",UMAP_n.iter,"_",UMAPspread, "_", suffix, ".rds"))
    }
    envv$InfoBox_sub = "batch-removal UMAP with cell scores complete"
    envv$UMAP_CS_batch <- UMAP_CS_batch
    
  }
  
  
})


output$UMAP_CS_batch1 <- renderPlot({
  
  if(is.null(envv$UMAP_CS_batch)){
    plot(x=0, y=0, main="UMAP CS Batch not found")
    
  } else {
    
    UMAPDF <- as.data.frame(envv$UMAP_CS_batch)
    rownames(UMAPDF)  <- rownames(envv$SDAres$scores)
    colnames(UMAPDF) <- c("UMAP1_batch", "UMAP2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      UMAPDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      UMAPDF$SumScore <- UMAPDF$SumScore/mean(UMAPDF$SumScore)
      
      ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("UMAP SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      
      UMAPDF$Meta <- MetaDF[rownames(UMAPDF), input$Metaselect1]
      
      ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch, color=Meta)) +
        geom_point(size = 1, alpha=.4)+ theme_bw() +
        theme(legend.position = "bottom", aspect.ratio=1) +
        ggtitle(paste0("UMAP - batch removed cell scores\n", input$Metaselect1)) +
        scale_color_manual(values = col_vector
                           #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
        ) + 
        guides(colour = guide_legend(override.aes = list(size = 2, alpha=1), ncol=5))
      
      
    }
  }
  
})



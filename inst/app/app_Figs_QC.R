output$SDAqcMaxScorefilt <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    MaxScore.thr <- quantile(SDAres$component_statistics$max_score, c(.95))
    MaxScore.thr9 <- quantile(SDAres$component_statistics$max_score, c(.9))
    MaxScore.thr75 <- quantile(SDAres$component_statistics$max_score, c(.75))
    
    sum(SDAres$component_statistics$max_score > MaxScore.thr)
    
    ggplot(SDAres$component_statistics, aes(y=max_score)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                         outlier.size = 2) + 
      geom_abline(slope=0, intercept = MaxScore.thr, colour="red") + 
      geom_abline(slope=0, intercept = MaxScore.thr9, colour="dodgerblue")+ 
      geom_abline(slope=0, intercept = MaxScore.thr75, colour="navy")+ 
      theme_bw() + ggtitle(paste0("Max_Score quantile-thresholds remove:\n 95th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr),
                                  "\n 90th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr9),
                                  "\n 75th perc: ", sum(SDAres$component_statistics$max_score > MaxScore.thr75)))
    
    
  }
  
  
})

output$SDAqc1 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    ggplot(SDAres$component_statistics, aes(max_score, max_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw()+ ggtitle("")}
  
  
})

output$SDAqc2 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    ggplot(SDAres$component_statistics, aes(max_score, mean_score,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw()+ ggtitle("")}
  
  
})

output$SDAqc3 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    ggplot(SDAres$component_statistics, aes(mean_score, mean_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw() + ggtitle("")}
  
  
})

output$SDAqc4 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    ggplot(SDAres$component_statistics, aes(sd_score, sd_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw()}
  
  
})

output$SDAqc5 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    ggplot(SDAres$component_statistics, aes(ssqrd_score, ssqrd_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw() + ggtitle("")
  }
  
  
})

output$SDAqc6 <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    ggplot(SDAres$component_statistics, aes(max_loading, mean_loading,
                                            label = Component_name_plot)) +
      geom_point() + geom_label_repel() + theme_bw() + ggtitle("")
  }
  
  
  
})

output$convergence <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    hist(as.numeric(SDAres$loadings[[1]]), 
         main = "Overall distribution of gene loadings\nrep1", 
         xlab = "Gene Loading", breaks = 300, xlim = range(-.15,.15), col="skyblue")
  }
  
  
  
})

output$loadhist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    SDAtools::check_convergence(SDAres)
  }
  
  
  
})

output$scoredist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    qplot(asinh(as.numeric(SDAres$scores)), 
          binwidth = 0.01, main = "Overall distribution of individual scores", 
          xlab = "asinh(Score)")+ scale_y_log10()
  }
  
  
  
})

output$pipdist <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    SDAres <- envv$SDAres
    
    qplot(as.numeric(SDAres$pips[[1]]), geom = "histogram", 
          binwidth = 0.005) + xlab("PIP") + ylab("Count") + scale_y_log10()
  }
  
  
  
})

output$slackslabprior <- renderPlot({
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    SDAres <- envv$SDAres
    
    
    library(ggforce)
    
    density_data <- rbind(data.table(density=as.data.table(density(SDAres$loadings[[1]][SDAres$pips[[1]]>0.5], bw=1e-4, n=5000)[c("x","y")]), type="Slab (PIP>0.5)"),
                          data.table(density=as.data.table(density(SDAres$loadings[[1]][SDAres$pips[[1]]<0.5],bw=1e-4, n=5000)[c("x","y")]), type="Spike (PIP<0.5)"))
    
    setnames(density_data, c("gene_loading","density","type"))
    
    
    sparsity_plot <- ggplot(density_data, aes(gene_loading, density, colour=type)) +
      geom_line() +
      facet_zoom(xy = density<40 & abs(gene_loading)<0.1) +
      scale_color_brewer(palette = "Set1") + 
      theme_bw() +
      theme(legend.title=element_blank()) +
      labs(x="Gene Loading",y="Density") +
      scale_x_continuous(labels = function(x) as.character(x)); 
    sparsity_plot
    
  }
  
  
})

## get biomart
observeEvent(input$getGeneAnn, {
  
  envv$InfoBox_sub = "Getting gene annotations"
  
  envv = ShinySDA:::Run_GeneAnn_evv(envv, input = input)
  
  
  
})
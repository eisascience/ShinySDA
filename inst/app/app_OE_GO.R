## get GO
observeEvent(input$getSDAGo, {
  
  envv$InfoBox_sub <- "Stating GO"
  
  envv = ShinySDA:::Run_GO_evv(envv, input)
  
  
  
})
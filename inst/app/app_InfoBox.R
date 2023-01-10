output$InfoBox <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_sub,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})

output$InfoBox_Main <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_sub,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})

output$InfoBox_Prime <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_sub,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})

output$InfoBox_Folder <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_sub,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})


output$InfoBox_PP <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_sub,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})
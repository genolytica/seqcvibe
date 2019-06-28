genomeBrowserTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    return(NULL)
}

genomeBrowserTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    return(NULL)
}

genomeBrowserTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    output$genomeBrowser <- renderUI({
        tags$iframe(
            src=loadJBrowse(
                source=as.character(currentMetadata$source),
                dataset=as.character(currentMetadata$dataset),
                config=currentMetadata$final,
                org=currentMetadata$genome
            ),
            name="JBrowse",seamless=NA,
            height="800px",width="100%"
        )
    })

}

genomeBrowserTabPanelObserve <- function(input,output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    genomeBrowserTabPanelReactiveEvents <- 
        diffExprTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    genomeBrowserTabPanelReactive <- diffExprTabPanelEventReactive(input,output,
        session,allReactiveVars,allReactiveMsgs)
        
    genomeBrowserTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
}

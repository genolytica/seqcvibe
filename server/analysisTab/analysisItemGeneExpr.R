geneExprTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    currentBoxplot <- allReactiveVars$currentBoxplot
}

geneExprTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

geneExprTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentBoxplot <- allReactiveVars$currentBoxplot
}

geneExprTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}
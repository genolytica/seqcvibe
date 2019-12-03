# server.R
library(auth0)

# Load required libraries
source("config/init_server_globals.R")

#auth0_server(
function(input,output,session) {
    # Load init packages script
    source("config/init_packages.R")
    # Load SeqCVIBE libs
    source("server/reactiveVars.R",local=TRUE)
    source("server/dataSelectorTab/dataSelectorItem.R",local=TRUE)
    source("server/dataSelectorTab/sessionLoaderItem.R",local=TRUE)
    source("server/signalViewerTab/signalViewerItemGeneSignal.R",local=TRUE)
    source("server/signalViewerTab/signalViewerItemAreaSignal.R",local=TRUE)
    source("server/expressionViewerTab/expressionViewerItemKnownGene.R",
        local=TRUE)
    source("server/expressionViewerTab/expressionViewerItemCalculator.R",
        local=TRUE)
    source("server/analysisTab/analysisItemDiffExpr.R",local=TRUE)
    source("server/analysisTab/analysisItemClustering.R",local=TRUE)
    source("server/analysisTab/analysisItemCorrelation.R",local=TRUE)
    source("server/analysisTab/analysisItemMdsPca.R",local=TRUE)
    source("server/genomeBrowserTab/genomeBrowserItem.R",local=TRUE)
    
    # Init packages
    initPackages(session)
    
    #assign("session",session,envir=.GlobalEnv)
    #assign("input",input,envir=.GlobalEnv)
    #assign("output",output,envir=.GlobalEnv)
    
    # Make %#^%$^%$@( globals visible AND changeable
    makeReactiveBinding("loadedGenomes")
    makeReactiveBinding("loadedData")
    
    # Initialize all the reactive variables used...
    allReactiveVars <- initReactiveVars()
    # ...and reactive messages
    allReactiveMsgs <- initReactiveMsgs()
    
    # Data selector
    dataSelectorTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
        
    # Session loader
    sessionLoaderTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
    
    # Signal viewer - Gene signal
    geneSignalTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
   
    # Signal viewer - Area signal
    areaSignalTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
    
    # Expresion viewer - Known genes
    expressionExplorerTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
    
    # Expresion viewer - Calculator
    expressionCalculatorTabPanelObserve(input,output,session,
        allReactiveVars,allReactiveMsgs)        
    
    # Analysis - Differential expression
    diffExprTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
    
    # Analysis - Clustering
    clusteringTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
        
    # Analysis - Correlation
    correlationTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
        
    # Analysis - MDS/PCA
    mdsPcaTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
    
    # Genome browser
    genomeBrowserTabPanelObserve(input,output,session,allReactiveVars,
        allReactiveMsgs)
        
    onRestore(function(state) {
        #assign("state",state,envir=.GlobalEnv)
        
        shinyjs::show("spinnerContainer")
        
        query <- getQueryString()
        # If onRestore has fired, it means that query is not empty
        if (!is.null(query$code)) { # Fired with auth0, further check
            if (!is.null(query$`_state_id_`)) { # Then fire server script
                updateNavbarPage(session,"seqcnavbar",selected="Data selector")
                #showModal(modalDialog("Please wait while the session is being ",
                #    "restored to its bookmarked state. Remember To 'Clear ",
                #    "Dataset'if you want to start over!",
                #    title="Session Bookmark",
                #    easyClose=TRUE,
                #    footer=NULL,
                #    fade=TRUE
                #))
                
                # Reruning the whole server script to restore Bookmarked session
                shinyjs::delay(3000,{
                dataSelectorTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                geneSignalTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                areaSignalTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)  
                expressionExplorerTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                expressionCalculatorTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)   
                diffExprTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                clusteringTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                correlationTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                mdsPcaTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                genomeBrowserTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                #shinyjs::hide("spinnerContainer")
                #showNotification('Restoring session...',duration=10,
                #    type='message')
                })
            }
        }
        else { # Just check _state_id_
            if (!is.null(query$`_state_id_`)) { # Then fire server script
                updateNavbarPage(session,"seqcnavbar",selected="Data selector")
                #showModal(modalDialog("Please wait while the session is being ",
                #    "restored to its bookmarked state. Remember To 'Clear ",
                #    "Dataset'if you want to start over!",
                #    title="Session Bookmark",
                #    easyClose=TRUE,
                #    footer=NULL,
                #    fade=TRUE
                #))
                
                # Reruning the whole server script to restore Bookmarked session
                shinyjs::delay(3000,{
                dataSelectorTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                geneSignalTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                areaSignalTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)  
                expressionExplorerTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                expressionCalculatorTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)   
                diffExprTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                clusteringTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                correlationTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                mdsPcaTabPanelObserve(state$input,state$output,session,
                    allReactiveVars,allReactiveMsgs)
                genomeBrowserTabPanelObserve(state$input,output,session,
                    allReactiveVars,allReactiveMsgs)
                #shinyjs::hide("spinnerContainer")
                #showNotification('Restoring session...',duration=10,
                #    type='message')
                })
            }
        }
    })
    
    onRestored(function(state) {
        observe({
            # Will always be filled after a session restore, as session creation
            # is not allowed if a dataset has not been created
            if (!is.null(allReactiveVars$currentMetadata$final)) {
                shinyjs::hide("spinnerContainer")
                showModal(modalDialog(HTML("The selected session has been ",
                    "restored! Remember To hit <strong>Clear Dataset</strong> ",
                    "if you want to start over!"),
                    title="Session restored",
                    easyClose=TRUE,
                    size="s"
                ))
            }
        })
    })
    
    # Excluding unnecessary actions from bookmarking
    setBookmarkExclude(c("bookmarkBtn","showFastqc","clearDataset","deleteBM"))
        
    onStop(function() {
        if (!is(metadata,"SQLiteConnection"))
            dbDisconnect(metadata)
    })
}
#,info=auth0_info("config/_auth0.yml"))

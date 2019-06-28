# server.R

# Load required libraries
source("config/init_server_globals.R")

shinyServer(
    function(input,output,session) {
        # Load init packages script
        source("config/init_packages.R")
        # Load SeqCVIBE libs
        source("server/reactiveVars.R",local=TRUE)
        source("server/dataSelectorTab/dataSelectorItem.R",local=TRUE)
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
        
    }
)

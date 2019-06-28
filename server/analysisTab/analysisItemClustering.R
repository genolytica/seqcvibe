clusteringTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    currentHeatmap <- allReactiveVars$currentHeatmap
        
    performRnaClustering <- eventReactive(input$performRnaClustering,{
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        D <- loadedData[[s]][[d]]
        genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
        
        # Determine which genes to use
        switch(input$rnaClusterGeneList,
            select = {
                g <- input$selectClusteringGeneName
            },
            custom = {
                g <- input$rnaClusteringCustomList
                g <- strsplit(g,split="\n")[[1]]
                genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(g,names(genes))
                na <- which(is.na(m))
                if (length(na)>0)
                    m <- m[-na]
                g <- genes[m]
            },
            degenes = {
                deTable <- currentRnaDeTable$totalTable
                if (is.null(deTable))
                    return()
                g <- rownames(deTable)
            }
        )
        
        # Determine the measurements table
        switch(input$rnaClusteringMeasureRadio,
            counts = {
                tab <- D$norm[g,samples,drop=FALSE]
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    tab <- rbind(tab,A)
                }
            },
            rpkm = {
                tab <- round(edgeR::rpkm(
                    D$counts[g,samples,drop=FALSE],
                    gene.length=D$length[g],
                    lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(edgeR::rpkm(A,
                        gene.length=currentCustomRnaTables$lengths,
                        lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                    tab <- rbind(tab,A)
                }
            },
            rpgm = {
                tab <- round(D$norm[g,
                    samples,drop=FALSE]/D$length[g],
                        digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(A/currentCustomRnaTables$lengths,
                        digits=6)
                    tab <- rbind(tab,A)
                }
            }
        )
        
        # We now have the tab, transformations based on other selections
        switch(input$rnaClusteringScaleRadio,
            natural = {
                tab <- tab
            },
            log2 = {
                tab <- round(log2(tab+1),digits=6)
            }
        )
        
        switch(input$rnaClusteringVariablesRadio,
            replicates = {
                X <- tab
                if (!is.null(meta$alt_id))
                    colnames(X) <- as.character(meta$alt_id)
            },
            averages = {
                cc <- unique(as.character(meta$class))
                classList <- vector("list",length(cc))
                names(classList) <- cc
                for (cl in cc)
                    classList[[cl]] <- 
                        as.character(meta$sample_id[which(meta$class==cl)])
                X <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,"mean","rpkm"))
                colnames(X) <- names(classList)
            },
            fc = {
                cc <- unique(as.character(meta$class))
                classList <- vector("list",length(cc))
                names(classList) <- cc
                for (cl in cc)
                    classList[[cl]] <- 
                        as.character(meta$sample_id[which(meta$class==cl)])
                control <- input$rnaClusteringControl
                treatments <- setdiff(cc,control)
                contrast <- c(treatments,control)
                contrast <- paste(contrast,collapse="_vs_")
                X <- round(makeFoldChange(contrast,classList,tab,
                    input$rnaClusteringScaleRadio),6)
            }
        )
        
        rownames(X) <- names(genes[match(rownames(X),genes)])
        
        # Check validity of cluster numbers
        kr <- as.integer(input$kRowInput)
        if (is.na(kr) || kr<=0 || kr>nrow(X)) {
            output$rnaClusteringSettingsError <- renderUI({
                div(class="error-message",paste("The number of row clusters (",
                    kr,") cannot exceed the number of rows in the data matrix ",
                    "(",nrow(X),")!",sep=""))
            })
            return()
        }
        kc <- as.integer(input$kColInput)
        if (is.na(kc) || kc<=0 || kc>ncol(X)) {
            output$rnaClusteringSettingsError <- renderUI({
                div(class="error-message",paste("The number of column clusters",
                    " (",kc,") cannot exceed the number of columns in the ",
                    "data matrix (",ncol(X),")!",sep=""))
            })
            return()
        }
        
        # Check validity of color saturation quantiles
        if (input$checkColorSaturation) {
            ql <- as.numeric(input$colorSaturationLowQuantile)
            if (is.na(ql) || ql<0 || ql>1) {
                output$rnaClusteringSettingsError <- renderUI({
                    div(class="error-message",paste("The lower quantile for ",
                        "color saturation must be a number from 0 to 1.",
                        sep=""))
                })
                return()
            }
            qh <- as.numeric(input$colorSaturationLowQuantile)
            if (is.na(qh) || qh<0 || qh>1) {
                output$rnaClusteringSettingsError <- renderUI({
                    div(class="error-message",paste("The upper quantile for ",
                        "color saturation must be a number from 0 to 1.",
                        sep=""))
                })
                return()
            }
            theq <- quantile(X,c(ql,qh))
            X[X<theq[1]] <- theq[1]
            X[X>theq[2]] <- theq[2]   
        }
        
        switch(input$rnaClusterWhat,
            both = {
                rowv <- TRUE
                colv <- TRUE
            },
            row = {
                rowv <- TRUE
                colv <- FALSE
            },
            column = {
                rowv <- FALSE
                colv <- TRUE
            },
            none = {
                rowv <- FALSE
                colv <- FALSE
            }
        )
        currentHeatmap$data <- X
        currentHeatmap$opts$dendrogram <- input$rnaClusterWhat
        currentHeatmap$opts$distfun <- input$selectClusteringDistance
        currentHeatmap$opts$hclustfun <- input$selectClusteringLinkage
        currentHeatmap$opts$k_row <- as.integer(input$kRowInput)
        currentHeatmap$opts$k_col <- as.integer(input$kColInput)
        currentHeatmap$opts$Rowv <- rowv
        currentHeatmap$opts$Colv <- colv
        currentHeatmap$opts$colors <- c(
            input$heatmapColourExtremeDown,
            input$heatmapColourMildDown,
            input$heatmapColourMiddle,
            input$heatmapColourMildUp,
            input$heatmapColourExtremeUp
        )
    })
    
    return(list(
        performRnaClustering=performRnaClustering
    ))
}

clusteringTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

clusteringTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentHeatmap <- allReactiveVars$currentHeatmap
    
    output$rnaClusteringFcControl <- renderUI({
        if (is.null(currentMetadata$final))
            list(
                div(style="font-weight:bold","Select control condition"),
                helpText("Please create a dataset first from the 'Data ",
                    "selector' menu on the top")
            )
        else {
            cls <- unique(as.character(currentMetadata$final$class))
            if (length(cls==2))
                list(
                    div(style="font-weight:bold","Select control condition"),
                    div(class="small",
                    helpText("Fold changes cannot be used for clustering with ",
                        "only two conditions in the dataset as this would ",
                        "produce a 1-column matrix which cannot be clustered.")
                    )
                )
            else
                selectInput(
                    inputId="rnaClusteringControl",
                    label="Select control condition",
                    choices=cls
                )
        }
    })
    
    output$heatmapOutput <- renderUI({
        if (is.null(currentHeatmap$data)) {
            output$heatmap <- renderPlot({
                currentHeatmap$entry
            })
            plotOutput("heatmap",height="640px")
        }
        else {
            require(R.utils)
            distFun <- distFuns()
            hclustFun <- hclustFuns()
            tryCatch({
                evalWithTimeout({
                    output$heatmap <- renderD3heatmap({
                        d3heatmap(
                            currentHeatmap$data,
                            dendrogram=currentHeatmap$opts$dendrogram,
                            Rowv=currentHeatmap$opts$Rowv,
                            Colv=currentHeatmap$opts$Colv,
                            distfun=distFun[[currentHeatmap$opts$distfun]],
                            hclustfun=
                                hclustFun[[currentHeatmap$opts$hclustfun]],
                            k_row=currentHeatmap$opts$k_row,
                            k_col=currentHeatmap$opts$k_col,
                            colors=currentHeatmap$opts$colors
                        )
                    })
                    div(
                        class="heatmap-container",
                        d3heatmapOutput("heatmap",height="800px")
                    )
                },
                timeout=300)
            },
            TimeoutException=function(ex) {
                output$heatmap <- renderPlot({
                    currentHeatmap$timeout
                })
                plotOutput("heatmap",height="640px")
            },
            error=function(e) {
                output$heatmap <- renderPlot({
                    currentHeatmap$error
                })
                plotOutput("heatmap",height="640px")
            })
        }
    })
    
    output$exportRnaHeatmapPDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_heatmap_",tt,".pdf", sep='')
        },
        content=function(con) {
            require(gplots)
            distFun <- distFuns()
            hclustFun <- hclustFuns()
            pdf(file=con,width=6,height=10)
            heatmap.2(
                currentHeatmap$data,
                dendrogram=currentHeatmap$opts$dendrogram,
                Rowv=currentHeatmap$opts$Rowv,
                Colv=currentHeatmap$opts$Colv,
                distfun=distFun[[currentHeatmap$opts$distfun]],
                hclustfun=hclustFun[[currentHeatmap$opts$hclustfun]],
                col=colorRampPalette(currentHeatmap$opts$colors),
                trace="none"
            )
            dev.off()
        }
    )
    
    output$exportRnaHeatmapPNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_heatmap_",tt,".png", sep='')
        },
        content=function(con) {
            require(gplots)
            distFun <- distFuns()
            hclustFun <- hclustFuns()
            png(filename=con,width=480,height=800)
            heatmap.2(
                currentHeatmap$data,
                dendrogram=currentHeatmap$opts$dendrogram,
                Rowv=currentHeatmap$opts$Rowv,
                Colv=currentHeatmap$opts$Colv,
                distfun=distFun[[currentHeatmap$opts$distfun]],
                hclustfun=hclustFun[[currentHeatmap$opts$hclustfun]],
                col=colorRampPalette(currentHeatmap$opts$colors),
                trace="none"
            )
            dev.off()
        }
    )
}

clusteringTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    currentHeatmap <- allReactiveVars$currentHeatmap
    
    clusteringTabPanelReactiveEvents <- 
        clusteringTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    performRnaClustering <- 
        clusteringTabPanelReactiveEvents$performRnaClustering
    
    clusteringTabPanelReactiveExprs <- 
        clusteringTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    clusteringTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"selectClusteringGeneName",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectClusteringGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectClusteringGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (isEmpty(currentMetadata$final))
            shinyjs::disable("performRnaClustering")
        else {
            if (input$rnaClusterGeneList %in% c("select","custom") 
                && isEmpty(input$selectClusteringGeneName) 
                && isEmpty(input$rnaClusteringCustomList))
                shinyjs::disable("performRnaClustering")
            else if (input$rnaClusterGeneList=="degenes" 
                && is.null(currentRnaDeTable$totalTable))
                shinyjs::disable("performRnaClustering")
            else if (input$rnaClusteringVariablesRadio=="fc" 
                && length(unique(currentMetadata$final$class))<=2)
                shinyjs::disable("performRnaClustering")
            else {
                kr <- as.integer(input$kRowInput)
                kc <- as.integer(input$kColInput)
                if (is.na(kr) || kr<=0 || is.na(kc) || kc<=0) {
                    output$rnaClusteringSettingsError <- renderUI({
                        div(class="error-message",paste("The number of ",
                            "row and column clusters must be a positive ",
                            "integer!",sep=""))
                    })
                    shinyjs::disable("performRnaClustering")
                }
                else {
                    output$rnaClusteringSettingsError <- renderUI({div()})
                    shinyjs::enable("performRnaClustering")
                }
            }
        }
    })
    
    observe({
        tryCatch({
            shinyjs::disable("performRnaClustering")
            performRnaClustering()
        },error=function(e) {
            #print(e)
        },
        finally={
            if (isEmpty(currentMetadata$final) || is.null(currentHeatmap$data))
                shinyjs::disable("performRnaClustering")
            else {
                if (input$rnaClusterGeneList %in% c("select","custom") 
                    && isEmpty(input$selectClusteringGeneName) 
                    && isEmpty(input$rnaClusteringCustomList))
                    shinyjs::disable("performRnaClustering")
                else if (input$rnaClusterGeneList=="degenes" 
                    && is.null(currentRnaDeTable$totalTable))
                    shinyjs::disable("performRnaClustering")
                else if (input$rnaClusteringVariablesRadio=="fc" 
                    && length(unique(currentMetadata$final$class))<=2)
                    shinyjs::disable("performRnaClustering")
                else
                    shinyjs::enable("performRnaClustering")
            }
        })
    })
}

